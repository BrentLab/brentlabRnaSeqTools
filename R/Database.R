# TODO use inheritParams etc to reduce redundancy in documentation

# TODO reformulate the dplyr stuff as SQL -- can use show_query() along with
# dplyr interactively to help with constructing the query. All the pivoting
# will be hard in sql

#' Get the combined metadata as a tibble from a remote database
#'
#' @description Join the biosample, rnasample, s1sample, s2sample, library, fastqFiles and qualityAssessment tables (in that order, left joins) and return the result as a tibble
#'
#' @importFrom RPostgres dbDisconnect
#' @importFrom dplyr left_join tbl
#' @importFrom stringr str_remove
#' @importFrom tidyr replace_na
#'
#' @description Use the RPostgres package to connect to a remote postgresql database, do the table joining, and return the joined metadata as a tibble. The database connection is closed
#' @param database_host if connecting to a database hosted on AWS, it might be something like ec2-54-83-201-96.compute-1.amazonaws.com
#' @param database_name name of the database, eg for cryptococcus kn99, the database might be named kn99_database. Check with the documentation, whoever set up the database, or get into the server and check directly
#' @param database_user a user of the actual database, with some level of permissions. You'll need to check with the database maintainer for this. It is suggested that you use a .Renviron file in your local project (make sure it is completely ignored by git, R, etc) to store this info
#' @param database_password password to the database user. You'll need to check with the database maintainer for this. It is suggested that you use a .Renviron file in your local project (make sure it is completely ignored by git, R, etc) to store this info
#' @note for information on using R environmental files, see \url{https://support.rstudio.com/hc/en-us/articles/360047157094-Managing-R-with-Rprofile-Renviron-Rprofile-site-Renviron-site-rsession-conf-and-repos-conf}
#' @source \url{https://rpostgres.r-dbi.org/}
#' @return A DBI connection to the remote database
#'
#' @export
getMetadata = function(database_host, database_name, database_user, database_password){

  # connect to the db and pull the tables
  db = connectToDatabase(database_host, database_name, database_user, database_password)

  baseNutMix         = tbl(db, "baseNutrientMix")
  nutMixMod          = tbl(db, "NutrientMixMod")
  aa                 = tbl(db, "aminoAcids")
  gc                 = tbl(db, "growthConditions")
  biosample          = tbl(db, 'bioSample') %>% collect()
  rnasample          = tbl(db, 'rnaSample') %>% collect()
  s1sample           = tbl(db, 's1cDNASample') %>% collect()
  s2sample           = tbl(db, 's2cDNASample') %>% collect()
  library            = tbl(db, 'library') %>% collect()
  fastqFiles         = tbl(db, 'fastqFiles') %>% collect()
  quality            = tbl(db, 'qualityAssessment') %>% collect()
  replicateAgreement = tbl(db, 'replicateAgreement') %>% collect()

  combined_nutrient_df =
    # join the growth conditions nad nutrient tables
    gc %>%
    left_join(nutMixMod, by = c('nutrientMixModName_id' = 'nutrientMixModName')) %>%
    left_join(baseNutMix, by = c('baseNutrientMixName_id' = 'baseNutrientName')) %>%
    # pull from database into memory
    collect()  %>%
    # pivot the columns with suffixes .x and .y -- these are the duplicated
    # nutrient columns
    pivot_longer(cols = ends_with(c(".x", ".y"))) %>%
    # separate the eg nutrient.x into two columns, the nut column (nutrient) and source column (y)
    separate(name, into = c("nut", "source"), sep = "\\.") %>%
    # group by the growth condition and the nut
    group_by(gcid, nut) %>%
    # and add the nut values, eg adenine.x and adenine.y for gcid 8 are added together
    summarise(mod_value = sum(value), .groups = 'keep')

  conditions_df = gc %>%
    left_join(aa, by = c('aaid_id' = 'aaid')) %>%
    collect() %>%
    left_join(combined_nutrient_df) %>%
    pivot_wider(names_from = nut, values_from = mod_value)

  joined_meta_tables = biosample %>%
    left_join(conditions_df, by = c('gcid_id' = 'gcid')) %>%
    left_join(rnasample, by = c('bioSampleNumber' = 'bioSampleNumber_id'))%>%
    left_join(s1sample, by = c('rnaSampleNumber' = 'rnaSampleNumber_id'))%>%
    left_join(s2sample, by = c('s1cDNASampleNumber' = 's1cDNASampleNumber_id'))%>%
    left_join(library, by = c('s2cDNASampleNumber' = 's2cDNASampleNumber_id'))%>%
    left_join(fastqFiles, by = c('librarySampleNumber' = 'librarySampleNumber_id'))%>%
    left_join(quality, by = c('fastqFileNumber' = 'fastqFileNumber_id')) %>%
    left_join(replicateAgreement, by = c('fastqFileNumber' = 'fastqFileNumber_id'))

  metadata_df = as_tibble(joined_meta_tables) %>%
    # cast timePoint from integer64 to integer
    mutate(timePoint = as.integer(timePoint)) %>%
    # replace empty strings in the following columns with NA
    mutate(across(c("timePoint"),
                  ~ifelse(.=="", NA, .) )) %>%
    # replace NA with defined value
    mutate(timePoint = replace_na(as.character(timePoint), 'noTimepoint')) %>%
    mutate(fastqFileName = str_remove(fastqFileName, ".fastq.gz")) %>%
    mutate(libraryProtocol = as.factor(libraryProtocol)) %>%
    mutate(libraryDate = as.factor(libraryDate)) %>%
    mutate(genotype1 = as.factor(genotype1))

  dbDisconnect(db)

  # cast integer64 to integers. see Utils is_integer64
  metadata_df <- metadata_df %>%
    mutate_if(is_integer64, as.integer)

  return(metadata_df)
}

# TODO: storing the counts the 'wide' way in the database is WRONG. a long table
# is created -- needs to be filled and this can then be re-written

#' Get combined raw counts
#'
#' @importFrom RPostgres dbGetQuery dbDisconnect
#' @importFrom dplyr bind_cols
#' @importFrom jsonlite fromJSON
#' @importFrom purrr map
#'
#' @param database_host if connecting to a database hosted on AWS,
#'                      it might be something like ec2-54-83-201-96.compute-1.amazonaws.com
#' @param database_name name of the database, eg for cryptococcus kn99, the database might be named kn99_database.
#'                      Check with the documentation, whoever set up the database, or get into the server and check
#'                      directly
#' @param database_user a user of the actual database, with some level of permissions. You'll need to check with the
#'                      database maintainer for this. It is suggested that you use a .Renviron file in your
#'                      local project (make sure it is completely ignored by git, R, etc) to store this info
#' @param database_password password to the database user. You'll need to check with the database maintainer for this.
#'                          It is suggested that you use a .Renviron file in your local project
#'                          (make sure it is completely ignored by git, R, etc) to store this info
#'
#' @return a gene by samples dataframe of all counts
#'
#' @export
getRawCounts = function(database_host, database_name, database_user, database_password){
  db = connectToDatabase(database_host, database_name, database_user, database_password)

  counts = dbGetQuery(db, 'select "rawCounts" from counts')

  counts_df = map(counts$rawCounts,
      ~head(as.data.frame(fromJSON(.), check.names = FALSE), 6967)) %>%
    bind_cols()
  dbDisconnect(db)

  message("WARNING: genes have been subset down to the first 1:6967 gene indicies")
  counts_df
}

#' Get gene names
#'
#' @importFrom RPostgres dbGetQuery dbDisconnect
#' @importFrom dplyr bind_cols
#' @importFrom jsonlite fromJSON
#'
#' @param database_host if connecting to a database hosted on AWS,
#'                      it might be something like ec2-54-83-201-96.compute-1.amazonaws.com
#' @param database_name name of the database, eg for cryptococcus kn99, the database might be named kn99_database.
#'                      Check with the documentation, whoever set up the database, or get into the server and check
#'                      directly
#' @param database_user a user of the actual database, with some level of permissions. You'll need to check with the
#'                      database maintainer for this. It is suggested that you use a .Renviron file in your
#'                      local project (make sure it is completely ignored by git, R, etc) to store this info
#' @param database_password password to the database user. You'll need to check with the database maintainer for this.
#'                          It is suggested that you use a .Renviron file in your local project
#'                          (make sure it is completely ignored by git, R, etc) to store this info
#'
#' @return a dataframe of gene_ids
#'
#' @export
getGeneNames = function(database_host,
                        database_name,
                        database_user,
                        database_password){



  db = connectToDatabase(database_host, database_name, database_user, database_password)

  df = as_tibble(dbGetQuery(db, 'select * from "genes"'))

  dbDisconnect(db)

  return (df)
}


#' pull entire database (not counts) and save to output_dir for archival purposes
#'
#' @description saves both the individual tables, including counts, and the combined_df
#'
#' @importFrom RPostgres dbDisconnect
#' @importFrom readr write_csv
#' @importFrom dplyr tbl
#'
#' @param database_host if connecting to a database hosted on AWS, it might be something like ec2-54-83-201-96.compute-1.amazonaws.com
#' @param database_name name of the database, eg for cryptococcus kn99, the database might be named kn99_database. Check with the documentation, whoever set up the database, or get into the server and check directly
#' @param database_user a user of the actual database, with some level of permissions. You'll need to check with the database maintainer for this. It is suggested that you use a .Renviron file in your local project (make sure it is completely ignored by git, R, etc) to store this info
#' @param database_password password to the database user. You'll need to check with the database maintainer for this. It is suggested that you use a .Renviron file in your local project (make sure it is completely ignored by git, R, etc) to store this info
#' @param output_dir where to deposit a subdirectory, named by todays date in this format: 20210407, with the tables and combined_df inside. eg a mounted local directory /mnt/htcf_lts/crypto_database_archive/ --> /lts/mblab/Crypto/rnaseq_data/crypto_database_archive
#' @param archive_counts_flag boolean indicating whether or not to save the counts. default is TRUE
#' @return None, writes a directory called <today's date> with tables and combined_df as .csv to output_dir
#'
#' @export
archiveDatabase = function(database_host, database_name, database_user, database_password, output_dir, archive_counts_flag = TRUE){

  today_date = format(Sys.Date(), "%Y%m%d")
  current_output_path = file.path(output_dir, today_date)
  dir.create(current_output_path)

  db = connectToDatabase(database_host, database_name, database_user, database_password)

  tbl_list = list()

  tbl_list[['biosample']] = tbl(db, 'bioSample')
  tbl_list[['rnasample']] = tbl(db, 'rnaSample')
  tbl_list[['s1sample']] = tbl(db, 's1cDNASample')
  tbl_list[['s2sample']] = tbl(db, 's2cDNASample')
  tbl_list[['library']] = tbl(db, 'library')
  tbl_list[['fastqFiles']] = tbl(db, 'fastqFiles')
  tbl_list[['qualityAssessment']] = tbl(db, 'qualityAssessment')

  tbl_list_names = names(tbl_list)

  tbl_list = lapply(tbl_list, as_tibble)
  names(tbl_list) = tbl_list_names # maybe not necessary

  lapply(names(tbl_list), function(x) write_csv(tbl_list[[x]], file.path(current_output_path, paste0(x, ".csv") )))

  combined_df = getMetadata(database_host, database_name, database_user, database_password)
  write_csv(combined_df, file.path(current_output_path, paste0("combined_df_", today_date,'.csv')))

  if(archive_counts_flag){
    counts_df = getRawCounts(database_host, database_name, database_user, database_password)
    write_csv(counts_df, file.path(current_output_path, "counts.csv"))
  }

  dbDisconnect(db)

}

#'
#' Connect to a remote postgresql database
#'
#' @importFrom RPostgres Postgres dbConnect
#'
#' @description Use the RPostgres package to connect to a remote postgresql database
#' @param database_host if connecting to a database hosted on AWS, it might be something like ec2-54-83-201-96.compute-1.amazonaws.com
#' @param database_name name of the database, eg for cryptococcus kn99, the database might be named kn99_database. Check with the documentation, whoever set up the database, or get into the server and check directly
#' @param database_user a user of the actual database, with some level of permissions. You'll need to check with the database maintainer for this. It is suggested that you use a .Renviron file in your local project (make sure it is completely ignored by git, R, etc) to store this info
#' @param database_password password to the database user. You'll need to check with the database maintainer for this. It is suggested that you use a .Renviron file in your local project (make sure it is completely ignored by git, R, etc) to store this info
#' @note for information on using R environmental files, see \url{https://support.rstudio.com/hc/en-us/articles/360047157094-Managing-R-with-Rprofile-Renviron-Rprofile-site-Renviron-site-rsession-conf-and-repos-conf}
#' @source \url{https://rpostgres.r-dbi.org/}
#' @return A DBI connection to the remote database
#'
#' @export
connectToDatabase = function(database_host, database_name, database_user, database_password){

  dbConnect(RPostgres::Postgres(),dbname = database_name,
            host = database_host, # i.e. 'ec2-54-83-201-96.compute-1.amazonaws.com'
            port = 5432, # or any other port specified by your DBA
            user = database_user,
            password = database_password)

}

#'
#' list tables in databse
#'
#' @importFrom RPostgres dbGetQuery
#'
#' @param db a connection to the database
#'
#' @seealso \url{https://www.postgresqltutorial.com/postgresql-show-tables/}
#' @return all tables in database
#'
#' @export
listTables = function(db){
  dbGetQuery(db, "SELECT *
                  FROM pg_catalog.pg_tables
                  WHERE schemaname != 'pg_catalog' AND
                  schemaname != 'information_schema';")
}

#'
#' get (via a http POST request) your user authentication token from the database
#'
#' @importFrom httr http_status content POST
#'
#' @param url check the database_info variable. for configured organisms, you can find this under database_info$organism$token_auth
#' @param username a valid username for the database. If you don't have one, then you'll need to ask for one to be created
#' @param password password associated with your username
#'
#' @note do not save your auth token in a public repository. For example, you might put it in your .Renviron and then make sure
#'       that your .Renviron is in your .gitignore. Otherwise, save it outside of a github tracked directory or otherwise ensure
#'       that it will not be pushed up to github
#'
#' @return the auth token associated with the username and password
#'
#' @export
getUserAuthToken = function(url, username, password){

  # see package httr for help
  token_response = POST(url=url,
                        body=list(username=username,
                                  password=password),
                        encode='json')

  if( http_status(token_response)$category == "Success" ){
    message("You might want to put your token in your .Renviron. If you do, please make sure the .Renviron file is in your .gitignore")
    httr::content(token_response)$token
  } else{
    message("There was a problem getting your token:")
    message(http_status(token_response)$message)
  }
}

#'
#' post new fastq sheet to database
#'
#' @importFrom httr content_type add_headers POST
#' @importFrom jsonlite toJSON
#' @importFrom dplyr select mutate
#'
#' @param database_fastq_url eg database_info$kn99_urls$FastqFiles. See see \code{\link{database_info}}
#' @param auth_token see \code{\link{getUserAuthToken}}
#' @param new_fastq_path path to new fastq sheet
#'
#' @export
postFastqSheet = function(database_fastq_url, auth_token, new_fastq_path){

  # see utils
  fastq_df = brentlabRnaSeqTools::readInData(new_fastq_path)

  fastq_colnames = colnames(fastq_df)

  augment_fastq_df = fastq_df %>%
    {if("libraryDate" %in% fastq_colnames)
      select(.,-libraryDate)
      else .} %>%
    {if("libraryPreparer" %in% fastq_colnames)
      select(.,-libraryPreparer)
      else .} %>%
    {if(!"fastqObservations" %in% fastq_colnames)
      mutate(.,fastqObservations = "")
      else .} %>%
    {if(!"laneNumber" %in% fastq_colnames)
      mutate(.,laneNumber = "")
      else .} %>%
    mutate(volumePooled = round(as.numeric(volumePooled), 15)) %>%
    mutate(tapestationConc = round(as.numeric(tapestationConc), 4))

  post_body = jsonlite::toJSON(augment_fastq_df, auto_unbox = TRUE)

  POST(url = database_fastq_url,
       add_headers(Authorization = paste("token" , auth_token, sep=" ")),
       content_type("application/json"),
       body = post_body,
       encode = 'json')
}


# TODO add dry run option
# TODO make the counts long, key on fastqFileNumber and gene_id

#'
#' post counts to database
#'
#' @importFrom httr content_type add_headers POST
#' @importFrom jsonlite toJSON
#' @importFrom dplyr filter pull
#' @importFrom stringr str_remove
#'
#' @description using the package httr, post the raw count .csv, which is the
#'   compiled counts for a given run, to the database
#'
#' @param database_counts_url eg database_info$kn99_urls$Counts.
#'   See \code{\link{database_info}}
#' @param run_number the run number of this counts sheet -- this is important
#'   b/c fastqFileNames aren't necessarily unique outside of their runs
#' @param auth_token see \code{\link{getUserAuthToken}}
#' @param count_df the counts read in as a dataframe. Make sure there is no
#'   'feature' or 'gene_id' column. Colnames need to be the sample identifier
#' @param fastq_table a recent pull of the database fastq table
#' @param count_file_suffix the suffix appended to the fastqFileName in the
#'   count file column headings. default is "_read_count.tsv"
#'
#' @return a list of httr::response() objects
#'
#' @export
postCounts = function(database_counts_url, run_number, auth_token,
                      count_df, fastq_table,
                      count_file_suffix = "_read_count.tsv"){

  # fastqFileNames may not be unique outside of their run
  fastq_table = filter(fastq_table, runNumber == run_number)

  # ensure that count_df is a dataframe of some sort
  stopifnot(sum(class(count_df) == "data.frame") > 0)

  # remove the gene_ids
  count_df = count_df[colnames(count_df) != 'gene_id']
  # remove filename suffix from colnames, leaving just the fastqFileName behind
  colnames(count_df) = str_remove(colnames(count_df), count_file_suffix)

  # remove suffixes from the fastqfiles if they exist
  fastq_table$fastqFileName = str_remove(fastq_table$fastqFileName, ".fastq.gz")

  # halt if there are sample names in the count_df that are not in the database
  stopifnot(setdiff(colnames(count_df), fastq_table$fastqFileName) == 0)

  # create named list with structure list(fastqFileName = fastqFileNumber, ...)
  fastqFileNumber_lookup_list = pull(fastq_table, fastqFileNumber)
  names(fastqFileNumber_lookup_list) = pull(fastq_table, fastqFileName)
  # filter to just those in count_df
  fastqFileNumber_lookup_list = fastqFileNumber_lookup_list[names(fastqFileNumber_lookup_list) %in% colnames(count_df)]

  # check that we still have the same number of samples
  stopifnot(length(fastqFileNumber_lookup_list) == length(colnames(count_df)))

  # send each column to the count table of the database
  res_list = list()
  for (column in colnames(count_df)){
    counts = list(as.integer(pull(count_df, column)))
    names(counts) = column

    post_body = jsonlite::toJSON(list(fastqFileNumber = fastqFileNumber_lookup_list[[column]],
                                      rawCounts = counts), auto_unbox = TRUE)

    res = POST(url=database_counts_url,
               add_headers(Authorization = paste("token" , auth_token, sep=" ")),
               content_type("application/json"),
               body=post_body,
               encode='json')

    res_list[[column]] = res

  }
  res_list
}

# TODO add dry run option

#'
#' post new qc sheet to database
#'
#' @importFrom httr content_type add_headers POST
#' @importFrom jsonlite toJSON
#' @importFrom dplyr rename filter left_join select
#' @importFrom stringr str_remove
#'
#' @description using the package httr, post the new qc sheet to the database
#'
#' @note there can be problems with dependencies and the rename function. this is working for now,
#'       but see here for more info \url{https://statisticsglobe.com/r-error-cant-rename-columns-that-dont-exist}
#'
#' @param database_qc_url eg database_info$kn99_urls$QualityAssess. see \code{\link{database_info}}.
#' @param auth_token \code{\link{getUserAuthToken}}
#' @param run_number the run number of this qc sheet -- this is important b/c fastqFileNames aren't necessarily unique
#'                   outside of their runs
#' @param new_qc_path path to the new counts csv
#' @param fastq_table_path path to a recent pull of the database fastq table
#'
#' @return a list of httr::response() objects
#'
#' @export
postQcSheet = function(database_qc_url, auth_token,
                       run_number, new_qc_path, fastq_table_path) {

  fastq_df = brentlabRnaSeqTools::readInData(fastq_table_path)

  # fastqFileNames may not be unique outside of their run
  fastq_df = filter(fastq_df, runNumber == run_number)

  new_qc_df = brentlabRnaSeqTools::readInData(new_qc_path)

  new_qc_df = new_qc_df %>%
    left_join(fastq_df %>%
                filter(runNumber == run_number) %>%
                select(fastqFileNumber))

  # check that we still have the same number of samples
  stopifnot(sum(is.na(new_qc_df$fastqFileNumber))==0)

  post_body = jsonlite::toJSON(new_qc_df, auto_unbox = TRUE)

  POST(url = database_qc_url,
       add_headers(Authorization = paste("token" , auth_token, sep=" ")),
       content_type("application/json"),
       body = post_body,
       encode = 'json')
}

# TODO add dry run option

#'
#' PATCH entries in database table
#'
#' @importFrom httr content_type add_headers PATCH
#' @importFrom jsonlite toJSON
#' @importFrom dplyr select
#'
#' @description using the package httr, update entries in certain fields in given rows of a table
#'
#' @param database_table_url NO TRAILING '/'. eg "http://18.224.181.136/api/v1/QualityAssess"
#' @param auth_token see brentlabRnaSeqTools::getUserAuthToken()
#' @param update_df a dataframe, preferrably a tibble, already read in, subsetted. Columns must be correct data type for db table
#' @param id_col name of the id column of the table. this number will be appended to the url to create the uri for the record
#'
#' @return a list of httr::response() objects
#'
#' @export
patchTable = function(database_table_url, auth_token, update_df, id_col){

  # send each column to the count table of the database
  res_list = list()
  for (i in seq(1,nrow(update_df))){

    id = update_df[[i,id_col]]
    row_as_list = jsonlite::toJSON(as.list(dplyr::select(update_df[i,], -id_col)), auto_unbox = TRUE, pretty=TRUE)

    base_url = str_remove(database_table_url, "/$")

    url = paste(base_url, paste0(as.character(id), "/"), sep="/")

    res = PATCH(url=url,
                add_headers(Authorization = paste("token" , auth_token, sep=" ")),
                content_type("application/json"),
                body=row_as_list,
                encode="json")

    res_list[[as.character(id)]] = res

  }
  res_list
}

#'
#' Post a table to the database
#'
#' @importFrom httr content_type add_headers POST
#' @importFrom jsonlite toJSON
#'
#' @param database_table_url see \code{\link{database_info}}. Use one of the URLS in the url slot
#' @param auth_token see \code{\link{getUserAuthToken}}
#' @param df a dataframe read in with, for example read_csv or vroom
#'
#' @return POST results object
#'
#' @export
postTable = function(database_table_url, auth_token, df){

  # TODO error check the df -- is it a data frame? does it have more than 0 rows?

  post_body = jsonlite::toJSON(df, auto_unbox = TRUE)

  POST(url = database_table_url,
       add_headers(Authorization = paste("token" , auth_token, sep=" ")),
       content_type("application/json"),
       body = post_body,
       encode = 'json')

}

#'
#' Send Counts to Database (legacy pipeline)
#'
#' @description for use in the novoalign + htseq brentlab rnaseq pipeline.
#'
#' @importFrom stringr str_remove
#' @importFrom dplyr select
#'
#' @param htseq_path path to htseq-counts output
#' @param run_number run number from which the sample originates
#' @param auth_token your user authentication token to the database. See
#'   \code{\link{getUserAuthToken}}
#' @param fastq_table_path path to recent fastq table reflecting current state of
#'   database
#'
#' @return http response
#'
#' @export
sendCountsToDatabase = function(htseq_path, run_number, auth_token,
                                fastq_table_path){

  if(!file.exists(fastq_table_path)){
    stop(paste0("'fastq_table_path': ", fastq_table_path, " does not exist."))
  } else if(tools::file_ext(fastq_table_path) != 'csv'){
    stop(paste0("'fastq_table_path': ", fastq_table_path, " must have extension
                .csv (and be a csv file)"))
  }

  fastq_table = read_csv(fastq_table_path)

  sample_name = str_remove(basename(htseq_path), "_read_count.tsv")

  htseq_out = readHTSeqFile(htseq_path, sample_name) %>%
    select(-feature)

  postCounts(database_info$kn99$urls$counts,
             run_number,
             auth_token,
             htseq_out,
             fastq_table)
}
