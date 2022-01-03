#' test if argument is numeric
#'
#' @description copied directly from the limma codebase
#' @param x any R object
#' @details copied from the limma docs: This function is used to check the validity of arguments for numeric functions. It is an attempt to emulate the behavior of internal generic math functions. IsNumeric differs from is.numeric in that data.frames with all columns numeric are accepted as numeric.
#'
#' @export
isNumeric = function (x) {
  is.numeric(x) || (is.data.frame(x) && length(x) > 0 && all(unlist(lapply(x, is.numeric))))
}

#'
#' correct run number to add leading zero where approprirate
#' @param run_number a run number, most likely from the metadata runNumber field
#'
#' @export
getRunNumberLeadingZero = function(run_number){
  run_number = as.character(run_number)
  if(!is.null(run_numbers_with_leading_zero[[run_number]])){
    brentlabRnaSeqTools::run_numbers_with_leading_zero[[run_number]]
  } else{
    run_number
  }
}

#'
#' test bam path
#' @param metadata_df path to metadata dataframe
#'
#' @export
testBamPath = function(metadata_df){
  for (i in seq(1,nrow(metadata_df))){
    run_number = getRunNumberLeadingZero(metadata_df[[1,'runNumber']])
    run_dir = paste0("run_", run_number, "_samples")

    sample_name = str_remove(metadata_df[[1,'fastqFileName']], ".fastq.gz")
    bam_name = paste0(sample_name, Sys.getenv("BAM_SUFFIX"))

    bam_path = file.path(Sys.getenv("LTS_ALIGN_EXPR_PREFIX"), run_dir, "align", bam_name)
    message(paste0("checking: ", bam_path))
    stopifnot(file.exists(bam_path))
  }
}

#'
#' read in columnar data
#'
#' @importFrom readr read_csv read_tsv
#' @importFrom readxl read_excel
#'
#' @description given a csv, tsv or excel sheet, use the right function to read in the data
#'
#' @param path path to a csv, tsv or xlsx
#'
#' @export
readInData = function(path){

  switch(tools::file_ext(path),
         "csv" = read_csv(path),
         "tsv" = read_tsv(path),
         "xlsx" = readxl::read_excel(path),
         message("File extension not recognized. Must be one of {csv, tsv, xlsx}"))
}

#'
#' Test if value is datatype integer64
#'
#' @description integer64 can cause some problems in filtering functions, eg
#'              if timePoint is an integer64 datatype column, then
#'              timePoint %in% c(0L, 1440L) and
#'              timePoint == 0 | timePoint == 1440 return different subsets
#'
#' @references \url{https://community.rstudio.com/t/how-to-handle-the-integer64-type/50024}
#'
#' @param x a numeric value
#'
#' @return boolean with truth value determined by whether input is integer64
#'
#' @export
is_integer64 <- function(x){
  class(x)=="integer64"
}

#'
#' parse text comparative sentence
#' @description return logical for sentence such as "3 > 4")
#'
#' @param value1 a number represented as a string, eg "3"
#' @param comparative one of c(">", ">=", "<", "<=", "==")
#' @param value2 a number represented as a string, eg "4"
#'
#' @return a logical as if the characters were evaluated such as 3 > 4
#'
#' @export
parseComparatives = function(value1, comparative, value2){
  switch (comparative,
          ">"  = value1 > value2,
          ">=" = value1 >= value2,
          "<"  = value1 < value2,
          "<=" = value1 <= value2,
          "==" = value1 == value2,
          stop("comparative not recognized")
  )
}

#' To view a data.frame / data.table in LibreOffice Calc
#'
#' @importFrom openxlsx write.xlsx
#'
#' @description Copied from
#' \url{https://gitlab.com/zauster/ormisc/-/blob/master/R/view.R}.
#' The purpose of this function is to write the dataframe to at tmp file
#' and print the path. Use the clipr package to write the output to the
#' clipboard for easy pasting into a terminal. See example.
#'
#' @references \url{https://gitlab.com/zauster/ormisc/-/blob/master/R/view.R}
#'
#' @param df the data.frame (or data.table or tibble)
#'
#' @examples
#' \dontrun{
#' install.packages(clipr)
#' library(clipr)
#' # df is some dataframe in your environment
#' write_clip(localView(df))
#' # in a terminal, you can how paste in the line that the cmd above put
#' # in your clipboard and hit enter. your computer will know what to do next
#' # hopefully
#' }
#'
#' @export
localView <- function(df) {
  open_command <- switch(Sys.info()[['sysname']],
                         Windows= 'open',
                         Linux  = 'xdg-open',
                         Darwin = 'open')

  temp_file <- paste0(tempfile(), '.xlsx')

  # data <- as.data.frame(df)
  write.xlsx(df, file = temp_file)
  open_cmd = paste(open_command, temp_file)
  open_cmd

  # this following does not work on my computer, but would automate this if it did.
  # this error:
  # Warning: failed to read path from javaldx
  # /usr/lib/libreoffice/program/soffice.bin: error while loading shared libraries:
  # libreglo.so: cannot open shared object file: No such file or directory
  # possible solution:
  # https://askubuntu.com/a/977080
  # invisible(system(open_cmd))
}
