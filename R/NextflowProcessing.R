#'
#' create nf-co sample sheet
#'
#' @description create a samplesheet for the nf-co rnaseq pipeline
#'              \url{https://nf-co.re/rnaseq/3.2/usage#samplesheet-input}
#'
#'
#' @note: currently hard coded to set strandedness to 'reverse' if libraryProtocol is E7420L, or 'unstranded',
#'        AND for single strand libraries (fastq2 is an empty string)
#'
#'
#' @importFrom tidyr unite
#' @importFrom dplyr mutate select
#'
#' @param metadata most likely from the database. must include runNumber, fastqFileName, libraryProtocol,
#'                 and the columns used to create the sample name (see param sample_columns)
#' @param sample_columns a character vector of columns to concat into a sample identifier.
#'                       eg c("genotype1", "fastqFileNumber")
#' @param sequence_dir_prefix path to directory where the sequence run directories are stored.
#'                            eg /mnt/htcf_lts/lts_sequence
#' @param check_files_flag default FALSE. Set to true to check that files exist
#'
#' @return a dataframe with the columns which correspond to those required for the nf-co rnaseq pipeline
#'         sample sheet. See \url{https://nf-co.re/rnaseq/3.2/usage#samplesheet-input}
#'
#' @export
createNfCorePipelineSampleSheet = function(metadata, sample_columns,
                                           sequence_dir_prefix,
                                           check_files_flag=FALSE){

  # remove trailing / if it exists
  sequence_dir_prefix = gsub("/$", "", sequence_dir_prefix)

  sample_sheet_df = metadata %>%
    mutate(strandedness = ifelse(libraryProtocol == "E7420L",
                                 "reverse",
                                 "unstranded")) %>%
    mutate(runNumber = unlist(lapply(runNumber,
                                     getRunNumberLeadingZero))) %>%
    mutate(fastq_1 = file.path(sequence_dir_prefix,
                               paste0("run_", runNumber, "_samples"),
                               fastqFileName)) %>%
    mutate(fastq_2 = "") %>%
    unite(sample, sample_columns, sep="_") %>%
    dplyr::select(sample, fastq_1, fastq_2, strandedness)

  unique_sample_check_flag =
    length(unique(sample_sheet_df$sample)) == nrow(sample_sheet_df)

  if(!unique_sample_check_flag){
    message(
      "The columns identified to create the sample column entries
      do not result in unique identifiers for the sample.
      The effect is that two samples with the same name will be collapsed.
      If this is not what you want to happen, then run this command again
      and include a column which will make the samples unique.
      Adding fastqFileNumber will always work."
    )
  }

  if(check_files_flag){
    message("Checking whether fastq_1 files exist...")
    missing_files = unlist(lapply(sample_sheet_df$fastq_1, file.exists))

    # if there are any FALSE, return a summary
    if(sum(!missing_files)>0){
      message("Some files could not be located. See output")
      return(sample_sheet_df[['fastq_1']][!missing_files])
    }
  }

  sample_sheet_df

}

#'
#' move fastq files from lts to scratch for nf-co/rnaseq pipeline
#' @description using the samplesheet created by createNfCorePipelineSampleSheet,
#' create a file to move fastq files from lts to scratch
#'
#' @importFrom dplyr mutate select bind_rows
#'
#' @param nf_co_samplesheet a samplesheet created by
#' createNfCorePipelineSamplesheet() with the prefixes set to what the path will
#' be once the files are moved to scratch
#' @param from_prefix filepath from source up to the filename of the file
#'
#' @return a dataframe with two columns, source and destination
#'
#' @export
moveNfCoFastqFiles = function(nf_co_samplesheet, from_prefix){

  if(unique(nf_co_samplesheet$fastq_2) == ''){

    fastq_1 = nf_co_samplesheet %>%
      mutate(source = file.path(from_prefix,
                                basename(dirname(fastq_1)),
                                basename(fastq_1))) %>%
      mutate(destination = fastq_1) %>%
      dplyr::select(source, destination)

    return(fastq_1)

  } else{

    fastq_2 = nf_co_samplesheet %>%
      mutate(source = file.path(from_prefix,
                                basename(dirname(fastq_1)),
                                basename(fastq_1))) %>%
      mutate(destination = fastq_2) %>%
      dplyr::select(source, destination)

    return(bind_rows(fastq_1, fastq_2))

  }


}
