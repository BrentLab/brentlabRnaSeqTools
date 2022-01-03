#' Parse HTSeq output to tibble
#'
#' @description parse the htseq output as a tibble with columns (feature, sample).
#'   This code is copied from the National Cancer Institute GenomicDataCommons
#'
#' @references
#'   \href{https://bioconductor.org/packages/release/bioc/html/GenomicDataCommons.html}{GenomicDataCommons::readHTSeqFile}
#'
#' @param fname filepath to htseq output
#' @param samplename in the legacy pipeline, the fastqFileName with .fastq.gz
#'   replaced with _read_count.tsv
#'
#' @importFrom readr read_tsv
#'
#' @return a dataframe with columns feature, 'sample'
#'
#'@export
readHTSeqFile = function (fname, samplename = "sample"){

  if (!file.exists(fname)){
    stop(sprintf("The specified file, %s, does not exist", fname))
  }

  if (!((length(fname) == 1) & (is.character(fname)))){
    stop("fname must be of type character(1)")
  }

  tmp = read_tsv(fname, col_names = FALSE)

  if (ncol(tmp) != 2){
    stop(sprintf("%s had %d columns, expected 2 columns",
                 fname, ncol(tmp)))
  }

  colnames(tmp) = c("feature", samplename)

  tmp
}
