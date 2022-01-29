#' @rdname brentlabRnaSeqSet
#' @export
setClass("brentlabRnaSeqSet",
         contains = "DESeqDataSet")


#' brentlabRnaSeqSet object and constructor
#'
#' @description A container for expression data in brentlab rnaseq databases
#'
#' @rdname brentlabRnaSeqSet
#'
#' @docType class
#'
#' @seealso
#' * [extending DESeqDataSet](https://github.com/mikelove/DESeq2)
#' * [extending SummarizedExperiment](https://www.bioconductor.org/packages/devel/bioc/vignettes/SummarizedExperiment/inst/doc/Extensions.html).
#'    Note that A [SummarizedExperiment] object becomes a [RangedSummarizedExperiment]
#'    object when rowRanges are added.
#'
#' @aliases brentlabRnaSeqSet brentlabRnaSeqSet-class brentlabRnaSeqSetFromDatabase
#'
#' @import methods
#' @importClassesFrom DESeq2 DESeqDataSet
#'
#' @param dds a deseq data object. Passing a null object will instantiate with
#'   the [DESeq2::makeExampleDESeqDataSet]
#' @param organism currently, only kn99
#' @param username db username
#' @param password db password
#' @param full_rnaseq_only filter out spikeins before returning. Default TRUE
#'
#' @return a brentlabRnaSeqSet for the given organism
#'
#' @export
brentlabRnaSeqSet = function(dds = NULL) {

  if(is.null(dds)){
    dds = DESeq2::makeExampleDESeqDataSet()
  }

  new("brentlabRnaSeqSet", dds)
}

#' @rdname brentlabRnaSeqSet
#' @export
brentlabRnaSeqSetFromDatabase =  function(organism,
                                          username,
                                          password,
                                          full_rnaseq_only = TRUE){

  metadata = getMetadata(database_info[[organism]]$db_host,
                         database_info[[organism]]$db_name,
                         username,
                         password)

  # filter out any NA fastqFileNames
  metadata = metadata %>%
    filter(!is.na(fastqFileName)) %>%
    mutate(fastqFileName = str_remove(fastqFileName, ".fastq.gz"))

  raw_counts = getRawCounts(database_info[[organism]]$db_host,
                            database_info[[organism]]$db_name,
                            username,
                            password)

  gene_ids = getGeneNames(database_info[[organism]]$db_host,
                          database_info[[organism]]$db_name,
                          username,
                          password)

  gene_ids = gene_ids$gene_id

  if(full_rnaseq_only){
    # note -- this run was missing during lts re-org.
    metadata = filter(metadata, purpose == "fullRNASeq", runNumber != 5361)
  }

  if(length(setdiff(metadata$fastqFileName, colnames(raw_counts))) != 0){
    message(
      paste0("WARNING: there are fullRun fastqFileNames in the metadata ",
             "which are not present in the counts. This means that there is ",
             "an unprocessed run(s?). Proceeding with only those samples for ",
             "which counts exist.")
      )
    metadata = metadata %>%
      filter(fastqFileName %in% colnames(raw_counts))
  }

  raw_counts = raw_counts[, metadata$fastqFileName]

  rownames(raw_counts) = gene_ids

  # check that columns are in same order as rows of colData
  stopifnot(identical(colnames(raw_counts), metadata$fastqFileName))

  dds = DESeqDataSetFromMatrix(colData = metadata,
                               countData = raw_counts,
                               design = ~1)
  # for some reason, setting rownames on the count data didn't set the rownames
  # of the object. For SummarizedExperiment, it should be the case that it does.
  # maybe something to do with DESeqDataSetFromMatrix
  rownames(dds) = gene_ids

  brentlabRnaSeqSet(dds)

}


#' @rdname brentlabRnaSeqExperiment
#' @export
setClass("brentlabRnaSeqExperiment",
         contains="brentlabRnaSeqSet",
         slots = representation(
           set_name = "character"
         ))

setValidity2("brentlabRnaSeqExperiment", function(x) {
  msg = NULL

  # check that the set name is recognized
  if (!x@set_name %in%
        c('ninetyMin_2016Grant', 'ninetyMin_2016GrantWithDoubles',
          'ninetyMin_non2016grant', 'ninetyMin_all', 'envPert_epWT',
          'envPert_perturbed', 'envPert_titrationWT')) {
    msg = c(msg, "'set_name', the name of the experiment set,
            must be one of {'envPert', 'ninetyMin'}")
  }

  if (is.null(msg)) {
    TRUE
  } else msg
})

#' brentlabRnaSeqExperiment constructor.
#'
#' @description Extends brentlabRnaSeqSet to use with pre-specified experiment
#'   sets, eg Environmental Perturbation. Intended for lower level use, primarily.
#'
#' @rdname brentlabRnaSeqExperiment
#'
#' @docType class
#' @aliases brentlabRnaSeqExperiment brentlabRnaSeqExperiment-class
#'
#' @param brentlabRnaSeqSet A brentlabRnaSeqSet
#' @param set_name one of {'ninetyMin_2016Grant', 'ninetyMin_2016GrantWithDoubles',
#'  'ninetyMin_non2016grant', 'ninetyMin_all',  'envPert_epWT',
#'  'envPert_perturbed', 'envPert_titrationWT'}
#'
#' @return a brentlabRnaSeqExperiment object
#'
#' @export
brentlabRnaSeqExperiment <- function(brentlabRnaSeqSet,
                                                set_name = "") {
  blrs <- brentlabRnaSeqSet
  if (!is(blrs, "brentlabRnaSeqSet")) {
      stop("'brentlabRnaSeqSet' must be a brentlabRnaSeqSet object")
  }

  new("brentlabRnaSeqExperiment", blrs, set_name = set_name)
}

#' @rdname brentlabRnaSeqSetTransform
#' @export
setClass("brentlabRnaSeqSetTransform", contains="DESeqTransform")

#' brentlabRnaSeqSetTransform object and constructor. This extends DESeqTransform.
#'
#' Extends DESeqTransform for use in rleTransform
#' It is used by \code{\link{rleTransform}}
#'
#' @param SummarizedExperiment a RangedSummarizedExperiment
#'
#' @return a brentlabRnaSeqSetTransform object
#' @docType class
#' @aliases brentlabRnaSeqSetTransform-class
#' @rdname brentlabRnaSeqSetTransform
#' @export
brentlabRnaSeqSetTransform <- function(SummarizedExperiment) {
  se <- SummarizedExperiment
  if (!is(se, "RangedSummarizedExperiment")) {
    if (is(se, "SummarizedExperiment")) {
      se <- as(se, "RangedSummarizedExperiment")
    } else {
      stop("'SummarizedExperiment' must be a RangedSummarizedExperiment object")
    }
  }
  new("brentlabRnaSeqSetTransform", se)
}
