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

  message("WARNING: subsetting gene ids down to only those with prefix CKF44")
  gene_ids = gene_ids %>%
    filter(startsWith(gene_id, "CKF44")) %>%
    pull(gene_id)

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

# Variant Explorer -------------------------------------------------------------

#' @title VariantExplorer
#' @description Use this to coordinate the analysis of the snpEff summary,
#' alignment files and bam files.
#'
#' @rdname VariantExplorer
#' @export
setClass("VariantExplorer",
         slots = list(
           metadata = "data.frame",
           gff = "GenomicRanges",
           bsgenome = "BSgenome",
           igv_genome = "character",
           variants = "data.frame",
           expected_metadata_fields = "character"
         ),
         prototype = prototype(
           expected_metadata_fields = c("sample_id","bam","snpeff","variant_caller")
         ),
         contains = c("data.frame"))

# variant explorer validity check
setValidity2("VariantExplorer", function(self){
  msg = NULL
  if ( length(setdiff(self@expected_metadata_fields, colnames(self@metadata))) > 0){
    msg = sprintf("Colnames for metadata frame must include %s",
                  paste0(self@expected_metadata_fields, col=","))
  } else if(is.null(self@gff)){
    msg = paste0("The gff slot must be filled with a GenomicRanges object. ",
                 "Use rtracklayer::import to read in the gff")
  } else if(is.null(self@bsgenome)){
    msg = "You must provide a BSgenome object"
  } else if(identical(self@igv_genome, character(0))){
    warning("WARNING: no igv genome set!")
  }

  if (!is.null(msg)){
    msg
  } else{
    TRUE
  }
})


#' @rdname VariantExplorer
#' @param sample_metadata sample metadata with at least the columns sample_id,
#'   snpeff, variant_caller. Suggested columns are bam and vcf.
#' @param gff_ranges the gff for this organism imported with rtracklayer::import()
#' @param bsgenome path to a BSgenome object. KN99 has one -- check bioconductor
#' @param igv_genome path to a genome object created by igvtools
#'
#' @return A VariantExplorer object
#'
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#'
#' @seealso
#'  \code{\link[Biostrings]{character(0)}}
#'  \code{\link[GenomicRanges]{character(0)}}
#'
#' @export
create_variant_explorer = function(sample_metadata,
                                   gff_ranges,
                                   bsgenome,
                                   igv_genome = ""){

  required_cols = c('sample_id', "snpeff", "variant_caller")
  if (length(setdiff(required_cols, colnames(sample_metadata))) > 0){
    stop(sprintf("The following columns must be in  the colnames of sample_metadata: %s",
                 paste(required_cols, collapse=",")))
  }



  variant_frame = apply(sample_metadata,
                        1,
                        function(x) snpeff_summary_to_long(x[['snpeff']],
                                                           x[['sample_id']],
                                                           x[['variant_caller']])) %>%
    do.call('rbind',.)

  # todo: check that these fields exists
  gff_ranges[is.na(gff_ranges$gene)]$gene = gff_ranges[is.na(gff_ranges$gene)]$ID
  gff = gff_ranges[gff_ranges$gene %in% unique(variant_frame[['GeneId']])]


  stopifnot(length(unique(gff$gene)) == length(unique(variant_frame[['GeneId']])))

  new("VariantExplorer",
      metadata = sample_metadata,
      gff = gff,
      bsgenome = bsgenome,
      variants = variant_frame)

  }






















