#' Parse HTSeq output to tibble
#'
#' @description parse the htseq output as a tibble with columns (feature, sample).
#'   This code is copied from the National Cancer Institute GenomicDataCommons
#'
#' @references
#'   \href{https://bioconductor.org/packages/release/bioc/html/GenomicDataCommons.html}{GenomicDataCommons::readHTSeqFile}
#'
#' @param htseq_filename filepath to htseq output
#' @param samplename in the legacy pipeline, the fastqFileName with .fastq.gz
#'   replaced with _read_count.tsv
#'
#' @importFrom readr read_tsv
#'
#' @return a dataframe with columns feature, 'sample'
#'
#'@export
readHTSeqFile = function (htseq_filename, samplename = "sample"){

  if (!file.exists(htseq_filename)){
    stop(sprintf("The specified file, %s, does not exist", htseq_filename))
  }

  if (!((length(htseq_filename) == 1) & (is.character(htseq_filename)))){
    stop("htseq_filename must be of type character(1)")
  }

  tmp = read_tsv(htseq_filename, col_names = FALSE, show_col_types = FALSE)

  if (ncol(tmp) != 2){
    stop(sprintf("%s had %d columns, expected 2 columns",
                 htseq_filename, ncol(tmp)))
  }

  colnames(tmp) = c("feature", samplename)

  tmp
}

#' Parse a novoalign log
#'
#' @importFrom dplyr filter mutate
#' @importFrom readr read_csv
#' @importFrom stringr str_detect str_remove str_remove_all
#' @importFrom tidyr separate
#'
#' @param novoalign_log_path path to a given sample's novoalign log
#'
#' @return a dataframe of with columns metric, reads, percent. Metrics will
#'   include: Read Sequences, Unique Alignment, Multi Mapped, No Mapping Found,
#'   Homopolymer Filter
#'
#' @export
parseNovoalignLog = function(novoalign_log_path){

  if(!file.exists(novoalign_log_path)){
    stop(paste0("path to novoalign log DNE: ", novoalign_log_path))
  }

  novoalign_log_metrics =
    paste0("Read Sequences|Unique Alignment|Multi Mapped|",
           "No Mapping Found|Homopolymer Filter")

  # TODO: re-direct warnings to a log file?
  suppressWarnings(
    read_csv(novoalign_log_path, col_names = "metric", show_col_types = FALSE) %>%
    filter(str_detect(metric, novoalign_log_metrics)) %>%
    separate(metric, into = c("metric", "reads"), sep = ":") %>%
    mutate(metric = str_remove(metric, "# +"),
           reads = trimws(reads)) %>%
    separate(reads, into = c('reads', 'percent'), sep = " ", extra = "merge") %>%
    mutate(percent = trimws(str_remove_all(percent, "\\(|%\\)"))) %>%
    mutate(reads = as.integer(reads),
           percent = as.double(percent)))


}

#' Given a GenomicFeatures annotation_db and a gene_id, extract an GRanges
#'   object of the cds
#'
#' @importFrom stringr str_detect regex
#'
#' @param annote_obj_path path to an annotation file parsed by rtracklayer::import
#' @param gene_id the ID of a gene in the db. Eg, for cryptococcus CKF44_05222
#' @param id_col gene feature column. Default is 'ID'
#' @param feature_col feature (col 3) column of annote_obj. Default is 'type'
#' @param feature what feature to select. Default is 'cds'
#' @references rtracklayer::import
#'
#' @return an IRanges object of a given feature (eg, a gene's cds features)
#'
#' @export
geneGRanges = function(annote_obj_path, gene_id, id_col = "ID",
                       feature_col = "type", feature = "cds"){

  annot_obj = readRDS(annote_obj_path)

  if(!class(annot_obj) == "GRanges"){
    stop(paste0("annote_obj_path must lead to a GRanges object. ",
                "It should be created with a command similar to the ",
                "following: rtracklayer::import('path/to/gff')"))
  }

  regions = tryCatch(
    expr = {
      annot_obj[grepl(gene_id, annot_obj@elementMetadata[,id_col]) &
                  grepl(feature, annot_obj@elementMetadata[,feature_col],
                        ignore.case = TRUE),]
    },
    error = function(e){
      message(
        paste0('geneGRanges() Error: cannot create GRanges for gene_id: ',
               gene_id))
      print(e)
    },
    warning = function(w){
      message("geneGRanges() warning: ")
      print(w)
    },
    finally = {
      # none
    }
  )

  return(regions)
}

#'
#' A helper function to create a ScanBamParam object
#'
#' @description helper function to create ScanBamParam object with appropriate
#'   strandedness information
#'
#' @importFrom Rsamtools ScanBamParam scanBamFlag
#'
#' @param locus_granges a granges object for a given gene
#'   (or some other feature on only one strand)
#' @param library_strandedness one of c("reverse", "same", "unstranded")
#' @param quality_threshold quality threshold above which reads will be
#'   considered. 20l is default, which is chosen b/c it is the default for HTSeq
#'
#' @return a ScanBamParam object with certain configured options, as well as
#' some reasonable defaults, to filter a bam file for reads of interest based
#' on the strandedness (protocol) of the library prep.
#'
#' @seealso \code{\link[Rsamtools]{ScanBamParam}}
#'
#' @export
strandedScanBamParam = function(locus_granges, library_strandedness, quality_threshold=20L){

  # ensure the locus is entirely on the same strand, error out if not
  if(!length(gene_strand) == 1){
    message(paste0("The granges object contains features on both strands. ",
                   "Separate the granges object into sets so that any set is on ",
                   "the same strand and try again"))
  }

  # TODO add support for forward stranded libraries
  # set some information for the ScanBamParam object below. gene_strand extracts
  # the +/- strand from the GRanges object
  # and the conditional sets the minus_strand_flag used in the ScanBamParam
  # constructor. This determines what reads are
  # returned -- either from a given strand, or from both
  gene_strand = as.character(unique(data.frame(locus_granges)$strand))

  minus_strand_flag = switch (paste(strandedness, gene_strand, sep="_"),
                              "reverse_+" = TRUE,
                              "reverse_-" = FALSE,
                              "same_+" = FALSE,
                              "same_-" = TRUE,
                              NA
  )

  ScanBamParam(which = locus_granges,
               mapqFilter = quality_threshold,
               flag = scanBamFlag(isMinusStrand=minus_strand_flag,
                                  isSecondaryAlignment=FALSE,
                                  isNotPassingQualityControls=FALSE,
                                  isSupplementaryAlignment=FALSE,
                                  isDuplicate=FALSE,
                                  isUnmappedQuery=FALSE,
               ))
}

#'
#' Calculate KN99 protein coding total
#'
#' @description Specifically for the 'old' pipeline with the 'old' annotations.
#'
#' @importFrom dplyr filter summarize pull
#'
#' @param htseq_filename path to the htseq output for a given sample
#'
#' @return The summed total of reads over protein coding genes (no nctr RNA)
#'
#' @export
htseq_proteinCodingTotal = function(htseq_filename){

  htseq_df = readHTSeqFile(htseq_filename)

  # note: there are other reads, numbered likely in the 100s, which are
  # ambiguously mapped, meaning they overlap features, to protein coding
  # sequences. Extracting these sequences is more trouble than it is worth,
  # though was a feature of a previous iteration of the pipeline.
  htseq_df %>%
    filter(startsWith(feature, "CKF44")) %>%
    summarize(total = sum(sample)) %>%
    pull(total)

}

#'
#' Extract notAlignedTotalPercent from novoalign log
#'
#' @importFrom dplyr filter select pull
#'
#' @inheritDotParams parseNovoalignLog
#'
#' @return unaligned reads / total reads
#'
#' @export
htseq_notAlignedTotalPercent = function(...){

  # parse the novoalign log file
  novo_df = parseNovoalignLog(...)

  # extract the not_aligned_total
  filter(novo_df, metric %in% c("No Mapping Found", "Homopolymer Filter")) %>%
    select(reads) %>%
    sum() ->
    not_aligned_total

  # divide not_aligned_total by the total read sequences,
  # round to 4 decimal places
  round(not_aligned_total /
          pull(filter(novo_df, metric == "Read Sequences"), reads), 4)

}

#' Calculate percentage of library accounted for by top expressed genes
#'
#' @importFrom utils head
#' @importFrom dplyr filter summarize pull arrange
#'
#' @inheritParams htseq_proteinCodingTotal
#'
#' @param num_genes used to select the top expressed n genes. Default is 25
#'
#' @return sum(top_n_genes) / total_counted_reads, rounded to 4 decimal places
#'
#' @export
htseq_libraryComplexity = function(htseq_filename, num_genes = 25){

  htseq_df = readHTSeqFile(htseq_filename)

  counted_total = htseq_df %>%
    filter(!startsWith(feature, "__")) %>%
    summarize(total = sum(sample)) %>%
    pull(total)

  top_5_sequence_total = htseq_df %>%
    filter(!startsWith(feature, "__")) %>%
    arrange(desc(sample)) %>%
    head(num_genes) %>%
    summarize(total = sum(sample)) %>%
    pull(total)

  round(top_5_sequence_total / counted_total, 4)

}

#' Calculate coverage over a given feature
#'
#' @importFrom GenomicRanges width
#' @importFrom Rsamtools BamFile PileupParam pileup
#'
#' @inheritParams strandedScanBamParam
#' @inheritParams geneGRanges
#' @inheritDotParams strandedScanBamParam
#'
#' @param bam_path path to a given samples alignment file (.bam)
#' @param index_suffix suffix to append to bampath. Default .bai
#' @param coverage_threshold minimum threshold to consider. Default is 0
#'
#' @return coverage of feature
#'
#' @export
locusCoverage = function(bam_path, locus_granges, library_strandedness,
                         index_suffix = ".bai",
                         coverage_threshold = 0, ...){

  # construct index path
  bam_index_path = paste0(bam_path, index_suffix)
  # check arguments
  for(arg in c(bam_path, bam_index_path)){
    if(!file.exists(arg)){
      stop(paste0("path not valid: ", arg, ". If this ends in .bai, then
                  it means you must first index the bam file. Make sure the
                  index is saved in the same directory as the bam"))
    }
  }

  sbp = strandedScanBamParam(locus_granges, library_strandedness, ...)

  # set parameter options
  p_param <- PileupParam(min_mapq = sbp@mapqFilter,
                         min_nucleotide_depth = coverage_threshold,
                         distinguish_strands=TRUE,
                         distinguish_nucleotides=FALSE)

  message("reading bam file...")
  bamfile = BamFile(bam_path, bam_index_path)

  message("calculating coverage...")
  coverage_df = pileup(bamfile,
                       scanBamParam = sbp,
                       pileupParam = p_param)

  message("done")
  length(unique(coverage_df$pos))/sum(width(locus_granges))
}

#' Get the log2cpm of a given locus in the htseq output
#'
#' @importFrom dplyr filter
#' @importFrom edgeR cpm
#'
#' @inheritParams readHTSeqFile
#'
#' @param gene_id gene id for which to return the log2cpm
#'
#' @return log2cpm of a given gene
#'
#' @export
htseq_locusLog2cpm = function(htseq_filename, gene_id){

  # calculate log2cpm on protein coding genes only
  htseq_df = suppressMessages(readHTSeqFile(htseq_filename) %>%
    filter(startsWith(feature, "CKF44")))

  log2cpm_mat = cpm(matrix(htseq_df$sample), log = TRUE)

  names(log2cpm_mat) = htseq_df$feature

  round(as.numeric(log2cpm_mat[gene_id]), 2)

}


# overexpressionFOW = function(){
#
# }
