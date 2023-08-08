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

  # TODO add support for forward stranded libraries
  # set some information for the ScanBamParam object below. gene_strand extracts
  # the +/- strand from the GRanges object
  # and the conditional sets the minus_strand_flag used in the ScanBamParam
  # constructor. This determines what reads are
  # returned -- either from a given strand, or from both
  gene_strand = as.character(unique(data.frame(locus_granges)$strand))

  # ensure the locus is entirely on the same strand, error out if not
  if(!length(gene_strand) == 1){
    message(paste0("The granges object contains features on both strands. ",
                   "Separate the granges object into sets so that any set is on ",
                   "the same strand and try again"))
  }

  minus_strand_flag = switch (paste(library_strandedness, gene_strand, sep="_"),
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
#' @param markers markers in library. Default CNAG_NAT, CNAG_G418
#'
#' @return log2cpm of a given gene
#'
#' @export
htseq_locusLog2cpm = function(htseq_filename, gene_id,
                              markers = c("CNAG_NAT", "CNAG_G418")){

  # calculate log2cpm on protein coding genes only
  htseq_df = suppressMessages(
    readHTSeqFile(htseq_filename) %>%
      filter(feature %in% markers | startsWith(feature, "CKF44")))

  log2cpm_mat = cpm(matrix(htseq_df$sample), log = TRUE)

  names(log2cpm_mat) = htseq_df$feature

  round(as.numeric(log2cpm_mat[gene_id]), 2)

}


# overexpressionFOW = function(){
#
# }

#'
#' QC a crypto run output by the novoalign+htseq pipeline
#'
#' @description coverage and log2cpm are both over the annotated CDS
#'
#' @importFrom dplyr filter pull select rename mutate left_join
#' @importFrom parallel makeForkCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom iterators iter
#' @import foreach
#'
#' @inheritParams geneGRanges
#' @inheritParams createNovoalignPipelineSamplesheet
#'
#' @param pipeline_output_dirpath path to the directory which stores the
#'   subdirectories align, count and logs,
#'   eg /mnt/scratch/rnaseq_pipeline/pipeline_out/run_5500
#' @param markers a list of markers. must be in the counts and genome
#'   annotations. default is c("NAT", "G418")
#' @param bam_suffix suffix appended to the bam files. default is
#'   "_sorted_aligned_reads_with_annote.bam"
#' @param novolog_suffix suffix appended to log files. default is
#'   "_novoalign.log"
#' @param exon_counts_suffix suffix appended to exon count files. default is
#'   '_read_count.tsv'
#' @param cds_counts_suffix suffix appended to cds count files. default is
#'   '_read_count.tsv'
#' @param num_nodes number of cpus(by slurm definition)/threads(on your local).
#'   the argument in the parallel function is nnodes, hence the name of the
#'   argument. Default is 10
#'
#' @return a dataframe, long format, with columns fastqFileNumber, perturbed
#'   locus coverage/log2cpm, marker coverage/log2cpm and the library
#'   quality metrics
#'
#' @export
novoalignPipelineQC = function(meta_df,
                               pipeline_output_dirpath,
                               annote_obj_path,
                               markers = c("NAT", "G418"),
                               bam_suffix = "_sorted_aligned_reads_with_annote.bam",
                               novolog_suffix = "_novoalign.log",
                               exon_counts_suffix = "_read_count.tsv",
                               cds_counts_suffix = '_read_count_cds.tsv',
                               num_nodes = 10){

  message(paste0("WARNING: this function will not work on a windows machine ",
                 "due to the parallelization method currently implemented. ",
                 "Additionally, threads/cpus must be on the same machine ",
                 "(eg on slurm nodes_per_task=1, cpus_per_task=8)"))

  # TODO removed the filtering on genotype1 having cnag/ckf44 in prefix
  # b/c of change to strain table. check that samples form eg daniel do not
  # slip through now

  # extract all perturbed loci in run
  geno1_loci = meta_df %>%
    filter(strain != 'TDY451', strain != 'TDY450') %>%
    pull(genotype1) %>%
    as.character() %>%
    unique()

  geno2_loci = meta_df %>%
    filter(!is.na(genotype2), genotype2 != '') %>%
    pull(genotype2) %>%
    as.character() %>%
    unique()

  perturbed_loci_df = expand.grid(unique(meta_df$fastqFileNumber),
                                  c(geno1_loci, geno2_loci)) %>%
    dplyr::rename(fastqFileNumber = Var1, locus = Var2)

  # create qc df
  qc_df = meta_df %>%
    dplyr::select(fastqFileNumber, fastqFileName, libraryProtocol) %>%
    dplyr::rename(strandedness = libraryProtocol) %>%
    mutate(strandedness =
             ifelse(as.character(strandedness) == "E7420L",
                    "reverse", "unstranded")) %>%
    left_join(perturbed_loci_df)

  # initiate parallelization
  cl = parallel::makeForkCluster(nnodes = num_nodes)
  doParallel::registerDoParallel(cl)

  qc_df_mod =
    foreach(row=iter(qc_df, by='row'), .combine=rbind) %dopar% {

      # get paths to bam/count/log
      path_list = list(
        bam = file.path(pipeline_output_dirpath,
                        "align",
                        paste0(row$fastqFileName, bam_suffix)),
        exon_count = file.path(pipeline_output_dirpath,
                               "count",
                               paste0(row$fastqFileName, exon_counts_suffix)),
        cds_count = file.path(pipeline_output_dirpath,
                          "count",
                          paste0(row$fastqFileName, cds_counts_suffix)),
        novolog = file.path(pipeline_output_dirpath,
                            "logs",
                            paste0(row$fastqFileName, novolog_suffix))
      )

      # make granges for markers
      marker_granges = list(
        nat = geneGRanges(annote_obj_path,
                          "NAT",
                          feature = 'cds'),
        g418 = geneGRanges(annote_obj_path,
                           "G418",
                           feature = 'cds')
      )

      # calculate perturbed locus metrics
      if(!is.na(row$locus)){

        # note this replacement is unnecessary after the strain table addition
        # on 20220830. kept for safety and markers
        locus_granges = geneGRanges(annote_obj_path,
                                    str_replace(row$locus, "CNAG", "CKF44"),
                                    feature = 'cds')

        row$perturbedCoverage = locusCoverage(path_list$bam,
                                              locus_granges,
                                              row$strandedness)
        row$perturbedLog2cpm = htseq_locusLog2cpm(path_list$cds_count,
                                                  str_replace(row$locus,
                                                              "CNAG", "CKF44"))
      }

      # nat metrics
      row$natCoverage = locusCoverage(path_list$bam,
                                      marker_granges$nat,
                                      row$strandedness)
      row$natLog2cpm = htseq_locusLog2cpm(path_list$cds_count, "CNAG_NAT")

      # g418 metrics
      row$g418Coverage = locusCoverage(path_list$bam,
                                       marker_granges$g418,
                                       row$strandedness)
      row$g418Log2cpm = htseq_locusLog2cpm(path_list$cds_count, "CNAG_G418")

      # lib quality metrics
      row$proteinCodingCounted = htseq_proteinCodingTotal(path_list$exon_count)
      row$notAlignedTotalPercent = htseq_notAlignedTotalPercent(path_list$novolog)
      row$libraryComplexity = htseq_libraryComplexity(path_list$exon_count)

      # return the modified row (will be rbind together)
      row
    }

  parallel::stopCluster(cl)

  # return the qc_df with filled entries
  qc_df_mod %>%
    dplyr::select(-c(fastqFileName, strandedness))

}

#'
#' add qc metrics to metadata
#'
#' @importFrom tidyr as_tibble
#' @importFrom iterators iter
#' @importFrom dplyr distinct filter pull select rename mutate left_join
#' @import foreach
#'
#' @description take the long form of the qc_df and add the appropriate metrics
#'   to the metadata in a 'wide' format
#'
#' @inheritParams createNovoalignPipelineSamplesheet
#' @param qc_df output of the function [brentlabRnaSeqTools::novoalignPipelineQC]
#'
#' @return a subset of the meta_df columns suitable for auto/manual auditing
#'
#' @export
addQcColsToMeta = function(meta_df, qc_df){

  if(length(setdiff(meta_df$fastqFileNumber, qc_df$fastqFileNumber)) != 0){
    stop(paste0("There are fastqFileNumbers in meta_df that are not in qc_df. ",
          "Filter one or the other and resubmit -- the same samples and no
          more should be in each."))
  }

  meta_qc_select = meta_df %>%
    dplyr::select(fastqFileNumber,
                  genotype1,
                  genotype2,
                  marker1,
                  marker2)

  meta_qc_filled = foreach(
    row = iterators::iter(meta_qc_select, by = 'row'),
    .combine = 'rbind') %do% {

      ffn = row$fastqFileNumber
      genotype1 = as.character(row$genotype1)
      genotype2 = as.character(row$genotype2)

      genotype1_df = qc_df %>%
        filter(fastqFileNumber == ffn,
               locus == genotype1) %>%
        dplyr::select(fastqFileNumber, perturbedCoverage, perturbedLog2cpm) %>%
        dplyr::rename(genotype1Coverage = perturbedCoverage,
                      genotype1Log2cpm = perturbedLog2cpm)

      genotype2_df = qc_df %>%
        filter(fastqFileNumber == ffn,
               locus == genotype2) %>%
        dplyr::select(fastqFileNumber, perturbedCoverage, perturbedLog2cpm) %>%
        dplyr::rename(genotype2Coverage = perturbedCoverage,
                      genotype2Log2cpm = perturbedLog2cpm)

      marker_libqual_df = qc_df %>%
        filter(fastqFileNumber == ffn) %>%
        dplyr::select(fastqFileNumber, natCoverage, natLog2cpm, g418Coverage,
                      g418Log2cpm, proteinCodingCounted, notAlignedTotalPercent) %>%
        distinct(fastqFileNumber, .keep_all = TRUE)

      library_qual_df = qc_df %>%
        filter(fastqFileNumber == ffn)

      qc_select = marker_libqual_df %>%
        left_join(genotype1_df) %>%
        left_join(genotype2_df)

      row %>%
        as_tibble() %>%
        left_join(qc_select)
    }

  meta_qc_filled
}

#'
#' Audit a qc table with metrics added
#' @description audits the output of [brentlabRnaSeqTools::addQcColsToMeta]
#'
#' @param qc_table output of [brentlabRnaSeqTools::addQcColsToMeta]
#'
#' @return the qc_table with autoStatus and autoStatusDecomp columns added
#'
#' @export
autoAuditQcTable = function(qc_table){

  qc_table %>%
    mutate(autoStatus = 0) %>%
    mutate(autoStatus =
             ifelse(proteinCodingCounted <
                      kn99_novo_htseq_thresholds$proteinCodingCounted,
                    autoStatus + kn99_novo_htseq_status$proteinCodingCounted,
                    autoStatus)) %>%

    mutate(autoStatus =
             ifelse(notAlignedTotalPercent >
                      kn99_novo_htseq_thresholds$notAlignedTotalPercent,
                    autoStatus + kn99_novo_htseq_status$notAlignedTotalPercent,
                    autoStatus)) %>%

    mutate(autoStatus =
             ifelse((genotype1 != "CKF44_00000" &
                       genotype1Coverage >
                       kn99_novo_htseq_thresholds$perturbedCoverage) |
                      ((!is.na(genotype2) | genotype2 != "") &
                         !is.na(genotype2Coverage) & genotype2Coverage >
                         kn99_novo_htseq_thresholds$perturbedCoverage),
                    autoStatus + kn99_novo_htseq_status$perturbedCoverage,
                    autoStatus)) %>%

    mutate(autoStatus =
             ifelse((marker1 == "NAT" | marker2 == "NAT") &
                      (natCoverage < kn99_novo_htseq_thresholds$natExpectedCoverage |
                         natLog2cpm < kn99_novo_htseq_thresholds$natExpectedLog2cpm),
                    autoStatus + kn99_novo_htseq_status$natExpected,
                    autoStatus)) %>%

    mutate(autoStatus =
             ifelse((marker1 == "G418" | marker2 == "G418") &
                      g418Log2cpm < kn99_novo_htseq_thresholds$g418ExpectedLog2cpm,
                    autoStatus + kn99_novo_htseq_status$g418Expected,
                    autoStatus)) %>%

    mutate(autoStatus =
             ifelse((genotype1 == "CKF44_00000" |
                       (marker1 == "G418" & (is.na(marker2) | marker2 == ""))) &
                      (natCoverage > kn99_novo_htseq_thresholds$natUnexpectedCoverage &
                         natLog2cpm > kn99_novo_htseq_thresholds$natUnexpectedLog2cpm),
                    autoStatus + kn99_novo_htseq_status$natUnxpected,
                    autoStatus)) %>%

    mutate(autoStatus =
             ifelse((genotype1 == "CKF44_00000" |
                       (marker1 == "NAT" & (is.na(marker2) | marker2 == ""))) &
                      g418Log2cpm > kn99_novo_htseq_thresholds$g418UnxpectedLog2cpm,
                    autoStatus + kn99_novo_htseq_status$g418Unexpected,
                    autoStatus)) %>%

    mutate(autoStatusDecomp = unlist(map(autoStatus, decomposeStatus2Bit))) %>%
    # remove passing_sample -- current database doesn't have this as an option
    mutate(autoStatusDecomp = ifelse(autoStatusDecomp == 'passing_sample',
                                     NA, autoStatusDecomp)) %>%
    # add brackets to prevent database from considering this a number
    mutate(autoStatusDecomp = paste0("[", autoStatusDecomp, "]"))
}

#' decompose sums of powers of two to a list of the summed powers
#'
#' @description eg 18 = 2 + 16 decomposes to 1, 4
#'
#' @references yiming kang
#' \url{https://github.com/yiming-kang/rnaseq_pipe/blob/master/tools/utils.py}
#'
#' @param status an integer that represents the sum of powers of 2
#' @return a string of powers of 2 representing the bit, eg 1,4
#'
#' @export
decomposeStatus2Bit = function(status){

  #TODO can use negative numbers for "passing reasons"!!

  status_decomp = list()

  if(is.na(status)){
    status_decomp = "status_NA"
  } else if(status == 0){
    status_decomp = "passing_sample"
  } else if(status > 0){
    for(i in seq(floor(log2(status)),0)){
      if ((status -2**i) >= 0){
        status_decomp = append(status_decomp, i)
        status = status - 2**i
      }
    }
    status_decomp = paste(sort(unlist(status_decomp)), collapse=",")
  }
  status_decomp
}
