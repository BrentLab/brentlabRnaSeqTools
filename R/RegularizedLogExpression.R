#' transform normalized counts into regularized log expression
#'
#'
#' @import foreach
#' @importFrom SummarizedExperiment SummarizedExperiment rowRanges
#' @importFrom S4Vectors SimpleList
#'
#' @param object a brentlabRnaSeqSet object
#'
#' @return a brentlabRnaSeqSetTransform object
#'
#' @export
rleTransform = function(object){

  # TODO check that class is brentlabRnaSeqSet
  # TODO make this a generic
  # TODO check that replicate_group exists in colData
  # TODO remove parameter effect and add that as assay

  if(is.null(object$replicate_group)){
    stop("You must set replicate_group as a feature like so:
         object$replicate_group = replicate_group_vector, where
         replicate_group_vector is a factored vector of replicate groups")
  }

  rle_matrix = foreach(
    replicate_group = unique(colData(object)$replicate_group),
    .combine = 'cbind'
  ) %do% {
    replicate_group_norm_counts =
      counts(object[,colData(object)$replicate_group == replicate_group],
             normalized = TRUE)
    # if there are less than 3 samples in the replicate group, do not calculate.
    if(ncol(replicate_group_norm_counts) < 3){

      message(paste0("The following replicate set has less than 3 reps: ",
                     colnames(replicate_group_norm_counts),
                     ". It is being returned as a count matrix
                     of NAs (with appropriate column headings"))

      replicate_group_rle = replicate_group_norm_counts %>%
        mutate_all(~replace(., is.numeric(.), NA))

      # else, calculate rle for replicate group and summarize
    } else{
      replicate_group_rle = calculateRLE(replicate_group_norm_counts,
                                         log2_transformed_flag = FALSE)
    }
    replicate_group_rle
  }

  rle_matrix = rle_matrix[,colData(object)$fastqFileName]

  se = SummarizedExperiment(
    colData = as_tibble(colData((object))),
    rowRanges = rowRanges(object),
    assays  = SimpleList(rle = rle_matrix))

  brentlabRnaSeqSetTransform(se)

}

#' calculate RLE of a numeric dataframe
#'
#' @importFrom dplyr mutate_all
#'
#' @param counts_df gene by samples dataframe of raw counts or logged counts (see paramter logged)
#' @param log2_transformed_flag Default FALSE set to true if log2 transformed counts are passed
#'
#' @return rle dataframe with genes x samples. Values are the logged differences from the gene-wise medians
#'
#' @export
calculateRLE = function(counts_df, log2_transformed_flag = FALSE){

  if(!isNumeric(counts_df)){
    stop("counts_df must have all numeric columns")
  }

  if(ncol(counts_df) < 3){
    rle_table_full = counts_df %>%
      mutate_all(~replace(., is.numeric(.), NA))

    message(paste0("The following replicate set has less than 3 reps: ",
                   colnames(counts_df), ". It is being returned as a count matrix
                   of NAs (with appropriate column headings"))
  } else{

    # if log2_transformed_flag==TRUE, re-assign counts_df to log2_counts_df. else, add a pseudocount and log2
    ifelse(log2_transformed_flag,
           assign('log2_counts_df', counts_df),
           assign('log2_counts_df', log2(counts_df + 1)))

    # calculate median expression for each gene across samples
    gene_wise_medians = apply(log2_counts_df, 1, median, na.rm=TRUE)

    # calculate deviations from the median
    rle_table_full = sweep(log2_counts_df, 1, gene_wise_medians, '-')
  }

  return(rle_table_full)

} # end calculateRLE()

#' calculate medians across rows of dataframe
#'
#' @param count_df could be any numeric dataframe, but in this context it will typically be a count (raw or log2) df
#' @return a vector of row-wise medians (length == nrow of input df)
#'
#' @export
calculateGeneWiseMedians = function(count_df){

  # calculate median expression for each gene across samples
  all_gene_medians = apply(count_df, 1, median)

  return(all_gene_medians)

} # end calculateGeneWiseMedians

#' rleSummary calculates summary statistics of rleFullTable
#'
#' @importFrom matrixStats iqr
#' @importFrom stats median
#'
#' @param rle_table the output of rleFullTable
#'
#' @return a dataframe sample by rle summary statistics
#'
#' @export
rleSummary = function(rle_table){

  # assemble table
  rle_table_summary = data.frame(
    rle_median_deviation = apply(rle_table, 2, median, na.rm=TRUE),
    rle_iqr = apply(rle_table, 2, iqr, na.rm=TRUE))

  rownames(rle_table_summary) = colnames(rle_table)

  rle_table_summary

} # end rleSummary()

#'
#' calculate RLE by replicate groups
#'
#' @param object a brentlabRnaSeqSet with colData replicate_group
#'
#' @return a list of dataframes for each replicate group in replicateS_sample_list, each with dimensions gene x sample. values are RLE of the gene in a given sample
#'
#' @export
addRleSummaryStats = function(object){

  # TODO check that replicate_group exists in colData
  # TODO check that is class brentlabRnaSeqTransform
  # and has slot rleExpr

  rle_vector = foreach(
    replicate_group = unique(colData(object)$replicate_group),
    .combine = 'c'
  ) %dopar% {
    # replicate_group_rle =
    #   rleExpr(object[,colData(object)$relpicate_group == replicate_group])

    # if there are less than 3 samples in the replicate group, do not calculate.
    if(ncol(replicate_group_rle) < 3){

      message(paste0("The following replicate set has less than 3 reps: ",
                     colnames(replicate_group_rle),
                     ". It is being returned as a count matrix
                   of NAs (with appropriate column headings"))

      # rle_summary_vector = rep(list(list(rle_median_deviation = NA, rle_iqr = NA)), ncol(raw_counts))
      names(rle_summary_vector) = colnames(rle_summary_vector)

    # else, calculate rle for replicate group and summarize
    } else{

      # summarize and return nested list where each sample has two sublists
      # called rle_median_deviation and rle_iqr
      rle_summary_vector = as.list(rleSummary(replicate_group_rle))
    }
    rle_summary_vector
  }

  df = as.data.frame(do.call(rbind, rle_vector))
  # df[[unique_identifier]] = rownames(df)

  # addColData(object, df)

  object

}

#'
#' calculate RLE by replicate groups
#'
#' @importFrom tibble as_tibble
#' @importFrom purrr map
#' @importFrom dplyr select all_of
#'
#' @param replicates_vector a list of lists where each sublist represents a replicate group. Entries must be a metadata
#'                               parameter, such as fastqFileName, that corresponds to the columns of the counts.
#'                               Suggestion: use something like these dplyr functions to create the list of lists group_by() %>% group_split %>% pull(fastqFileName)
#' @param gene_quants a gene x sample dataframe with values as some sort of gene quantification (eg normed counts, or log2(norm_counts) with some effect removed), possibly already logged (@see already_logged_flag)
#' @param log2_transformed_flag a boolean where TRUE means the counts are already in log2 space
#'
#' @seealso \code{\link{rlePlotCompareEffectRemoved}} to plot the norm counts and removedEffect 'counts' on the same plot
#'
#' @return a list of dataframes for each replicate group in replicateS_sample_list, each with dimensions gene x sample. values are RLE of the gene in a given sample
#'
#' @export
rleByReplicateGroup = function(replicates_vector, gene_quants, log2_transformed_flag){

  gene_quants_df = as_tibble(gene_quants)

  map(replicates_vector,
      ~calculateRLE(select(gene_quants_df, all_of(.)),
                    log2_transformed_flag=log2_transformed_flag))

}

# TODO fix these plotting functions

#' plot RLE for a given column filter
#'
#' @import DESeq2
#'
#' @param deseq_object a deseq object with results from the DESeq() call
#' @param model_matrix the deseq_object model matrix
#' @param column_filter a vector of fastqFileNames (or whatever the columns -- samples -- are called)
#' @param title of the plots
#' @return list with slots norm_count_rle and effect_removed_rle
#'
#' @export
rlePlot = function(deseq_object, model_matrix, column_filter, title){

  norm_counts = counts(deseq_object, normalize=TRUE)

  fltr_norm_counts = norm_counts[ , column_filter]

  effect_removed_counts = removeParameterEffects(deseq_object, model_matrix)

  fltr_effect_removed_counts = effect_removed_counts[ ,column_filter]

  norm_count_rle = rlePlot_helper(fltr_norm_counts, log2_transformed_flag=FALSE, paste(title, 'Norm Counts', sep=" - "))
  effect_removed_rle = rlePlot_helper(fltr_effect_removed_counts, log2_transformed_flag=TRUE, paste(title, 'Effect Removed', sep=' - '))

  return (list('norm_count_rle' = norm_count_rle, 'effect_removed_rle' = effect_removed_rle))


}

#' the actual plotting function for rlePlot
#'
#' @importFrom readr read_tsv
#' @importFrom tidyr pivot_longer
#' @import ggplot2
#'
#' @param count_df counts in gene x sample
#' @param log2_transformed_flag boolean where TRUE indicates the counts are in log2 space
#' @param title title of the output plot
#' @param gene_id_path path to gene_ids. Default set to read from Sys.getenv("GENE_ID_PATH")
#'
#' @return a ggplot
#'
rlePlot_helper = function(count_df, log2_transformed_flag, title, gene_id_path = Sys.getenv("GENE_ID_PATH")){


  rle_full_table = calculateRLE(count_df, log2_transformed_flag=log2_transformed_flag)

  # TODO MAKE THIS FAR MORE RESILENT -- VERY ERROR PRONE
  gene_id = read_tsv(gene_id_path, col_names = 'gene_id')[1:6967,]
  rle_full_table$gene_id = gene_id

  rle_full_table %>%
    pivot_longer(!gene_id, names_to = 'FASTQFILENAME', values_to = "RLE") %>%
    ggplot()+
    geom_boxplot(aes(FASTQFILENAME, RLE), outlier.shape=NA)+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
    xlab('Samples')+
    ylab('Deviation from Median')+
    coord_cartesian(ylim = c(-3,3))+
    scale_y_continuous(n.breaks = 20)+
    geom_hline(yintercept = 0)+
    ggtitle(title)
}

#' plots output of rleSummaryByReplicateGroup
#'
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate left_join bind_rows
#' @import ggplot2
#'
#' @param norm_counts_rle output of calculateRLE (maybe one of the sublists in rleByReplicateGroup())
#' @param removed_effect_rle see norm_counts_rle, but after removing some batch effects
#' @param metadata_df metadata with at least FASTQFILENAME and LIBRARYDATE
#' @param title title of the plot
#'
#' @return a ggplot with both the norm counts (more transparent) and removedEffect 'counts' on the same plot
#'
#' @export
rlePlotCompareEffectRemoved = function(norm_counts_rle, removed_effect_rle, metadata_df, title){

  norm_counts_rle %>%
    pivot_longer(colnames(.), names_to="FASTQFILENAME", values_to="RLE") %>%
    mutate(quant_type="norm_counts") %>%
    bind_rows(
      (removed_effect_rle %>%
         pivot_longer(colnames(.), names_to="FASTQFILENAME", values_to="RLE") %>%
         mutate(quant_type="removed_effect"))) %>%
    left_join(metadata_df) %>%
    mutate(LIBRARYDATE = as.factor(LIBRARYDATE)) %>%
    ggplot() +
    geom_boxplot(aes(FASTQFILENAME, RLE, fill=LIBRARYDATE, color=quant_type, alpha=quant_type), outlier.shape = NA) +
    scale_color_manual(values=c("#999999", "#000000")) +
    scale_alpha_manual(values=c(.1, 1)) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
    ylim(-3,3)+
    ggtitle(title)

}

#'
#' composite plot of all rle_stats in one plot
#' @note no title on graphs. use the names of the returned object to title in presentation
#'
#' @param rle_df a samples x rle stats df with at minimum columns SAMPLE_DEVIATION_MEDIAN, ABS_SAMPLE_DEVIATION_MEDIAN, INTERQUARTILE_RANGE
#'
#' @return a plot with three horizontal panels for each of the rle stats
#'
# plotRLEhistograms = function(rle_df){
#   rle_df %>%
#     pivot_longer(-FASTQFILENAME, names_to="rle_stat", values_to="value") %>%
#     mutate(rle_stat = factor(rle_stat, levels=c("SAMPLE_DEVIATION_MEDIAN", "ABS_SAMPLE_DEVIATION_MEDIAN", "INTERQUARTILE_RANGE"))) %>%
#     ggplot() +
#     geom_histogram(aes(value))+
#     theme(axis.title.x=element_blank())+
#     facet_wrap('rle_stat', scales="free_x", dir='v')
# }

#'
#' Remove the libraryDate effects from KN99 dds objs
#'
#' @param dds a deseq data object with either a formula or matrix in the design
#'   note that if libraryDate is not included in the design, then this function
#'   doesn't remove anything
#' @param replicate_colname which column to use as the replicate grouping. If
#'   your replicates are described by more than 1 column, for this purpose,
#'   group them to create 1 column. does not need to be the same as what is in
#'   the design slot of the dds
#'
#' @return a list with both the unremoved and removed effect rle and summary
#'
#' @export
removeLibdateByReplicate = function(dds, replicate_colname){

  message(paste0("WARNING! This removes only the libraryDate effect. If there are ",
                 "no columns in the design matrix that start with library, this ",
                 "removes NOTHING."))

  mm = NA

  if(is_formula(design(dds))){
    mm = model.matrix(design(dds), as_tibble(colData(dds)))
  } else if(is.matrix(design(dds))){
    mm = design(dds)
  }

  if(!is.matrix(mm)){
    stop(paste0("Can't extract design matrix from DeseqDataSet -- ",
                "design slot filled ",
                "with neither a formula or a matrix."))
  }

  unwanted_indicies = which(startsWith(colnames(mm), "library"))

  replicate_split =
    as_tibble(colData(dds)) %>%
    group_by(!!rlang::sym(replicate_colname)) %>%
    group_split()

  replicate_names =
    lapply(replicate_split,
           function(x)
             as.character(unique(pull(x,!! rlang::sym(replicate_colname)))))

  replicate_sample_list =
    lapply(replicate_split, function(x)
      as.vector(pull(x,fastqFileName)))

  names(replicate_sample_list) = replicate_names


  removed_effect_log_norm = removeParameterEffects(dds, unwanted_indicies)

  removed_effect_rle_by_replicate_group =
    rleByReplicateGroup(replicate_sample_list,
         removed_effect_log_norm,
         log2_transformed_flag = TRUE)

  removed_effect_rle_summary =
    as_tibble(bind_rows(lapply(removed_effect_rle_by_replicate_group,
                               rleSummary)),
              rownames = "fastqFileName") %>%
    left_join(as_tibble(colData(dds)),
              by='fastqFileName')

  log_norm_rle = calculateRLE(counts(dds),
                              log2_transformed_flag = FALSE)

  log_norm_rle_by_replicate_group =
    rleByReplicateGroup(replicate_sample_list,
                        log_norm_rle,
                        log2_transformed_flag = TRUE)

  log_norm_rle_summary =
    as_tibble(bind_rows(lapply(log_norm_rle_by_replicate_group,
                               rleSummary)), rownames = "fastqFileName") %>%
    left_join(as_tibble(colData(dds)),
              by='fastqFileName')

  list(
    with_libdate_effect = list(
      rle = log_norm_rle,
      summary = log_norm_rle_summary
    ),
    without_libdate_effect = list(
      rle = removed_effect_log_norm,
      summary = removed_effect_rle_summary
    )
  )
}

