# TODO parameterize the column headeing names for the binding/expression dfs
# TODO check that the lists are named, and that the names are unique within/
# between lists

#'
#' Given expresion and binding signal, plot the rank response a la
#' doi:
#'
#' @import foreach
#' @importFrom purrr map
#' @importFrom dplyr mutate
#' @import ggplot2
#'
#' @param expression_df_list a NAMED list of gene expression dataframes. Each
#'   dataframe must have the following columns, at minimum:
#'   c('gene', "log2FoldChange", "padj"). Names of items in the list must
#'   correspond to the data source (eg, 'kemmeren', 'brentlab', etc)
#' @param binding_df_list a NAMED list of gene expression dataframes. Each
#'   dataframe must have the following columns, at minimum:
#'   c('gene', 'binding_signal'). Names of items in the list must
#'   correspond to the data source (eg, 'chip-chip', 'calling-cards', etc)
#' @param tf the name of the TF -- this will be the plot title
#' @param lfc_thres expression log2foldchange threshold. default is 0
#' @param padj_thres expression padj threshold. default is .05
#'
#' @return a ggplot object rank response plot
#'
#' @export
rank_response_plot = function(expression_df_list, binding_df_list,
                              tf, lfc_thres = 0, padj_thres = .05){

 message("Summarizing data by rank response...")
  rank_res_df = foreach(
   expr_name = names(expression_df_list),
   .combine = 'rbind'
  ) %do% {
    df_list = map(names(binding_df_list),
        ~sort_rank_mean_expr(expression_df_list[[expr_name]],
                             expr_name,
                             binding_df_list[[.]],
                             .,
                             lfc_thres,
                             padj_thres))

    do.call('rbind', df_list)
  }

 rank_res_df %>%
   mutate(expression_source = as.factor(expression_source),
          binding_source = as.factor(binding_source)) %>%
   ggplot() +
     geom_line(aes(rank,
                   response_ratio,
                   color = binding_source,
                   linetype = expression_source)) +
   ggtitle(tf) +
   scale_x_continuous(breaks = seq(0,150,10)) +
   scale_y_continuous(breaks = seq(0,1,.1)) +
   theme(text = element_text(size = 20),
         axis.text.x = element_text(angle = 45, hjust =1))
}

#'
#' Given a vector length adn the size of each parition, create a vector which
#' represents those partitions.
#' @description There may be one more partition of size less than the equally
#' divided parts. eg if the vector length is 14 and the desired partitions are
#' size 3, then there will be 4 partitions of length 3 and one of length 2,
#' 1 1 1 2 2 2 3 3 3 4 4 4 5 5.
#'
#' @import dplyr
#'
#' @param vector_length the total lenght of the partition vector
#' @param equal_parts how large should each partition be?
#'
#'
#' @return a vector of `vector_length` divded into `equal_parts` with possibly
#'   one additional vector of size less than `equal_parts`.
#'
#' @export
create_partitions = function(vector_length, equal_parts = 100){
  c(rep(seq(1,vector_length/equal_parts),
        each=equal_parts),
    rep(floor(vector_length/equal_parts)+1, vector_length%%equal_parts))
}

#'
#' Create a rank-response data frame
#'
#' @inheritParams rank_response_plot
#' @param expression_df a three column dataframe with the column names
#'   c('gene', 'log2FoldChange', 'padj') which describes expression data. It
#'   may have more columns, but those three must exist
#' @param binding_df a two column dataframe dataframe with the columns
#'   c('gene', 'binding_signal'). It may have more columns, but those must exist
#' @param expression_src the source of the data, eg 'kemmeren' or 'brentlab'
#' @param binding_src the source of the data, eg 'chip-chip' or 'calling-cards'
#' @param rank_resolution how to bin ranks. Default is 10, which means that
#'   ranks are binned into groups of 10 (eg the first group are the 10 highest)
#' @param num_rows the number of rows to return. Default is 15.
#'   For instance, if the rank_resolution is 10, and you only want to look at
#'   the first 150 ranked genes, then num_rows is 15.
#'
#' @return a rank_response dataframe
#'
#' @export
sort_rank_mean_expr = function(expression_df, expression_src,
                               binding_df, binding_src,
                               lfc_thres, padj_thres,
                               rank_resolution = 10,
                               num_rows = 15){

  expression_cols = c('gene', "log2FoldChange", "padj")
  binding_cols = c('gene', 'binding_signal')

  stopifnot(expression_cols %in% colnames(expression_df))
  expression_df = select(expression_df, all_of(expression_cols)) %>%
    filter(complete.cases(.))
  stopifnot(binding_cols %in% colnames(binding_df))
  binding_df = select(binding_df, all_of(binding_cols)) %>%
    filter(complete.cases(.))

  message(paste0("the overlap between the complete cases in the expression and binding ",
                 sprintf("gene sets is: %s out of %s in the expression data ",
                         length(intersect(expression_df$gene, binding_df$gene)),
                         nrow(expression_df)),
                 sprintf("and %s in the binding data", nrow(binding_df))))

  inner_join(expression_df, binding_df, by = "gene") %>%
    arrange(binding_signal) %>%
    mutate(rank = create_partitions(nrow(.), rank_resolution)*rank_resolution) %>%
    group_by(rank) %>%
    summarise(group_ratio = sum(abs(log2FoldChange) > lfc_thres & padj <= padj_thres)) %>%
    mutate(response_ratio = (cumsum(group_ratio)/rank)) %>%
    mutate(binding_source = binding_src) %>%
    mutate(expression_source = expression_src) %>%
    head(num_rows)
}

