#'
#' Given expresion and binding signal, plot the rank response a la
#' doi:
#'
#' @param expression_df a two column dataframe. One column must be
#'   called 'gene' and must store gene identifiers (doesn't matter what, as long
#'   as they are consistent both within the expression list and between the
#'   expression list and the binding list). The second column identifies the
#'   experiment and stores the expression value (eg, log2FC)
#'
#' @param binding_df a list of dataframes with two columns. One column must be
#'   called 'gene' and must correspond to the expression_df gene column. The
#'   other column should be binding data. The column name should reference the
#'   source, eg 'calling_cards' or 'chip_chip'.
#'
#' @param tf the name of the TF -- this will be the plot title
#'
#' @return a ggplot object rank response plot
#'
#' @export
rank_response_plot = function(expression_df_list, binding_df_list, tf){

 message("Summarizing data by rank response...")
 rank_res_df_list = foreach(
   expr_name = names(expression_df_list)
  ) %do% {
    map(names(binding_df_list),
        ~sort_rank_mean_expr(expression_df_list[[expr_name]],
                             expr_name,
                             binding_df_list[[.]],
                             .))
  }

 rank_res_df = do.call('rbind', rank_res_df_list)

 rank_res_df %>%
   mutate(expression_source = as.factor(expression_source),
          binding_source = as.factor(binding_source),
          data = as.factor(data),
          rank = factor(as.character(rank))) %>%
   ggplot() +
     geom_line(aes(rank,
                   response_ratio,
                   color = binding_source,
                   linetype = expression_source))
}

#'
#' Given a vector length adn the size of each parition, create a vector which
#' represents those partitions.
#' @description There may be one more partition of size less than the equally
#' divided parts. eg if the vector length is 14 and the desired partitions are
#' size 3, then there will be 4 partitions of length 3 and one of length 2,
#' 1 1 1 2 2 2 3 3 3 4 4 4 5 5.
#'
#' @param vector_length the total lenght of the partition vector
#' @param equal_parts how large should each partition be?
#'
#'
#' @return a vector of `vector_length` divded into `equal_parts` with possibly
#'   one additional vector of size less than `equal_parts`.
create_partitions = function(vector_length, equal_parts = 100){
  c(rep(seq(1,vector_length/equal_parts),
        each=equal_parts),
    rep(floor(vector_length/equal_parts)+1, vector_length%%equal_parts))
}

sort_rank_mean_expr = function(expression_df, expression_src,
                               binding_df, binding_src,
                               lfc_thres, padj_thres,
                               rank_resolution = 10){

  expression_cols = c('gene', "log2FoldChange", "padj")
  binding_cols = c('gene', 'binding_signal')

  stopifnot(expression_cols %in% colnames(expression_df))
  expression_df = select(expression_df, all_of(expression_cols))
  stopifnot(binding_cols %in% colnames(binding_df))
  binding_df = select(binding_df, all_of(binding_cols))

  message(paste0("the overlap between the expression and binding ",
                 sprintf("gene sets is: %s out of %s in the expression data ",
                         length(union(expression_df$gene, binding_df$gene)),
                         nrow(expression_df)),
                 sprintf("and %s in the binding data", nrow(binding_df))))

  inner_join(expression_df, binding_df, by = "gene") %>%
    arrange(binding_signal) %>%
    mutate(rank = create_partitions(nrow(.), rank_resolution)*rank_resolution) %>%
    group_by(rank) %>%
    summarise(group_ratio = sum(abs(log2FoldChange) > lfc_thres & padj <= padj_thres)) %>%
    mutate(response_ratio = (cumsum(group_ratio)/rank)) %>%
    mutate(binding_source = binding_src) %>%
    mutate(expression_source = expression_src)
}
