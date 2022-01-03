#' filter low replicate parameters from metadata
#'
#' @importFrom dplyr filter pull
#' @importFrom stringr str_remove str_trim str_split
#' @importFrom stats model.matrix
#'
#' @description given a model formula, remove samples with less than a specified number of replicates from the metadata
#'
#' @param metadata_df a data frame that contains at least the model paramters of interest
#' @param design_formula an R formula, eg ~libraryDate+treatment, of parameters contained in the metadata_df
#' @param replicate_threshold the number of replicates below which samples will be removed. Default is 2
#'
#' @return the input metadata with samples in replicate groups with less than the specified thershold filtered out
#'
#' @export
fltrLowReplicateParams = function(metadata_df, design_formula, replicate_threshold=2){

  # TODO write test and error handling for this function

  # create a character vector of the formula parameters
  # eg 'librarydate', 'treatment' if the original formula is ~librarydate+treatment
  formula_vector = str_trim(unlist(str_split(as.character(design_formula)[[2]], "\\+")))

  low_rep_params = -1
  while(length(low_rep_params != 0)){

    # create a model matrix and summary of number of replicates in each of the parameters
    model_matrix = model.matrix(design_formula, metadata_df)
    mm_summary_df = tibble(model_params = colnames(model_matrix), replicate_tally = colSums(model_matrix))

    # get the model parameters' factor levels with samples that fall in replicate groups below a specified number
    low_rep_params = mm_summary_df %>%
      filter(replicate_tally < replicate_threshold) %>%
      pull(model_params)

    low_rep_params = str_remove(low_rep_params, formula_vector)

    # remove samples with the factor levels in the parameters specified in low_rep_params
    samples_to_remove = unlist(lapply(formula_vector, function(x)
      filter(metadata_df, !!rlang::sym(x) %in% low_rep_params) %>% pull(FASTQFILENAME)))

    metadata_df = droplevels(filter(metadata_df, !FASTQFILENAME %in% samples_to_remove))
  }

  return(metadata_df)

} # end fltrLowReplicateParams()
