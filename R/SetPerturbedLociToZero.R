#' Zero out perturbed expression
#'
#' @description set expression of the gene in genotype1 for a particular
#'   sample to zero
#'
#' @note the error handling isn't good here yet -- it will error if there is a
#'   perturbed loci that is not in the rownames of the object (eg, if it is filtered
#'   due to low expression)
#'
#' @import foreach
#'
#' @param object a deseq data set column genotype1 in colData
#'   (a brentlabRnaSeqSet also works)
#'
#' @export
setPerturbedLociToZero = function(object){

  count_matrix = counts(object)

  if(!startsWith(rownames(counts(object))[1], "CKF44")){
    stop("This function is only set up for single perturbation 90 min induction
         currently. Rownames of object/counts must be the KN99 gene ids, with
         the correct CKF44 prefix")
  }

  mod_count_matrix = foreach(
    colname = colnames(object),
    .combine = 'cbind',
    .inorder = TRUE
  ) %dopar% {
    perturbed_locus = as.character(unique(object[,as.character(colname)]$genotype1))
    col_matrix = counts(object)[,as.character(colname), drop = FALSE]
    if (perturbed_locus != "CKF44_00000"){
      perturbed_locus = str_replace(perturbed_locus, "CNAG", "CKF44")
      col_matrix[perturbed_locus,] = 0
    }
    col_matrix
  }

  mode(mod_count_matrix) = "integer"

  mod_count_matrix = mod_count_matrix[,colnames(object)]

  assay(object) = mod_count_matrix

  stopifnot(identical(colnames(counts(object)), colData(object)$fastqFileName))

  object
}
