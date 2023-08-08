# TODO add example

#' @title snpEff Summary to Long
#'
#' @description Transform a snpEff tsv summary frame to long format
#'
#' @param snpeff_summary_path path to the .txt summary output of snpEff
#' @param sample_id unique identifier of this sample -- sample_id becomes a column
#' @param variant_caller the variant caller used to create the vcf
#' @param effect_threshold minimum count in a variant_effect factor level (see
#'  the seealso link below for the snpEff docs if variant_effect factor doesn't
#'  ring a bell)
#' @param impact_threshold minimum count in a variant_impact factor level (see
#' the command in effect_threshold)
#'
#' @return a tibble with the fields GeneId, TranscriptId, Biotype, effect,
#'   effect_score, impact_score, impact, sample_id, variant_caller
#'
#' @examples
#'
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#'
#' @seealso \url{https://pcingola.github.io/SnpEff/se_outputsummary/}
#'
#' @rdname snpeff_to_long
#'
#' @export
#'
#' @importFrom readr read_tsv
#' @importFrom dplyr select all_of starts_with mutate filter
#' @importFrom tidyr pivot_longer
#' @importFrom rlang sym
#' @importFrom stringr str_extract
snpeff_summary_to_long = function(snpeff_summary_path,
                                  sample_id, variant_caller,
                                  effect_threshold = 0,
                                  impact_threshold = 0){

  # final field order of the output dataframe
  final_col_order = c('GeneId','TranscriptId','BioType','effect',
                      'effect_score', 'impact_score', 'impact',
                      'sample_id', 'variant_caller')

  # read in the snpeff summary. The first line is a comment and is skipped
  df = readr::read_tsv(snpeff_summary_path, skip = 1)

  # transform the snpeff summary into a long format
  df %>%
    # first, pivot the variant_effect columns. The result will be a column
    # with the variant_effect levels and column with the count for that level
    tidyr::pivot_longer(cols = dplyr::starts_with('variants_effect'),
                 names_to = 'effect',
                 values_to = 'effect_score') %>%
    # pivot again to create a column, impact, with the variant_impact levels
    # and a column impact_score with the counts of variants which have a given
    # impact level
    tidyr::pivot_longer(cols = dplyr::starts_with("variants_impact"),
                 names_to = "impact",
                 values_to = "impact_score") %>%
    dplyr::filter(effect_score > effect_threshold,
           impact_score > impact_threshold) %>%
    # add columns sample_id and variant_caller
    dplyr::mutate(sample_id = sample_id,
           variant_caller = variant_caller) %>%
    # select/order the fields in the output
    dplyr::select(dplyr::all_of(final_col_order)) %>%
    group_by(GeneId) %>%
    arrange(desc(impact_score), desc(effect_score), .by_group = TRUE)

}
