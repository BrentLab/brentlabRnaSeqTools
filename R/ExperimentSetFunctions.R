# Environmental Perturbation Sets ----------------------------------------------

#' filter combined_df for environmental perturbation sample set
#'
#' @import magrittr
#' @importFrom dplyr mutate_if mutate filter across rename_with
#' @importFrom stringr str_remove
#' @importFrom tidyr replace_na
#'
#' @param blrs a brentlabRnaSeqSet object
#'
#' @return environmental pertubation set
#'
#' @export
createEnvPert_epWT = function(blrs){

  # filter
  ep_meta_fastqFileNumbers = extractColData(blrs) %>%
    filter(strain == 'TDY451',
           treatment %in% c("cAMP","noTreatment"),
           treatmentConc %in% c("20", 'noTreatmentConc'),
           experimentDesign == 'Environmental_Perturbation',
           purpose == "fullRNASeq",
           libraryProtocol == 'E7420L',
           !is.na(fastqFileName)) %>%
    pull(fastqFileNumber)

  ep_set = blrs[,colData(blrs)$fastqFileNumber %in%
                  ep_meta_fastqFileNumbers]

  experimental_conditions = alist(medium,
                                  atmosphere,
                                  temperature,
                                  treatment,
                                  treatmentConc,
                                  pH,
                                  timePoint)

  concat_treatment = extractColData(ep_set) %>%
    mutate(concat_treatment = as.factor(paste(!!!experimental_conditions,
                                              sep="_"))) %>%
    pull(concat_treatment)

  colData(ep_set)$concat_treatment = concat_treatment

  ep_set

}

#' create Environmental Perturbation Titration (TDY451 only) set
#'
#' @import magrittr
#' @importFrom dplyr mutate_if mutate filter across rename_with
#' @importFrom stringr str_remove
#' @importFrom tidyr replace_na
#'
#' @param blrs a brentlabRnaSeqSet object
#'
#' @return environmental pertubation set
#'
#' @export
createEnvPertSet_titrationWT = function(blrs){

  # filter
  ep_meta_fastqFileNumbers = extractColData(blrs) %>%
    filter(strain == 'TDY451',
           medium == 'RPMI',
           atmosphere == "noAtmosphere",
           temperature == 30,
           treatment %in% c("cAMP", "noTreatment"),
           pH == "noPH",
           experimentDesign %in%
             c('ep_cAMP_titration', 'Environmental_Perturbation'),
           purpose == "fullRNASeq",
           libraryProtocol == 'E7420L',
           !is.na(fastqFileName)) %>%
    pull(fastqFileNumber)

  ep_set = blrs[,colData(blrs)$fastqFileNumber %in% ep_meta_fastqFileNumbers]

  experimental_conditions = alist(medium,
                                  atmosphere,
                                  temperature,
                                  treatment,
                                  treatmentConc,
                                  pH,
                                  timePoint)

  concat_treatment = extractColData(ep_set) %>%
    mutate(concat_treatment = as.factor(paste(!!!experimental_conditions,
                                              sep="_"))) %>%
    pull(concat_treatment)

  colData(ep_set)$concat_treatment = concat_treatment

  ep_set

}

#' create Environmental Perturbation Perturbed set
#'
#' @import magrittr
#' @importFrom dplyr mutate_if mutate filter across rename_with
#' @importFrom stringr str_remove
#' @importFrom tidyr replace_na
#'
#' @param blrs a brentlabRnaSeqSet object
#'
#' @return environmental pertubation brentlabRnaSeqSet
#'
#' @export
createEnvPertSet_perturbed = function(blrs){

  # extract conditions in the perturbed samples ----
  non_wt_ep_meta = extractColData(blrs) %>%
    filter(strain != 'TDY451',
           treatment %in% c("cAMP", "noTreatment"),
           experimentDesign %in%
             c('Environmental_Perturbation', 'ep_cAMP_titration'),
           purpose == "fullRNASeq",
           libraryProtocol == 'E7420L',
           !is.na(fastqFileName))

  medium_conds = non_wt_ep_meta %>%
    pull(medium) %>%
    unique()

  temperature_conds = non_wt_ep_meta %>%
    pull(temperature) %>%
    unique()

  atmosphere_conds = non_wt_ep_meta %>%
    pull(atmosphere) %>%
    unique()

  treatment_conds = non_wt_ep_meta %>%
    pull(treatment) %>%
    unique()

  treatmentConc_conds = non_wt_ep_meta %>%
    pull(treatmentConc) %>%
    unique()

  timePoint_conds = non_wt_ep_meta %>%
    pull(timePoint) %>%
    unique()

  libraryDate_conds = non_wt_ep_meta %>%
    pull(libraryDate) %>%
    unique()

  # get all samples matching these conditions, regardless of geno ----
  set_ffn = extractColData(blrs) %>%
    filter(medium %in% medium_conds,
           temperature %in% temperature_conds,
           atmosphere %in% atmosphere_conds,
           treatment %in% treatment_conds,
           treatmentConc %in% treatmentConc_conds,
           timePoint %in% timePoint_conds,
           libraryDate %in% libraryDate_conds,
           experimentDesign %in%
             c('Environmental_Perturbation', 'ep_cAMP_titration'),
           purpose == "fullRNASeq",
           !is.na(fastqFileName)) %>%
    pull(fastqFileNumber)

  # filter the blrs obj ----
  ep_set = blrs[,colData(blrs)$fastqFileNumber %in% set_ffn]

  experimental_conditions = alist(medium,
                                  atmosphere,
                                  temperature,
                                  treatment,
                                  treatmentConc,
                                  pH,
                                  timePoint)

  concat_treatment = extractColData(ep_set) %>%
    mutate(concat_treatment = as.factor(paste(!!!experimental_conditions, sep="_"))) %>%
    pull(concat_treatment)

  colData(ep_set)$concat_treatment = concat_treatment

  ep_set

}

# 90 minuteInduction sets ------------------------------------------------------

#'
#' 90minuteInduction set Definition
#'
#' @description The current definition of the 90 minute induction dataset,
#'   according to the 2016 grant summary (loaded into environment, see head(grant_df)) -- single KO only
#'
#' @import magrittr
#' @importFrom dplyr mutate_if mutate filter bind_rows
#' @importFrom stringr str_remove str_detect
#' @importFrom tidyr replace_na
#'
#' @param blrs a brentlabRnaSeqSet object
#'
#' @return the set metadata -- single KO only
#'
#' @export
createNinetyMinuteInduction_2016grant = function(blrs){


  condition_fltr_metadata = extractColData(blrs) %>%
    filter(medium %in% c("DMEM"),
           temperature %in% c(37),
           atmosphere %in% c("CO2"),
           treatment %in% c("noTreatment"),
           otherConditions %in% c("noOtherConditions"),
           pH %in% c("noPH"),
           timePoint %in% c(90),
           purpose=="fullRNASeq",
           !is.na(fastqFileName),
           str_detect(genotype1, "CNAG"),
           is.na(genotype2) | genotype2 == "",
           perturbation1 =="deletion" | is.na(perturbation1) | perturbation1 == "")

  # TODO: COMBINE THE FILTER STATEMENTS INTO SINGLE FILTER

  # get wildtypes
  wt_induction_set = condition_fltr_metadata %>%
    filter(genotype1 == "CNAG_00000", strain == 'TDY451')

  # filter for genotypes in the grant summary
  # exclude genotypes labelled 2 (which means 2 failed attempts)
  fltr_grant_df = grant_df %>%
    filter(strain_status == "done")

  perturbed_induction_set = condition_fltr_metadata %>%
    filter(genotype1 %in% fltr_grant_df$genotype1 &
             (is.na(genotype2) | genotype2 == ""))

  # put the wt and filtered genotypes together
  set = bind_rows(wt_induction_set, perturbed_induction_set)

  out = blrs[,colData(blrs)$fastqFileNumber %in% set$fastqFileNumber]

  filterWtByExperimentalLibdate_90min(out)
}

#' ninetyMin set where both doubles are in the grant_df
#'
#' @import magrittr
#' @importFrom dplyr mutate_if mutate filter bind_rows
#' @importFrom stringr str_remove str_detect
#' @importFrom tidyr replace_na
#'
#' @param blrs a brentlabRnaSeqSet object
#'
#' @return the set metadata
#'
#' @export
createNinetyMinuteInduction_2016grantWithDoubles = function(blrs){

  condition_fltr_metadata = extractColData(blrs) %>%
    filter(medium %in% c("DMEM"),
           temperature %in% c(37),
           atmosphere %in% c("CO2"),
           treatment %in% c("noTreatment"),
           otherConditions %in% c("noOtherConditions"),
           pH %in% c("noPH"),
           timePoint %in% c(90),
           purpose=="fullRNASeq",
           !is.na(fastqFileName),
           str_detect(genotype1, "CNAG"),
           perturbation1 =="deletion" | is.na(perturbation1) | perturbation1 == "",
           perturbation2 =="deletion" | is.na(perturbation2) | perturbation2 == "")

  # TODO: COMBINE THE FILTER STATEMENTS INTO SINGLE FILTER
  # get wildtypes
  wt_induction_set = condition_fltr_metadata %>%
    filter(genotype1 == "CNAG_00000", strain == 'TDY451')

  double_ko_df = condition_fltr_metadata %>%
    filter(!is.na(genotype2) & genotype2 != "",
           (genotype1 %in% grant_df$genotype1 |
             genotype2 %in% grant_df$genotype1))

  double_genos = unique(c(unique(as.character(double_ko_df$genotype1)),
                          unique(as.character(double_ko_df$genotype2))))

  single_ko_for_doubles = condition_fltr_metadata %>%
    filter(is.na(genotype2) | genotype2 == "" & genotype1 %in% double_genos)

  # put the wt and filtered genotypes together
  set = bind_rows(wt_induction_set, single_ko_for_doubles, double_ko_df)


  out = blrs[,colData(blrs)$fastqFileNumber %in% set$fastqFileNumber]

  filterWtByExperimentalLibdate_90min(out)

}

#'
#' Samples with genotypes not in the grant_df
#' @description return a set with only those genotypes (including doubles) not in the grant_df.
#'
#'
#' @import magrittr
#' @importFrom dplyr mutate_if mutate filter bind_rows
#' @importFrom stringr str_remove str_detect
#' @importFrom tidyr replace_na
#'
#' @param blrs a brentlabRnaSeqSet object
#'
#' @return the set metadata
#'
#' @export
createNinetyMinuteInduction_non2016grant = function(blrs){

  condition_fltr_metadata = extractColData(blrs) %>%
    filter(medium %in% c("DMEM"),
           temperature %in% c(37),
           atmosphere %in% c("CO2"),
           treatment %in% c("noTreatment"),
           otherConditions %in% c("noOtherConditions"),
           pH %in% c("noPH"),
           timePoint %in% c(90),
           purpose =="fullRNASeq",
           !is.na(fastqFileName),
           str_detect(genotype1, "CNAG"),
           perturbation1 == "deletion" | is.na(perturbation1) | perturbation1 == "",
           perturbation2 == "deletion" | is.na(perturbation2) | perturbation2 == "")

  # TODO: COMBINE THE FILTER STATEMENTS INTO SINGLE FILTER
  # get wildtypes
  wt_induction_set = condition_fltr_metadata %>%
    filter(genotype1 == "CNAG_00000", strain == 'TDY451')

  # filter for genotypes for those not in the grant -- include doubles
  perturbed_induction_set = condition_fltr_metadata %>%
    filter(genotype1 != "CNAG_00000",
           !genotype1 %in% grant_df$genotype1 |
             (!is.na(genotype2) & !genotype2 %in% grant_df$genotype1))

  # put the wt and filtered genotypes together
  set = bind_rows(wt_induction_set, perturbed_induction_set)


  out = blrs[,colData(blrs)$fastqFileNumber %in% set$fastqFileNumber]

  filterWtByExperimentalLibdate_90min(out)

}

#'
#' Samples with genotypes not in the grant_df
#' @description return a set with only those genotypes (including doubles) not in the grant_df.
#'
#'
#' @import magrittr
#' @importFrom dplyr mutate_if mutate filter bind_rows
#' @importFrom stringr str_remove str_detect
#' @importFrom tidyr replace_na
#'
#' @param blrs a brentlabRnaSeqSet object
#'
#' @return the set metadata
#'
#' @export
createNinetyMinuteInduction_all = function(blrs){

  condition_fltr_metadata = extractColData(blrs) %>%
    filter(medium %in% c("DMEM"),
           temperature %in% c(37),
           atmosphere %in% c("CO2"),
           treatment %in% c("noTreatment"),
           otherConditions %in% c("noOtherConditions"),
           pH %in% c("noPH"),
           timePoint %in% c(90),
           purpose =="fullRNASeq",
           !is.na(fastqFileName),
           str_detect(genotype1, "CNAG"),
           perturbation1 == "deletion" | is.na(perturbation1) | perturbation1 == "",
           perturbation2 == "deletion" | is.na(perturbation2) | perturbation2 == "")

  # TODO: COMBINE THE FILTER STATEMENTS INTO SINGLE FILTER
  # get wildtypes
  wt_set = condition_fltr_metadata %>%
    filter(genotype1 == "CNAG_00000", strain == 'TDY451')

  # filter for genotypes for those not in the grant -- include doubles
  perturbed_set = condition_fltr_metadata %>%
    filter(genotype1 != "CNAG_00000")

  # put the wt and filtered genotypes together
  set = bind_rows(wt_set, perturbed_set)

  out = blrs[,colData(blrs)$fastqFileNumber %in% set$fastqFileNumber]

  filterWtByExperimentalLibdate_90min(out)

}

#' create 90 minute induction set tally
#'
#' @importFrom dplyr group_by tally left_join mutate bind_rows coalesce
#' @importFrom tidyr replace_na
#'
#' @param induction_meta_qual the metadata of the entire set, unfiltered
#' @param sorted_passing_induction_meta_qual metadata (with quality columns)
#'   filtered for manual/auto status
#' @param iqr_fltr_rle_summary sorted_passing_meta_qual filtered for IQR
#' @param grant_df the definition of the 90minuteInduction set. This object is
#'   available in the brentlabRnaSeqTools package
#'
#' @export
createInductionSetTally = function(induction_meta_qual,
                                   sorted_passing_induction_meta_qual,
                                   iqr_fltr_rle_summary,
                                   grant_df){

  # add genotype column -- this is remnant of old system,
  # but kept b/c it might help with concat double KO
  induction_meta_qual$genotype = induction_meta_qual$genotype1
  sorted_passing_induction_meta_qual$genotype = sorted_passing_induction_meta_qual$genotype1
  iqr_fltr_rle_summary$genotype = iqr_fltr_rle_summary$genotype1

  induction_samples_genotype_tally = induction_meta_qual %>%
    group_by(genotype) %>%
    tally()
  colnames(induction_samples_genotype_tally)[2] = "complete_set_no_fltr"

  qc1_passing_tally = sorted_passing_induction_meta_qual %>%
    group_by(genotype) %>%
    tally()
  colnames(qc1_passing_tally)[2] = "qc_passing"

  qc1_iqr_passing_tally = iqr_fltr_rle_summary %>%
    group_by(genotype) %>%
    tally()
  colnames(qc1_iqr_passing_tally)[2] = "qc_passing_iqr_filtered"

  # add genotypes from grant that have no replicates in the database at all
  done_grant_strains_df = grant_df %>% filter(strain_status !=2)

  # TODO: CHECK THAT GRANT_DF COLUMN IS genotype1
  no_geno_df = unique(done_grant_strains_df %>%
                        filter(!genotype1 %in%
                                 unique(induction_meta_qual$genotype1)) %>%
                        select(genotype1))

  no_replicate_genotypes = no_geno_df$genotype1

  no_replicate_genotypes_tally_df = tibble(genotype=no_replicate_genotypes,
                                           complete_set_no_fltr = rep(0, length(no_replicate_genotypes)),
                                           qc_passing = rep(0, length(no_replicate_genotypes)),
                                           qc_passing_iqr_filtered = rep(0, length(no_replicate_genotypes))
  )

  genotype_tally_summary = induction_samples_genotype_tally %>%
    # join the passing and iqr filter tables
    left_join(qc1_passing_tally) %>%
    left_join(qc1_iqr_passing_tally) %>%
    # pass values from previous column forward, if there are NA in preceeding columns
    # mutate(qc_passing = coalesce(qc_passing, complete_set_no_fltr)) %>%  # don't do this, because we don't want to pass a value from complete_set_no_fltr to qc_passing if qc_passing is 0
    mutate(qc_passing_iqr_filtered = coalesce(qc_passing_iqr_filtered, qc_passing)) %>%
    bind_rows(no_replicate_genotypes_tally_df)

  genotype_tally_summary = genotype_tally_summary %>%
    dplyr::mutate(complete_set_no_fltr = replace_na(complete_set_no_fltr, 0)) %>%
    dplyr::mutate(qc_passing = replace_na(qc_passing, 0)) %>%
    dplyr::mutate(qc_passing_iqr_filtered = replace_na(qc_passing_iqr_filtered, 0))

  return(genotype_tally_summary)

}

#' 90min replicate by protocol tally
#'
#' tally 90min induction genotype1 replicates by libraryProtocol
#'
#' @import dplyr
#'
#' @param metadata_df sample metadata from the database or colData(blrs)
#'
replicateByProtocolTally_90min = function(metadata_df){

  # tally replicates by genotype1, protocol
  num_replicate_in_single_vs_both_protocol = metadata_df %>%
    group_by(genotype1, libraryProtocol) %>%
    tally() %>%
    group_by(genotype1) %>%
    tally() %>%
    group_by(n) %>%
    tally %>%
    dplyr::rename(number_of_protocol_by_replicate = n) %>%
    dplyr::rename(number_genotypes = nn)

  # genotypes with samples in single protocol
  genotypes_in_single_protocol = metadata_df %>%
    group_by(genotype1, libraryProtocol) %>%
    tally() %>%
    select(-n) %>%
    group_by(genotype1) %>%
    tally() %>%
    filter(n==1) %>%
    pull(genotype1)

  # genotypes with samples in both protocols
  genotypes_in_both_protocol = metadata_df %>%
    group_by(genotype1, libraryProtocol) %>%
    tally() %>%
    select(-n) %>%
    group_by(genotype1) %>%
    tally() %>%
    filter(n==2) %>%
    pull(genotype1)

  # of the genotypes in the single protocol category, how many in eithe protocol
  num_genotypes_in_single_category_by_nums_in_both_protocols = metadata_df %>%
    filter(genotype1 %in% genotypes_in_single_protocol) %>%
    group_by(genotype1, libraryProtocol) %>%
    tally() %>%
    group_by(libraryProtocol) %>%
    tally()

  # genotype1 replicates with less than 4 samples in both protocols
  replicates_split_btwn_protocols_with_less_than_x_replciates_in_both = metadata_df %>%
    filter(genotype1 %in% genotypes_in_both_protocol) %>%
    group_by(genotype1, libraryProtocol) %>%
    tally() %>%
    pivot_wider(genotype1, names_from=libraryProtocol, values_from=n) %>%
    filter((E7420L < 4 & SolexaPrep < 4))

  list(
    num_replicates_in_old_new_protocols = num_replicate_in_single_vs_both_protocol,
    num_replicates_in_single_protocol =
      num_genotypes_in_single_category_by_nums_in_both_protocols,
    replicates_with_less_than_four_in_both_old_or_new =
      replicates_split_btwn_protocols_with_less_than_x_replciates_in_both
  )

}

#'
#' WT filtering helper for 90minInduction sets
#'
#' @description filter WT samples on libdates without perturbed samples
#'
#' @param blrs_90min a brentlabRnaSeqExperiment object. this function is
#'   written with the 90minInduction sets in mind.
#'
#' @note TODO: Move this to method and select function based on set_name
#'
#' @export
filterWtByExperimentalLibdate_90min = function(blrs_90min){

  # get unique dates of WT
  wt_dates = unique(pull(filter(
    extractColData(blrs_90min), genotype1 == "CNAG_00000"),
    libraryDate))

  # get unique dates of perturbed samples
  perturbed_dates = unique(pull(filter(
    extractColData(blrs_90min), genotype1 != "CNAG_00000"),
    libraryDate))

  # find the wt dates which are not shared by perturbed samples
  wt_drop_dates = setdiff(wt_dates, perturbed_dates)

  # return a filtered object without those WT samples on dates without
  # any perturbed samples
  blrs_90min[,!colData(blrs_90min)$libraryDate %in% wt_drop_dates]
}
