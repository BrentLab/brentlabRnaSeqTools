
#'
#' set global variables, mostly as a hack to make the tidyverse data binding pass the R CMD Check
#'
#' @import tibble
#'
database_colnames = c("bioSampleNumber","harvestDate","harvester","experimentDesign","experimentObservations",
                      "bioSampleObservations","baseStrain","strain","genotype1","perturbation1","marker1","genotype2",
                      "perturbation2", "marker2","medium","temperature","atmosphere","treatment","treatmentConc",
                      "treatmentConcUnit","otherConditions","pH","timePoint","replicate","floodmedia","inductionDelay",
                      "rnaSampleNumber","rnaDate","rnaPreparer","rnaSampleObservations","rnaPrepMethod",
                      "roboticRNAPrep","ribosomalBand","ribosomalBandShape","smallRNABands","rin",
                      "s1cDNASampleNumber","s1cDNADate","s1cDNAPreparer","s1Observations","polyAIsolationProtocol",
                      "s1cDNAProtocol","roboticS1Prep","s1PrimerSeq","s2cDNASampleNumber","s2cDNADate",
                      "s2cDNAPreparer","s2Observations","s2cDNAProtocol","roboticS2Prep","PooledSecondStrand",
                      "librarySampleNumber","libraryDate","libraryPreparer","libraryObservations","index1Sequence",
                      "libraryProtocol","index1Name","index2Name","index2Sequence","roboticLibraryPrep",
                      "fastqFileNumber","fastqFileName","fastqObservations","runNumber","purpose","laneNumber",
                      "sequencerModel","flowcellType","tapestationConc","volumePooled","readsObtained",
                      "qualityAssessmentNumber","manualAudit","manualStatus","librarySize","effectiveLibrarySize",
                      "uniqueAlignment","effectiveUniqueAlignment","effectiveUniqueAlignmentPercent","multiMap",
                      "proteinCodingTotal","proteinCodingTotalPercent","noMap","proteinCodingCounted",
                      "proteinCodingCountedPercent","homopolymerFilter","multiMapPercent","readLengthFilter",
                      "ambiguousFeaturePercent","noFeature","noFeaturePercent","ambiguous","intergenicCoverage",
                      "tooLowAqual","notAlignedTotalPercent","notAligned","noMapPercent","alignmentNotUnique",
                      "homopolyFilterPercent","readLengthFilterPercent","tooLowAqualPercent","rRnaPercent",
                      "nctrRnaPercent","autoAudit","autoStatus","autoStatusDecomp","genotype1Coverage",
                      "genotype1Log2cpm","genotype2Coverage","genotype2Log2cpm","overexpressionFOW","natCoverage",
                      "natLog2cpm","g418Coverage","g418Log2cpm","interquartileRange","totaleDeduplicatedPercentage")

# remove these when the "old qc colnames" no longer matter (after the old rnaseq_pipeline package is taken down)
old_qc_colnames = c("LIBRARY_SIZE","EFFECTIVE_LIBRARY_SIZE","EFFECTIVE_UNIQUE_ALIGNMENT",
                    "EFFECTIVE_UNIQUE_ALIGNMENT_PERCENT","MULTI_MAP_PERCENT","PROTEIN_CODING_TOTAL",
                    "PROTEIN_CODING_TOTAL_PERCENT","PROTEIN_CODING_COUNTED","PROTEIN_CODING_COUNTED_PERCENT",
                    "AMBIGUOUS_FEATURE_PERCENT","NO_FEATURE_PERCENT","INTERGENIC_COVERAGE",
                    "NOT_ALIGNED_TOTAL_PERCENT","GENOTYPE1_COVERAGE","GENOTYPE1_LOG2CPM",
                    "GENOTYPE2_COVERAGE","GENOTYPE2_LOG2CPM","OVEREXPRESSION_FOW","NAT_COVERAGE",
                    "NAT_LOG2CPM","G418_COVERAGE",'G418_LOG2CPM',"NO_MAP_PERCENT","HOMOPOLY_FILTER_PERCENT",
                    "READ_LENGTH_FILTER_PERCENT","TOO_LOW_AQUAL_PERCENT","rRNA_PERCENT","nctrRNA_PERCENT",
                    "STATUS","AUTO_AUDIT","STATUS_DECOMP", "EFFECTIVEUNIQUE_ALIGNMENT")

# TODO all global variables should be capitalized from now on -- uncapitalized variables should be made capital
misc = c("strain_status", "qc_passing_iqr_filtered",
         "qc_passing", "complete_set_no_fltr", "replicate_tally", "model_params",
         "norm_count", ".x", "quant_type", ".", "Bam_file", "fastq_1", "fastq_2", "strandedness", "quantile",
         "IQR", "Sample", "Expression\ Profiling\ Efficiency", "Exonic\ Rate", "Intronic\ Rate", "Intergenic\ Rate",
         "Intragenic\ Rate", "rRNA\ Rate", "Estimated\ Library\ Complexity", "expression_profiling_efficiency",
         "exonic_rate", "intronic_rate", "intergenic_rate", "intragenic_rate", "rRna_rate", "estimated_library_complexity",
         "V1", "V2", "destination", "replicate_group","genotype", "colData<-", 'nn', 'pivot_wider', 'E7420L', 'SolexaPrep',
         'RLE', 'replicate_group_counts', 'feature', 'colname', 'assay<-',
         'unique_identifier', 'replicate_group_rle', 'total', 'head', 'metric',
         'reads', 'percent')

package_data_variables = c("test_metadata",
                           "run_numbers_with_leading_zero",
                           'database_info',
                           "grant_df")

global_variables_vector = c(database_colnames, old_qc_colnames, misc,
                            package_data_variables) #shiny_app_vars

utils::globalVariables(global_variables_vector)
utils::globalVariables(toupper(global_variables_vector))
