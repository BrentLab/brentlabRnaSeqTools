# library(foreach)
# library(iterators)
# library(rtracklayer)
# library(brentlabRnaSeqTools)
# library(tidyverse)
#
# meta = getMetadata(
#   database_info$kn99$db_host,
#   database_info$kn99$db_name,
#   Sys.getenv("db_username"),
#   Sys.getenv("db_password")
# )
#
#
#
# run_df = meta %>%
#   filter(runNumber == 5500)
#
# novoalignPipelineQC = function(meta_df,
#                                pipeline_output_dirpath,
#                                gff_path,
#                                markers = c("NAT", "G418"),
#                                bam_suffix = "_sorted_aligned_reads_with_annote.bam",
#                                novolog_suffix = "_novoalign.log",
#                                count_suffix = '_read_count.tsv',
#                                num_nodes = 10){
#
#   message(paste0("WARNING: this function will not work on a windows machine ",
#           "due to the parallelization method currently implemented. ",
#           "Additionally, threads/cpus must be on the same machine ",
#           "(eg on slurm nodes_per_task=1, cpus_per_task=8)"))
#
#   # extract all perturbed loci in run ----
#   geno1_loci = meta_df %>%
#     filter(genotype1 != 'CNAG_00000') %>%
#     pull(genotype1) %>%
#     as.character() %>%
#     unique()
#
#   geno2_loci = meta_df %>%
#     filter(!is.na(genotype2), genotype2 != '') %>%
#     pull(genotype2) %>%
#     as.character() %>%
#     unique()
#
#   perturbed_loci_df = expand.grid(unique(meta_df$fastqFileNumber),
#                                   c(geno1_loci, geno2_loci)) %>%
#     dplyr::rename(fastqFileNumber = Var1, locus = Var2)
#
#   # create qc df ----
#   qc_df = meta_df %>%
#     dplyr::select(fastqFileNumber, fastqFileName, libraryProtocol) %>%
#     dplyr::rename(strandedness = libraryProtocol) %>%
#     mutate(strandedness =
#              ifelse(as.character(strandedness) == "E7420L",
#                     "reverse", "unstranded")) %>%
#     left_join(perturbed_loci_df)
#
#   # initiate
#   cl = parallel::makeForkCluster(nnodes = num_nodes)
#   doParallel::registerDoParallel(cl)
#
#   perturbed_loci_df_mod =
#     foreach(row=iter(qc_df, by='row'), .combine=rbind) %dopar% {
#       path_list = list(
#         bam = file.path(pipeline_output_dirpath,
#                             "align",
#                             paste0(row$fastqFileName, bam_suffix)),
#         count = file.path(pipeline_output_dirpath,
#                               "count",
#                               paste0(row$fastqFileName, count_suffix)),
#         novolog = file.path(pipeline_output_dirpath,
#                             "logs",
#                             paste0(row$fastqFileName, novolog_suffix))
#       )
#
#
#       marker_granges = list(
#         nat = geneGRanges(gff_path,
#                           "NAT",
#                           feature = 'cds'),
#         g418 = geneGRanges(gff_path,
#                       "G418",
#                       feature = 'cds')
#       )
#
#       # calculate perturbed locus metrics
#       if(!is.na(row$locus)){
#
#         locus_granges = geneGRanges(gff_path,
#                                     str_replace(row$locus, "CNAG", "CKF44"),
#                                                     feature = 'cds')
#
#         row$perturbedCoverage = locusCoverage(path_list$bam,
#                                               locus_granges,
#                                                row$strandedness)
#         row$perturbedLog2cpm = htseq_locusLog2cpm(path_list$count,
#                                                   str_replace(row$locus,
#                                                               "CNAG", "CKF44"))
#       }
#
#       # nat metrics
#       row$natCoverage = locusCoverage(path_list$bam,
#                                       marker_granges$nat,
#                                       row$strandedness)
#       row$natLog2cpm = htseq_locusLog2cpm(path_list$count, "CNAG_NAT")
#
#       # g418 metrics
#       row$g418Coverage = locusCoverage(path_list$bam,
#                                        marker_granges$g418,
#                                        row$strandedness)
#       row$g418Log2cpm = htseq_locusLog2cpm(path_list$count, "CNAG_G418")
#
#       # lib quality metrics
#       row$proteinCodingCounted = htseq_proteinCodingTotal(path_list$count)
#       row$notAlignedTotalPercent = htseq_notAlignedTotalPercent(path_list$novolog)
#       row$libraryComplexity = htseq_libraryComplexity(path_list$count)
#
#       row
#     }
#
#   parallel::stopCluster(cl)
#
#   perturbed_loci_df_mod %>%
#     dplyr::select(-c(fastqFileName, strandedness))
#
#
#
#   # qc_df = meta_df %>%
#   #   dplyr::select(fastqFileNumber,
#   #                 proteinCodingCounted,
#   #                 genotype1Coverage,
#   #                 genotype1Log2cpm,
#   #                 genotype2Log2cpm,
#   #                 genotype2Coverage,
#   #                 natCoverage,
#   #                 natLog2cpm,
#   #                 g418Coverage,
#   #                 g418Log2cpm)
#
#
# }
#
# run_qc = novoalignPipelineQC(run_df, "/mnt/scratch/rnaseq_pipeline/align_count_results/run_5500/rnaseq_pipeline_results/run_5500_samples", Sys.getenv("kn99_stranded_gff_rds"))

# qc_df = run_df %>%
# dplyr::select(fastqFileNumber,
#               genotype1,
#               genotype2)
#
# qc_df = foreach(
#   row = iterators::iter(qc_df, by = 'row'),
#   .combine = 'rbind') %do% {
#     ffn = row$fastqFileNumber
#     genotype1 = as.character(row$genotype1)
#     genotype2 = as.character(row$genotype2)
#
#     genotype1_df = run_qc %>%
#       filter(fastqFileNumber == ffn,
#              locus == genotype1) %>%
#       dplyr::select(fastqFileNumber, perturbedCoverage, perturbedLog2cpm) %>%
#       dplyr::rename(genotype1Coverage = perturbedCoverage,
#                     genotype1Log2cpm = perturbedLog2cpm)
#
#     genotype2_df = run_qc %>%
#       filter(fastqFileNumber == ffn,
#              locus == genotype2) %>%
#       dplyr::select(fastqFileNumber, perturbedCoverage, perturbedLog2cpm) %>%
#       dplyr::rename(genotype2Coverage = perturbedCoverage,
#                     genotype2Log2cpm = perturbedLog2cpm)
#
#     marker_df = run_qc %>%
#       filter(fastqFileNumber == ffn) %>%
#       dplyr::select(fastqFileNumber, natCoverage, natLog2cpm, g418Coverage, g418Log2cpm) %>%
#       distinct(fastqFileNumber, .keep_all = TRUE)
#
#     qc_select = marker_df %>%
#       left_join(genotype1_df) %>%
#       left_join(genotype2_df)
#
#     row %>%
#       as_tibble() %>%
#       left_join(qc_select)
#   }
#
# {"fastqFileNumber":4054,
#   "natCoverage":0,
#   "natLog2cpm":-1.33,
#   "g418Coverage":0,
#   "g418Log2cpm":-1.33},
# {"fastqFileNumber":4055,
#   "natCoverage":0,
#   "natLog2cpm":-1.28,
#   "g418Coverage":0,
#   "g418Log2cpm":-1.28},
# {"fastqFileNumber":4056,
#   "natCoverage":0,
#   "natLog2cpm":-1.12,
#   "g41a8Coverage":0,
#   "g418Log2cpm":-1.12}
