#!/usr/bin/env Rscript

# this has to be here for now b/c HTCF can't handle shared software. The path
# is correct -- don't change it. Eventually it will be moved into a shared directory
.libPaths("/scratch/mblab/chasem/rnaseq_pipeline/experiments/zev_induction/r_library_for_mpi_deseq")

suppressMessages(library("optparse"))


main = function(args){

  message("loading dependencies")
  suppressMessages(library(DESeq2))
  suppressMessages(library(Rmpi))
  suppressMessages(library(BiocParallel))

  message("resgistering multicore parameters")
  param <- bpstart(SnowParam(mpi.universe.size() -1, "MPI")) # , log = TRUE, threshold = "DEBUG"))
  register(param)

  message("creating output path from input dds object")
  dds_basename = gsub("input|input_|_input$",
                      "",
                      tools::file_path_sans_ext(basename(args$deseq_data_set)))

  output_path = file.path(dirname(args$deseq_data_set),
                          paste0(dds_basename,
                                 "_", format(Sys.time(), "%Y%m%d"),
                                 "_output.rds"))

  message("reading in dds object")
  dds = readRDS(args$deseq_data_set)

  message("running DESeq()")
  deseq_model = DESeq(dds, parallel=TRUE)

  message("saving deseq output")
  saveRDS(deseq_model, output_path)

  message("quitting MPI")
  bpstop(param)
  mpi.quit()

  message(paste0("all done for: ", args$deseq_data_set))

} # end main()

parseArguments <- function() {
  # parse and return cmd line input

  option_list <- list(
    make_option(c('-d', '--deseq_data_set'),
                help='the dds(deseq data set). Output will be the same name in same directory with the word output and the date appended'),
    make_option(c('-l', '--lib_path'),
                help='path to library with deseq, Rmpi and biocparallel, at least'))

  args <- parse_args(OptionParser(option_list=option_list))
  return(args)
} # end parseAarguments

main(parseArguments()) # call main method

# for testing
# input_list = list()
# input_list['deseq_data_set'] = '/home/chase/code/cmatkhan/misc_scripts/deseq_model/data/test_2_counts.csv'
# input_list['output_directory'] = '/home/chase/Desktop/tmp/test_results'
# input_list['name'] = 'deseq_output_test'
# input_list['threads'] = '10'

# main(input_list)
