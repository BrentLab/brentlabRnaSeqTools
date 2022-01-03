#' create IGV viewer batch script (single range, as many tracks as you like)
#'
#' @description this will create the following:
#'   new                                # batchscript keyword new (new snapshot)
#'   snapshotDirectory IGV_Snapshots    # output directory
#'   maxPanelHeight 500                 # maximum height of igv browser viewer
#'   genome organism.genome             # igv .genome file
#'   load control.bam                   # load alignment files in bam_file_list
#'   load perturbed.bam
#'   goto chr1:1-100                    # load region of interest
#'   snapshot batchfilename_locus1.png  # saves a snapshot of the IGV window to an image file
#'   goto chr10:3-300                   # repeat at another locus
#'   snapshot batchfilename_locus2.png
#'   exit                               # quit session
#'
#' @importFrom stringr str_replace
#' @importFrom purrr map
#'
#' @param bam_list list of bamfiles -- multiple files means multiple tracks
#' @param granges granges describing the range you wish to visualize -- range
#'   must be on one chromosome/contig
#' @param igv_genome .genome file created by IGV tools
#' @param output_dir where to deposit both the script and the browser shots
#' @param output_file_basename this will serve as both the name of the browser
#'   shot after running IGV, and the name o the batchscript itself
#' @param exit_browser Default is TRUE, which means that the batchscript will be
#'   written in a way that the browser will take an image, save it, and exit. Set
#'   to false and this script may be used to launch an interactive session.
#' @param maxPanelHeight default 500. height of the IGV window
#'
#'
#' @return None. This prints a batch script to the output dir
#'
#' @examples
#' \dontrun{
#' # You need to have igv installed on the machine where you will run this. This
#' # is something you can run on your local. Create the browser shot like so:
#'
#' xvfb-run --auto-servernum igv.sh -b script.bat
#'
#' }
#'
#' @export
createIgvBatchscript = function(bam_list,
                                granges,
                                igv_genome,
                                output_dir,
                                output_file_basename,
                                exit_browser = TRUE,
                                maxPanelHeight=500){

  output_file_basename = tools::file_path_sans_ext(output_file_basename)
  output_img_filename = paste0(output_file_basename, ".png")

  load_samples = paste0(map(bam_list, ~sprintf("\tload %s\n", .)), collapse=" ")

  parsed_range = paste(as.character(granges@seqnames),
                       as.character(granges@ranges), sep=":")

  if(parsed_range == character(0)){
    stop("range object cannot be empty")
  }

  batch_script = paste0("new\nsnapshotDirectory %s\nmaxPanelHeight %s\ngenome %s\n",
                        load_samples, "goto %s")

  batch_script = sprintf(batch_script,
                         output_dir,
                         maxPanelHeight,
                         igv_genome,
                         parsed_range)

  if(exit_browser){
    batch_script = paste0(batch_script, "\nsnapshot %s\nexit")
    batch_script = sprintf(batch_script, output_img_filename)
  }

  script_output_dir = file.path(output_dir, "scripts")
  dir.create(script_output_dir, recursive = TRUE)

  output_filename = file.path(script_output_dir,
                              str_replace(output_file_basename, ".png", ".txt"))

  cat(batch_script, file = output_filename)
}
