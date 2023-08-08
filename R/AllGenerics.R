#' @rdname createExperimentSet
#' @param ... other arguments
#' @export
setGeneric("createExperimentSet",
           function(x, ...) standardGeneric("createExperimentSet"))

#' @rdname qaFilter
#' @param ... other arguments
#' @export
setGeneric("qaFilter",
           function(x, ...) standardGeneric("qaFilter"))

#' @rdname estimateSizeFactorsByProtocol
#' @param ... other arguments
#' @export
setGeneric("estimateSizeFactorsByProtocol",
           function(x, ...) standardGeneric("estimateSizeFactorsByProtocol"))

#' @rdname replicateByProtocolTally
#' @param ... other arguments
#' @export
setGeneric("replicateByProtocolTally",
           function(x, ...) standardGeneric("replicateByProtocolTally"))

#' @rdname splitProtocolGroups
#' @param ... other arguments
#' @export
setGeneric("splitProtocolGroups",
           function(x, ...) standardGeneric("splitProtocolGroups"))

#' @rdname test_train_partition
#' @param ... other arguments
#' @export
setGeneric("test_train_partition",
           function(x, ...) standardGeneric("test_train_partition"))

#' @rdname extractColData
#' @param ... other arguments
#' @export
setGeneric("extractColData",
           function(x, ...) standardGeneric("extractColData"))

#' @rdname extractDesignMatrix
#' @param ... other arguments
#' @export
setGeneric("extractDesignMatrix",
           function(x, ...) standardGeneric("extractDesignMatrix"))

#' @rdname coerceToDds
#' @param ... other arguments
#' @export
setGeneric("coerceToDds",
           function(x, ...) standardGeneric("coerceToDds"))

#' create IGV viewer batch script (single range, as many tracks as you like)
#'
#' @description this will create a batch script that may be run with the igv
#'   browser. An example command is given in the examples below. See the extern
#'   data for an example of what one of these scripts looks like.
#'
#' @rdname igv_script
#'
#' @importFrom stringr str_replace
#' @importFrom purrr map
#'
#' @param bam_list list of bamfiles -- multiple files means multiple tracks
#' @param output_dir directory to which igv snapshot images will be deposited
#'   by igv
#' @param image_basename basename of the igv snapshot files,
#'    eg ckf44_12345_run2. The output_dir and additional locus details, etc
#'    will be appended to this
#' @param locus the locus you wish to visualize. This will be used to filter
#'   the granges attribute of the object
#' @param exit_browser Default is TRUE, which means that the batchscript will be
#'   written in a way that the browser will take an image, save it, and exit. Set
#'   to false and this script may be used to launch an interactive session.
#' @param maxPanelHeight default 500. height of the IGV window
#'
#' @return None. This prints a batch script to the output dir
#'
#' @note run on the command line like so xvfb-run --auto-servernum igv.sh -b script.bat
#'
#' @export
setGeneric("igv_script",
           function(x,locus,output_dir = getwd(),
                    image_basename = '',
                    exit_browser = TRUE,
                    maxPanelHeight = 500,...) standardGeneric("igv_script"))
































