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

#' @rdname createTestTrainSet
#' @param ... other arguments
#' @export
setGeneric("createTestTrainSet",
           function(x, ...) standardGeneric("createTestTrainSet"))

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
