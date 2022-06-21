## Accessors -------------------------------------------------------------------

#' Get selectionMethod from MergedGInteractions object
#'
#' @param x MergedGInteractions object.
#' @param ... Additional arguments.
#'
#' @return A character vector describing which selection
#'  method was used for merging.
#'
#' @examples
#' ## Reference BEDPE files (loops called with SIP)
#' bedpeFiles <-
#'     system.file("extdata", package = "mariner") |>
#'     list.files(pattern = "Loops.txt", full.names = TRUE)
#'
#' x <- mergePairs(x = bedpeFiles,
#'                 binSize = 5e03,
#'                 radius = 2,
#'                 column = "APScoreAvg")
#' selectionMethod(x)
#'
#' @rdname selectionMethod
#' @export
setMethod("selectionMethod", "MergedGInteractions", function(x, ...) {
    x@selectionMethod
})
