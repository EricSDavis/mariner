#' Ensure pairs are binned and warn if not
#' @inheritParams pullHicMatrices
#' @importFrom rlang inform
#' @importFrom glue glue
#' @return GInteractions object that has been
#'  binned (or original object if already binned).
#' @noRd
.handleBinning <- function(x, binSize) {

    ## Check for binned GInteractions
    binned <- .checkBinnedPairs(x, binSize)

    ## Inform user and bin GInteractions
    if (!binned) {
        x <- binPairs(x, binSize)
        msg <- c("Pairs are not binned to `binSize`.",
                 'i' = glue("Binning with binSize={binSize}, ",
                            "pos1='center' and pos2='center'."),
                 'i' = glue("Use `binPairs()` for different binning."))
        inform(msg)
    }

    ## Return appropriately binned object
    return(x)
}

#' Internal pullHicMatrices
#' @inheritParams pullHicMatrices
#' @return Array of Hi-C submatrices.
#' @noRd
.pullHicMatrices <- function(x, binSize) {

    ## Input checking and processing ----

    ## Bin GInteractions if necessary
    x <- .handleBinning(x, binSize)





}

#' Pull matrices from a `.hic` file
#'
#' @param x GInteractions object.
#' @param binSize Integer (numeric) describing the
#'  resolution (range widths) of the paired data.
#' @return Array of Hi-C submatrices.
#' @noRd
# #' @rdname pullHicMatrices
# #' @export
# setMethod("pullHicMatrices",
#           signature(x = 'GInteractions',
#                     binSize = 'numeric'),
#           definition = .pullHicMatrices)
