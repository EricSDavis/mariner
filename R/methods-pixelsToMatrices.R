#' Internal pixelsToMatrices
#' @inheritParams pixelsToMatrices
#' @importFrom InteractionSet regions anchors GInteractions
#' @importFrom S4Vectors width mcols "mcols<-"
#' @importFrom rlang abort
#' @importFrom dplyr mutate
#' @noRd
.pixelsToMatrices <- function(x, buffer) {

    ## Check input types
    .checkTypes(c(buffer="number"))

    ## Get binSize for ranges in x &
    ## ensure x is binned correctly
    binSize <- .getBinSize(x)

    ## Resize each anchor
    a1 <-
        anchors(x, 'first') |>
        dplyr::mutate(end = start + binSize*(buffer+1),
                      start = start - binSize*buffer)

    a2 <-
        anchors(x, 'second') |>
        dplyr::mutate(end = start + binSize*(buffer+1),
                      start = start - binSize*buffer)

    ## Form delegate object & preserve mcols
    gi <- GInteractions(a1, a2)
    mcols(gi) <- mcols(x)

    ## Update x with new ranges
    x <- .updateGInteractions(.Object = x, delegate = gi)

    ## Return resized ranges
    return(x)
}

#' Expand pixels to submatrices
#'
#' Pixels are defined as paired-ranges with
#' starts & ends equal to their `binSize`.
#' This function takes GInteractions fitting
#' this description and expands the ranges
#' such that there is a `buffer` of pixels
#' around each range.
#'
#' For example, a buffer of 3 would return a
#' GInteractions object with 3 pixels surrounding
#' the original pixel ranges.
#'
#' After using `pullHicMatrices()`, the result will
#' return a matrix of row and column dimensions of
#' buffer*2+1.
#'
#' Note, this function does not handle out-of-bound
#' ranges.
#'
#' @param x GInteractions object.
#' @param buffer Number (length one numeric vector)
#'  of pixels around the pixels in `x`.
#'
#' @returns `x` with updated ranges.
#'
#' @examples
#' ## Define example 100bp pixel
#' library(InteractionSet)
#' pixel <- GInteractions(
#'     anchor1=GRanges("chr1:500-600"),
#'     anchor2=GRanges("chr1:2000-2100")
#' )
#'
#' ## Expand pixel to matrix with
#' ## 3 pixels surrounding the center
#' ## pixel
#' region <- pixelsToMatrices(x=pixel, buffer=3)
#' region
#'
#' @rdname pixelsToMatrices
#' @export
setMethod("pixelsToMatrices",
          signature(x = 'GInteractions',
                    buffer = 'numeric'),
          definition = .pixelsToMatrices)
