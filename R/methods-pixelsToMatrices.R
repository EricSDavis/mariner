#' Internal pixelsToMatrices
#' @inheritParams pixelsToMatrices
#' @importFrom InteractionSet regions anchors GInteractions
#' @importFrom S4Vectors width
#' @importFrom rlang abort
#' @importFrom plyranges mutate
#' @noRd
.pixelsToMatrices <- function(x, buffer) {

    ## Check input types
    .checkTypes(c(buffer="number"))

    ## Get binSize for ranges in x
    binSize <- unique(width(regions(x))) - 1

    ## Ensure x is binned correctly
    if (length(binSize) != 1L) {
        abort(c(
            "All ranges in `x` must be the same width.",
            "i"="Use `binPairs(x)` to set `binSize`."
        ))
    }

    ## Resize each anchor
    a1 <-
        anchors(x, 'first') |>
        mutate(end = start + binSize*buffer,
               start = start - binSize*buffer)

    a2 <-
        anchors(x, 'second') |>
        mutate(end = start + binSize*buffer,
               start = start - binSize*buffer,)

    ## Update x with new ranges
    gi <- GInteractions(a1, a2)
    x <- .updateGInteractions(.Object = x, delegate = gi)

    ## Return resized ranges
    return(x)
}

#'  Expand pixels to submatrices
#'
#'  Pixels are defined as paired-ranges with
#'  starts & ends equal to their `binSize`.
#'  This function takes GInteractions fitting
#'  this description and expands the ranges
#'  such that there is a `buffer` of pixels
#'  around each range.
#'
#'  For example, a buffer of 3 would return a
#'  GInteractions object with 3 pixels surrounding
#'  the original pixel ranges.
#'
#'  After using `pullHicMatrices()`, the result will
#'  return a matrix of row and column dimensions of
#'  buffer*2+1.
#'
#'  Note, this function does not handle out-of-bound
#'  ranges.
#'
#'  @param x GInteractions object.
#'  @param buffer Number (length one numeric vector)
#'   of pixels around the pixels in `x`.
#'
#'  @returns `x` with updated ranges.
#'
#'  @rdname pixelsToMatrices
#'  @export
setMethod("pixelsToMatrices",
          signature(x = 'GInteractions',
                    buffer = 'numeric'),
          definition = .pixelsToMatrices)
