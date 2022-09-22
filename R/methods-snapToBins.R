#' Internal snapRangesToBins function
#' @inheritParams snapRangesToBins
#' @importFrom plyranges mutate
#' @importFrom GenomicRanges trim
#' @return GRanges object snapped to the nearest `binSize`.
#' @noRd
.snapRangesToBins <- function(x, binSize) {

    ## Suppress R CMD CHECK NOTE
    start <- end <- NULL

    ## Check binSize
    if (binSize == 0) abort("`binSize` must be > 0")

    ## Snap ranges to nearest bin
    x |>
        mutate(start = round(start/binSize)*binSize,
               end = round(end/binSize)*binSize) |>
        mutate(start = ifelse(width <= binSize,
                              floor(start/binSize)*binSize,
                              start),
               end = ifelse(width <= binSize,
                            ceiling(end/binSize)*binSize,
                            end)) |>
        trim() |>
        suppressWarnings()

}

#' Snap GRanges or GInteractions to nearest bins
#'
#' @param x `GRanges` object.
#' @param binSize Integer (numeric) describing
#'  the new size of each range.
#'
#' @return GRanges object snapped to the nearest `binSize`.
#' @examples
#' ## Example GRanges object
#' x <- GRanges(seqnames = c("chr1"),
#'              ranges = IRanges(start = c(1, 1, 25, 19, 21),
#'                               end = c(15, 11, 31, 31, 39)))
#'
#' snapToBins(x, binSize = 5)
#' snapToBins(x, binSize = 10)
#' snapToBins(x, binSize = 20)
#'
#' @rdname snapToBins
#' @export
setMethod("snapToBins",
          signature(x = 'GRanges',
                    binSize = 'numeric'),
          definition = .snapRangesToBins)

#' Internal snapPairsToBins function
#' @inheritParams snapPairsToBin
#' @importFrom data.table data.table
#' @importFrom S4Vectors DataFrame mcols
#' @importFrom InteractionSet anchors
#'
#' @return Input object that has been snapped to `binSize`
#'
#' @noRd
.snapPairsToBins <- function(x, binSize) {

    ## Convert x to a GInteractions
    if (is(x, 'data.frame') | is(x, 'DFrame')) {
        x <- makeGInteractionsFromDataFrame(x)
    }

    ## Capture mcols
    metadataColumns <- mcols(x)

    ## Extract anchors
    a1 <- anchors(x, type = 'first')
    a2 <- anchors(x, type = 'second')

    ## Bin anchors
    a1 <- .snapRangesToBins(x = a1, binSize = binSize)
    a2 <- .snapRangesToBins(x = a2, binSize = binSize)

    ## Reconstruct binned GInteractions object
    gi <- GInteractions(a1, a2)

    ## Update x with new regions
    x <- .updateGInteractions(x, delegate = gi)

    ## Add back mcols
    mcols(x) <- metadataColumns

    ## Return binned object
    return(x)

}

#' Snap paired-objects to nearest bins
#'
#' @param x `GInteractions` object.
#' @param binSize Integer (numeric) describing
#'  the new size of each range.
#'
#' @return Input object snapped to the nearest `binSize`.
#' @examples
#' ## Sample GInteractions object
#' x <- GInteractions(anchor1 = c(GRanges("chr1:1-15"),
#'                                GRanges("chr1:1-11")),
#'                    anchor2 = c(GRanges("chr1:25-31"),
#'                                GRanges("chr1:19-31")))
#'
#' snapToBins(x, binSize = 5)
#' snapToBins(x, binSize = 10)
#' snapToBins(x, binSize = 20)
#'
#' @rdname snapToBins
#' @export
setMethod("snapToBins",
          signature(x = 'GInteractions',
                    binSize = 'numeric'),
          definition = .snapPairsToBins)
