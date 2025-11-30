#' Internal snap to correct bin
#'
#' Internal helper function to return
#' the starts & ends of ranges snapped
#' to the correct `binSize`.
#'
#' @param start Integer < `end`
#' @param end Integer > `start`
#' @param binSize Integer for size of bins.
#' @importFrom data.table data.table `:=`
#'
#' @returns List of snapped starts and ends.
#'
#' @noRd
.snap <- function(start, end, binSize, tolerance = 1e-08) {

    ## Suppress NSE notes in R CMD check
    coversBins <- coversOneBin <- crossMid <-
        equalSides <- leftSide <- rightSide <- NULL

    ## Get the starts and ends in terms of binSize
    s <- start/binSize
    e <- end/binSize

    ## Calculate the fraction covered by start and end
    sf <- 1 - s %% 1
    ef <- e %% 1

    ## Calculate
    binsCovered <- floor(e)-floor(s)

    ## Perform comparisons
    dt <-
        data.table(s = s,
                   e = e,
                   coversBins = binsCovered > 0,
                   coversOneBin = binsCovered == 1,
                   crossMid = s %% 1 <= 0.5 & e %% 1 >= 0.5,
                   equalSides = abs(sf - ef) <= tolerance,
                   leftSide = sf > ef,
                   rightSide = sf < ef)

    ## Update values
    ## Expand to the whole bin if
    ## the range doesn't cross a bin boundary
    dt[(!coversBins),
       c("s", "e") := .(floor(s), ceiling(e))]

    ## Span multiple bins if range
    ## crosses the bin midpoints
    dt[(coversBins & crossMid),
       c("s", "e") := .(floor(s), ceiling(e))]

    ## If range covers one bin nearly equally
    ## then push to lower bin
    dt[(coversOneBin & !crossMid & equalSides),
       c("s", "e") := .(floor(s), floor(e))]

    ## If range covers many bins nearly equally
    ## then round to nearest bin
    dt[(coversBins & !crossMid & equalSides),
       c("s", "e") := .(round(s), round(e))]

    ## If range is left-sided, push to lower bin
    dt[(coversBins & !crossMid & !equalSides & leftSide),
       c("s", "e") := .(floor(s), floor(e))]

    ## If range is right-sided, push to upper bin
    dt[(coversBins & !crossMid & !equalSides & !leftSide & rightSide),
       c("s", "e") := .(ceiling(s), ceiling(e))]

    ## Expand back to bin coordinates
    list(start = dt$s*binSize, end = dt$e*binSize)

}

#' Internal snapRangesToBins function
#' @inheritParams snapRangesToBins
#' @importFrom dplyr mutate
#' @importFrom GenomicRanges trim
#' @return GRanges object snapped to the nearest `binSize`.
#' @noRd
.snapRangesToBins <- function(x, binSize) {

    ## Suppress R CMD CHECK NOTE
    start <- end <- NULL

    ## Check binSize
    if (binSize == 0) abort("`binSize` must be > 0")

    ## Snap ranges to nearest bin
    snapped <- .snap(start(x), end(x), binSize)

    ## Update object and trim excess
    x |>
        dplyr::mutate(start = snapped$start,
                      end = snapped$end) |>
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
#' library(GenomicRanges)
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
#' library(InteractionSet)
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
