#' Internal binRanges function
#' @inheritParams binRanges
#' @importFrom dplyr mutate
#' @importFrom GenomicRanges trim
#' @return GRanges object binned to `binSize` from `pos`
#' @noRd
.binRanges <- function(x, binSize, pos) {

    ## Suppress R CMD CHECK NOTE
    start <- NULL

    ## Check binSize
    if (binSize == 0) abort("`binSize` must be > 0")

    ## Shift, bin, and trim ranges
    x |>
        shiftRanges(pos) |>
        dplyr::mutate(start = floor(start/binSize)*binSize,
               end = floor(start/binSize)*binSize + binSize) |>
        trim() |>
        suppressWarnings()

}

#' Flexibly bin ranges
#'
#' @param x `GRanges` object
#' @param binSize Integer (numeric) describing
#'  the new size of each range.
#' @param pos Position within range
#'  to resize the bin. Can be a character or
#'  integer vector of length 1 or `length(x)`
#'  designating the position for each element
#'  in `x`. Character options are "start", "end"
#'  and "center". Integers are referenced from
#'  the start position for '+' and '*' strands
#'  and from the end position for the '-' strand.
#'
#' @return `GRanges` object that has been shifted
#'  by `pos` and assigned to bins of `binSize`.
#'
#' @examples
#' library(GenomicRanges)
#'
#' ## Create example GRanges
#' gr1 <- GRanges(seqnames = "chr1",
#'                ranges = IRanges::IRanges(start = rep(5000,3),
#'                                          end = rep(6000,3)),
#'                strand = c('+', '-', '*'))
#'
#' gr2 <- gr1 |> promoters(upstream = 2000, downstream = 200)
#'
#' ## Binning the results
#' binRanges(x = gr1, binSize = 1000, pos = 'start')
#' binRanges(x = gr1, binSize = 1000, pos = 'end')
#' binRanges(x = gr1, binSize = 1000, pos = 'center')
#'
#' ## Bin after shifting back to TSS
#' binRanges(x = gr2, binSize = 1000, pos = 2000)
#'
#' @rdname binRanges
#' @export
setMethod("binRanges",
          signature(x = 'GRanges',
                    pos = 'character_OR_numeric_OR_missing',
                    binSize = 'numeric'),
          definition = .binRanges)
