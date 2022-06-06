#' Internal binRanges function
#' @inheritParams binRanges
#' @importFrom plyranges mutate
#' @importFrom GenomicRanges trim
#' @noRd
.binRanges <- function(x, binSize, pos) {

    ## Shift, bin, and trim ranges
    x |>
        shiftRanges(pos) |>
        mutate(start = floor(start/binSize)*binSize,
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
#'                ranges = IRanges(start = rep(5000,3),
#'                                 end = rep(6000,3)),
#'                strand = c('+', '-', '*'))
#'
#' gr2 <- gr1 |> promoters(upstream = 2000, downstream = 200)
#'
#' ## Binning the results
#' binRanges(gr1, 'start', 1000)
#' binRanges(gr1, 'end', 1000)
#' binRanges(gr1, 'center', 1000)
#'
#' ## Bin after shifting back to TSS
#' binRanges(gr2, 2000, 1000)
#'
#' @rdname binRanges
#' @export
setMethod("binRanges",
          signature(x = 'GRanges',
                    pos = 'character_OR_numeric_OR_missing',
                    binSize = 'numeric'),
          definition = .binRanges)
