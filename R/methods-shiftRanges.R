#' Internal shiftRanges function
#' @inheritParams shiftRanges
#' @importFrom plyranges mutate
#' @importFrom GenomicRanges trim resize shift
#' @importFrom magrittr `%<>%` `%>%`
#' @importFrom rlang arg_match
#' @noRd
.shiftRanges <- function(x, pos) {

    ## Use GenomicRanges::resize for character pos
    if (is(pos, "character")) {

        pos <- arg_match(pos, values = c('start', 'center', 'end'))
        x <- resize(x = x, width = 1, fix = pos)

    }

    ## Use numeric shifting by strand for numeric pos
    if (is(pos, 'numeric')) {

        ## Convert strand Rle to vector
        sx <- as.vector(strand(x))

        ## Subset for strand if length(pos) > 1
        if (length(pos) > 1) {

            ## Ensure pos matches length of x
            if (length(pos) != length(x)) {
                msg <- c(glue("Mismatched lengths."),
                         'i' = glue("`length(pos)` must equal `length(x)`"),
                         'x' = glue("`length(pos)` is {length(pos)} ",
                                    "but should be {length(x)}."))
                abort(msg)
            }
            pp <- pos[which(sx %in% c('+', '*'))]
            pm <- pos[which(sx == '-')] * -1
        } else {
            pp <- pos
            pm <- pos * -1
        }

        ## Shift '+' or '*' strand
        x[sx %in% c('+', '*')] %<>%
            resize(width = 1, fix = 'start') %>%
            shift(shift = pp)

        ## Shift '-' strand
        x[sx == '-'] %<>%
            resize(width = 1, fix = 'start') %>%
            shift(shift = pm)

    }

    ## Return shifted range
    return(x)

}

#' Flexibly shifting GRanges according to strand
#'
#' @param x GRanges object
#' @param pos Position within anchors to resize the bin.
#'  Can be a character or integer vector of length 1 or
#'  `length(x)` designating the position for each element
#'  in bedpe. Character options are "start", "end" and
#'  "center". Integers are referenced from the start position
#'  for '+' and '*' strands and from the end position
#'  for the '-' strand.
#'
#' @return GRanges object with a single position range
#'  that has been shifted appropriately.
#'
#' @rdname shiftRanges
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
#' ## Shifting anchors by keyword
#' shiftRanges(gr1, 'start')
#' shiftRanges(gr1, 'end')
#' shiftRanges(gr1, 'center')
#'
#' ## Shifting anchors by position
#' shiftRanges(gr1, 100)
#' shiftRanges(gr1, c(100, 200, 300))
#'
#' ## Shifting back to TSS
#' shiftRanges(gr2, 2000)
#'
#' @export
setMethod("shiftRanges",
          signature(x = 'GRanges',
                    pos = 'character_OR_numeric'),
          definition = .shiftRanges)
