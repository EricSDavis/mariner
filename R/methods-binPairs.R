#' Internal binPairs function
#' @inheritParams binPairs
#' @importFrom data.table data.table
#' @importFrom S4Vectors DataFrame mcols
#' @importFrom InteractionSet anchors
#' @noRd
.binPairs <- function(x, binSize, pos1, pos2) {

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
    a1 <- binRanges(x = a1, binSize = binSize, pos = pos1)
    a2 <- binRanges(x = a2, binSize = binSize, pos = pos2)

    ## Reconstruct binned GInteractions object
    gi <- GInteractions(a1, a2)

    ## Update x with new regions
    x <- .updateGInteractions(x, delegate = gi)

    ## Add back mcols
    mcols(x) <- metadataColumns

    ## Return binned object
    return(x)

}

#' Flexibly bin paired ranges
#'
#' Paired range objects (like `GInteractions`
#' or BEDPE-formatted `data.frame`-like objects)
#' can be binned separately for each set of
#' ranges.
#'
#' @param x `GInteractions` or `data.frame`-like
#'  object with paired interactions.
#' @param binSize Integer (numeric) describing
#'  the new size of each range.
#' @param pos1,pos2 Position within anchors
#'  to resize the bin. Can be a character or
#'  integer vector of length 1 or `length(x)`
#'  designating the position for each element
#'  in `x`. Character options are "start", "end"
#'  and "center". Integers are referenced from
#'  the start position for '+' and '*' strands
#'  and from the end position for the '-' strand.
#' @param ... Additional arguments.
#'
#' @return GInteractions object binned to `binSize`
#'  by `pos1` and `pos2`.
#'
#' @examples
#' ## Construct GInteractions
#' library(InteractionSet)
#' gi1 <-
#'     data.frame(chr1 = "chr1", x1 = 10000, x2 = 20000,
#'                chr2 = "chr1", y1 = 30000, y2 = 40000) |>
#'     as_ginteractions()
#'
#' ## Assign each range to 20-kb bins from the start positions
#' binPairs(x = gi1,
#'          binSize = 20000,
#'          pos1 = 'start',
#'          pos2 = 'start')
#'
#' @rdname binPairs
#' @export
setMethod("binPairs",
          signature(x = 'DF_OR_df_OR_dt',
                    binSize = 'numeric',
                    pos1 = 'character_OR_numeric_OR_missing',
                    pos2 = 'character_OR_numeric_OR_missing'),
          definition = .binPairs)

#' @rdname binPairs
#' @export
setMethod("binPairs",
          signature(x = 'GInteractions',
                    binSize = 'numeric',
                    pos1 = 'character_OR_numeric_OR_missing',
                    pos2 = 'character_OR_numeric_OR_missing'),
          definition = .binPairs)
