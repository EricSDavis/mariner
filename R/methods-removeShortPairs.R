#' Internal function for removeShortPairs
#' @inheritParams removeShortPairs
#' @importFrom data.table as.data.table
#' @returns A GInteractions object with
#'  the short pairs removed.
#' @noRd
.removeShortPairs <- function(x, padding) {

    ## Suppress NSE notes in R CMD check
    seqnames1 <- seqnames2 <- start1 <-
        start2 <- end1 <- end2 <- NULL

    ## Convert GInteractions to data.table
    dat <- as.data.table(x)

    ## Filter by distance to diagonal with padding
    keep <- dat[
        ## Keep interchromosomal
        seqnames1 != seqnames2 |

            ## Check upper triangular (lower left corner)
            (seqnames1 == seqnames2 &
                 start1 < start2 &
                 start2 > end1 & # ensure not overlapping
                 start2 - (end1 + padding) >= 0) |

            ## Check lower triangular (upper right corner)
            (seqnames1 == seqnames2 &
                 start1 > start2 &
                 start1 > end2 & # ensure not overlapping
                 start1 - (end2 + padding) >= 0),
        which=TRUE
    ]

    x[keep]
}

#' Remove interactions that would cross
#' the Hi-C diagonal or a specified
#' distance from the diagonal.
#'
#' Note this is only applies to
#' intrachromosomal pairs, as pair distance
#' is meaningless for interchromosomal
#' pairs. Therefore, all interchromosomal
#' pairs are kept.
#'
#' @param x A GInteractions object.
#' @param padding Minimum distance away
#'  from the diagonal.
#'
#' @examples
#' ## Example GInteractions object
#' gi <- as_ginteractions(read.table(
#'     text="
#'         seqnames1 start1 end1 seqnames2 start2 end2 keep
#'         chr1 300 400 chr1 300 400 'no'
#'         chr1 100 200 chr1 300 400 'yes'
#'         chr1 300 400 chr1 100 200 'yes'
#'         chr1 300 400 chr2 300 400 'yes'
#'         chr1 250 350 chr1 300 400 'only_with_padding_50'
#'         chr1 300 400 chr1 250 350 'only_with_padding_50'
#'         ",
#'     header=TRUE
#' ))
#'
#' ## Remove pairs that would cross the diagonal
#' removeShortPairs(gi)
#'
#' ## Add 50bp of padding
#' removeShortPairs(gi, padding=50)
#'
#' @returns A GInteractions object with
#'  the short pairs removed.
#' @rdname removeShortPairs
#' @export
setMethod("removeShortPairs",
          signature(x='GInteractions'),
          definition=.removeShortPairs)
