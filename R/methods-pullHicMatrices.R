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


#' Convert GInteractions to Extra short format
#'
#' Extra short format is the minimal information
#' needed to extract Hi-C contacts with `strawr`.
#' See the description for this format here:
#' https://github.com/aidenlab/juicer/wiki/Pre#extra-short-format-dev.
#'
#' It also orders the interactions and chromosomes
#' appropriately. For intrachromosomal interactions, the
#' anchors should be flipped such that start1 < start2.
#' For interchromosomal interactions, columns should be
#' flipped so that chr1 < chr2 (according to the .hic
#' file's internal chromosome map index).
#'
#' @inheritParams pullHicMatrices
#' @importFrom data.table as.data.table
#' @return `data.table` with columns:
#'  "chr1", "start1", "chr2", "start2".
#' @noRd
.GInteractionsToShortFormat <- function(x, file) {

    ## Convert to data.table format
    x <-
        as.data.table(x)[, c("seqnames1", "start1",
                             "seqnames2", "start2")] |>
        `colnames<-`(c("seqnames1", "pos1",
                       "seqnames2", "pos2"))

    return(x)

    ## Get strawr chromosome map index


    ## Interchromosomal: Flip column order
    ## so that chr1 < chr2

    ## Intrachromosmal: Flip column order
    ## so that start1 < start2
    # x[seqnames1 == seqnames2 & start1 > start2,
    #   `:=`(start1=start2, start2=start1)]

}

#' Internal pullHicMatrices
#' @inheritParams pullHicMatrices
#' @return Array of Hi-C submatrices.
#' @noRd
.pullHicMatrices <- function(x, binSize) {

    ## Bin GInteractions if necessary
    x <- .handleBinning(x, binSize)

    ## Convert to short format (seqnames1, pos1, seqnames2, pos2)
    x <- .GInteractionsToShortFormat(x)



}

#' Pull matrices from a `.hic` file
#'
#' @param x GInteractions object.
#' @param binSize Integer (numeric) describing the
#'  resolution (range widths) of the paired data.
#' @param file Character file path to `.hic` file.
#' @return Array of Hi-C submatrices.
#' @noRd
# #' @rdname pullHicMatrices
# #' @export
# setMethod("pullHicMatrices",
#           signature(x = 'GInteractions',
#                     binSize = 'numeric',
#                     file = 'character'),
#           definition = .pullHicMatrices)
