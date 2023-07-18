#' Internal fxn to add row/colnames to matrix
#' @param x DelayedArray object.
#' @param rows rownames assay from InteractionArray
#' @param cols colnames assay from InteractionArray
#' @noRd
.mat_with_dimnames <- function(x, rows, cols) {
    if (is.vector(x)) {
        x <- matrix(x, length(rows), length(cols))
        x <- DelayedArray(x)
    }
    dimnames(x) <- list(rows, cols)
    capture.output(show(x))[-1]
}

#' Visualize dimnames for a CountMatrix object
#'
#' @param object A CountMatrix object.
#' @returns A DelayedArray object with dimnames
#'  for the first two dimensions.
#' @examples
#' ## Load marinerData
#' if (!require("marinerData", quietly = TRUE))
#'     BiocManager::install("marinerData")
#'
#' ## Read .hic file paths
#' hicFiles <- c(
#'     marinerData::LEUK_HEK_PJA27_inter_30.hic(),
#'     marinerData::LEUK_HEK_PJA30_inter_30.hic()
#' )
#' names(hicFiles) <- c("FS", "WT")
#'
#' ## Read in loop pixels as GInteractions object
#' pixels <-
#'     WT_5kbLoops.txt() |>
#'     setNames("WT") |>
#'     read.table(header=TRUE) |>
#'     as_ginteractions(keep.extra.columns=FALSE) |>
#'     assignToBins(binSize=100e3)
#'
#' ## Removes the "chr" prefix for compatibility
#' ## with the preprocessed hic files
#' GenomeInfoDb::seqlevelsStyle(pixels) <- 'ENSEMBL'
#'
#' ## Expand pixels to regions for pulling
#' ## Hi-C submatrices
#' regions <- pixelsToMatrices(x=pixels, buffer=5)
#'
#' ## Extract 11x11 count matrices from the
#' ## first 100 regions and 2 Hi-C files
#' iarr <- pullHicMatrices(x=regions[1:100],
#'                         files=hicFiles,
#'                         binSize=100e3)
#'
#' ## Display the start bin of each
#' ## interaction in the count
#' ## matrices
#' counts(iarr, showDimnames=TRUE)
#' @rdname CountMatrix-class
#' @export
setMethod("show", "CountMatrix", function(object) {

    ## Pull out seed and InteractionArray object
    cnts <- object@seed
    object <- object@object

    ## Row/colnames
    if (!all(c("rownames", "colnames") %in%
             names(assays(object)))) {
        abort(c("Dimnames not available for this object.",
                "*"="Try again with `showDimnames=FALSE`."))
    }
    rows <- assay(object, 'rownames')
    cols <- assay(object, 'colnames')

    ## Get the header and dimensions
    fullHeader <- capture.output(show(cnts))[1]
    dims <- dim(cnts)
    dn <- dimnames(cnts)

    ## Construct first & last matrix header
    firstHeader <-
        paste0(c(",,"),
               ifelse(is.null(dn[[3]]), 1, dn[[1]][1]), ",",
               ifelse(is.null(dn[[4]]), 1, dn[[4]][1]))

    lastHeader <-
        paste0(c(",,"),
               ifelse(is.null(dn[[3]]),
                      dims[3], dn[[1]][dims[3]]), ",",
               ifelse(is.null(dn[[4]]),
                      dims[4], dn[[4]][dims[4]]))

    ## Construct new object to print
    ans <- fullHeader
    ans <- c(ans, firstHeader)
    ans <- c(ans, .mat_with_dimnames(cnts[,,1,1],
                                     rows[1,1,],
                                     cols[1,1,]))
    if (dims[3] != 1L | dims[4] != 1L) {
        ans <- c(ans, "\n...\n")
        ans <- c(ans, lastHeader)
        ans <- c(ans, .mat_with_dimnames(cnts[,,dims[3],dims[4]],
                                         rows[dims[3],dims[4],],
                                         cols[dims[3],dims[4],]))
    } else {
        dimnames(cnts)[c(1,2)] <- list(rows[1,1,], cols[1,1,])
    }
    cat(ans, sep="\n")
})
