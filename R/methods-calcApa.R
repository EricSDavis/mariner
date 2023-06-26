
#' Internal calcApa function
#' @inheritParams calcApa
#' @noRd
.calcApa <- function(x, files, binSize, buffer,
                     removeShort, minPairDist, normalize,
                     FUN, nBlocks, verbose, BPPARAM, ...) {
    
    ## Filter out short pairs
    if (removeShort) {
        x <- removeShortPairs(x, padding=minPairDist)
    }
    n <- length(x)
    
    ## Convert to matrix regions
    regions <- pixelsToMatrices(x, buffer)
    
    ## Extract counts
    mats <- pullHicMatrices(x=regions, files=files, binSize=binSize, ...)
    
    ## Aggregate
    agg <- aggHicMatrices(
        x=mats,
        FUN=FUN,
        nBlocks=nBlocks,
        verbose=verbose,
        BPPARAM=BPPARAM
    )
    
    ## Normalize to counts per loop
    if (normalize) {
        agg <- agg/n
    }
    return(agg)
}

#' Conduct aggregate peak analysis
#' 
#' calcApa optionally removes short interactions
#' that intersect the diagonal before converting
#' the interactions into square regions that are
#' extracted from Hi-C files then aggregated into
#' a single matrix.
#' 
#' @param x GInteractions object containing interactions
#'  to extract from Hi-C files.
#' @param files Character file paths to `.hic` files.
#' @param binSize Integer (numeric) describing the
#'  resolution (range widths) of the paired data.
#' @param buffer Integer indicating the
#'  buffer size, or number of pixels
#' @param removeShort Boolean, whether to remove
#'  short pairs (Default) or not.
#' @param minPairDist Pairs with a distance less than
#'  or equal to this value will be filtered out.
#' @param normalize Boolean, whether to normalize
#'  the aggregated values to the number of interactions
#'  (after filtering out short pairs - if applicable).
#' @param FUN Function to use for aggregating.
#' @param nBlocks Number of blocks for block-processing
#'  arrays. Default is 5. Increase this for large
#'  datasets. To read and process all data at once, set
#'  this value to 1.
#' @param verbose Boolean (TRUE or FALSE) describing
#'  whether to report block-processing progress.
#' @param BPPARAM Parallelization params (passed to
#'  `BiocParallel::bplapply()`). Default is the result
#'  of `BiocParallel::bpparams()`. Parallel processing
#'  is not available when `by=interactions`.
#' @param ... Additional arguments passed to 
#'  `pullHicMatrices()`.
#' 
#' @returns A DelayedMatrix of aggregated
#'  counts.
#'  
#' @examples
#' ## Load marinerData
#' if (!require("marinerData", quietly = TRUE))
#'     BiocManager::install("marinerData")
#' 
#' ## Read .hic file paths
#' hicFile <- marinerData::LEUK_HEK_PJA30_inter_30.hic()
#' names(hicFile) <- "WT"
#' 
#' ## Loops
#' loops <-
#'     WT_5kbLoops.txt() |>
#'     setNames("WT") |>
#'     read.table(header=TRUE, nrows=1000) |>
#'     as_ginteractions(keep.extra.columns=FALSE) |>
#'     binPairs(binSize=5e3)
#' 
#' ## Removes the "chr" prefix for compatibility
#' ## with the preprocessed hic files
#' GenomeInfoDb::seqlevelsStyle(loops) <- 'ENSEMBL'
#' 
#' ## APA
#' mat <- calcApa(
#'     x=loops,
#'     files=hicFile,
#'     binSize=5e3,
#'     minLoopSize=50e3,
#'     normalize=FALSE
#' )
#' 
#' 
#' @rdname calcApa
#' @export
setMethod("calcApa",
          signature(x='GInteractions',
                    files='character',
                    binSize='numeric'),
          definition=.calcApa)