#' Internal pileupBoundaries function
#' @inheritParams pileupBoundaries
#' @importFrom rlang inform
#' @importFrom glue glue
#' @importFrom InteractionSet reduceRegions regions GInteractions
#' @importFrom IRanges resize
#' @noRd
.pileupBoundaries <- function(x, files, binSize, width, normalize,
                              FUN, nBlocks, verbose, BPPARAM,
                              blockSize, ...) {
    
    ## Check inputs
    .checkTypes(c(width="number"))
    if (width < 0) abort("`width` must be >= 0")
    
    ## Warnings
    if (binSize <= 10e3) {
        inform(c(glue(
            "Caution: setting a small `binSize` ",
            "with large ranges may cause your ",
            "RSession to crash."),
            "i"="Try a larger `binSize`."))
    }
    
    ## Create redundant GInteractions
    if (is(x, "GInteractions")) {
        ## Get unique anchor regions used in object
        a <- reduceRegions(x) |> regions()
        gi <- GInteractions(a, a)
    } else {
        gi <- GInteractions(x, x)
    }
    
    ## Expand regions
    gi <- resize(gi, width=width, fix='center')
    n <- length(gi)
    
    ## Extract
    mat <- pullHicMatrices(
        x=gi,
        files=files,
        binSize=binSize,
        blockSize=blockSize,
        ...
    )
    
    ## Aggregate
    agg <- aggHicMatrices(
        x=mat,
        FUN=FUN,
        nBlocks=nBlocks,
        verbose=verbose,
        BPPARAM=BPPARAM
    )
    
    ## Normalize to counts per loop
    if (normalize) agg <- agg/n
    return(agg)
}

#' Pileup Hi-C contacts around boundary regions
#' 
#' pileupBoundaries expands input loci to the
#' specified `width`, extracts then aggregates them
#' into a single matrix. This can be used to aggregate
#' windows of interactions centered on a set of loci.
#' 
#' It may be necessary to adjust the `zrange` in 
#' `plotMatrix()` since the Hi-C diagonal will 
#' dominate the scale.
#' 
#' Using small `binSize` values with large ranges
#' may lead to pulling very large sections of a Hi-C map
#' that can crash your R session. If this happens try
#' increasing the `binSize` and `nBlocks` parameters,
#' while lower the `blockSize` parameter.
#' 
#' @param x GRanges or GInteractions object containing
#'  the loci to be aggregated. GInteractions will be
#'  split into unique anchors.
#' @param files Character file paths to `.hic` files.
#' @param binSize Integer (numeric) describing the
#'  resolution (range widths) of the paired data. Note
#'  that small values for this argument may lead to 
#'  R session crashes.
#' @param width Number of base pairs to expand the loci
#'  of interest in `x`.
#' @param normalize Boolean, whether to normalize
#'  the aggregated values to the number of interactions.
#' @param FUN Function to use for aggregating.
#' @param nBlocks Number of blocks for block-processing
#'  arrays. Default is 50. Increase this for large
#'  datasets. To read and process all data at once, set
#'  this value to 1.
#' @param verbose Boolean (TRUE or FALSE) describing
#'  whether to report block-processing progress.
#' @param BPPARAM Parallelization params (passed to
#'  `BiocParallel::bplapply()`). Default is the result
#'  of `BiocParallel::bpparams()`. Parallel processing
#'  is not available when `by=interactions`.
#' @param blockSize Number (length one numeric vector)
#'  describing the size in base-pairs to pull from each
#'  `.hic` file. Default is 1e6. For large `.hic` files
#'  `blockSize` can be reduced to conserve the amount of
#'   data read in at a time. Larger `blockSize` values
#'   speed up performance, but use more memory.
#' @param ... Additional arguments passed to 
#'  `pullHicMatrices()`.
#' 
#' @returns A DelayedArray of aggregated
#'  counts.
#' 
#' @examples
#' ## Load marinerData
#' if (!require("marinerData", quietly = TRUE))
#'     BiocManager::install("marinerData")
#' 
#' ## Read .hic file paths
#' hicFile <- marinerData::LEUK_HEK_PJA27_inter_30.hic()
#' names(hicFile) <- "FS"
#' 
#' ## Loops
#' loops <- 
#'     marinerData::FS_5kbLoops.txt() |>
#'     read.table(header=TRUE, nrows=100) |>
#'     as_ginteractions() |>
#'     GenomeInfoDb::`seqlevelsStyle<-`(value='ENSEMBL')
#' 
#' ## Warn about small binSize
#' pileupBoundaries(x=loops, files=hicFile, binSize=50e3)
#' 
#' @rdname pileupBoundaries
#' @export
setMethod("pileupBoundaries",
          signature(x='GRanges_OR_GInteractions',
                    files='character',
                    binSize='numeric'),
          definition=.pileupBoundaries)
