#' Internal pileupDomains function
#' @inheritParams pileupDomains
#' @importFrom rlang inform
#' @importFrom glue glue
#' @importFrom InteractionSet pairdist intrachr
#' @noRd
.pileupDomains <- function(x, files, binSize, buffer,
                       ndim, scale, normalize,
                       FUN, nBlocks, verbose, BPPARAM,
                       blockSize, ...) {
    
    ## Check inputs
    .checkTypes(c(buffer="number"))
    if (buffer < 0) abort("Buffer must be >= 0")
    
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
        
        ## Remove interchromosomal & warn
        if (any(!intrachr(x))) {
            interchr <- which(!intrachr(x))
            warn(c("Interchromosomal interactions are not allowed.",
                   "i"="The following indices in `x` have been removed:",
                   paste0(interchr, collapse=",")))
            x <- x[intrachr(x)]
        }
        
        ## Create new ranges
        df <- data.frame(
            seqnames=seqnames1(x),
            start=start1(x),
            end=end2(x)
        )
        gi <- as_ginteractions(cbind(df, df))
    } else {
        gi <- GInteractions(x, x)
    }
    
    ## Calculate expansion factor
    w <- width(first(gi))
    ifelse(buffer == 0, d <- 0, d <- w/(1/buffer))
    
    ## Expand
    gi <- resize(gi, width=w+d*2, fix='center')
    
    ## Remove out of bound ranges & warn
    # TODOs
    n <- length(gi)
    
    ## Extract
    mat <- pullHicMatrices(
        x=gi,
        files=files,
        binSize=binSize,
        blockSize=blockSize,
        ...
    )
    
    ## Regularize
    if (is(mat, "InteractionJaggedArray")) {
        mat <- regularize(
            x=mat,
            ndim=ndim,
            scale=scale,
            nBlocks=nBlocks,
            verbose=verbose
        )
    }
    
    ## Aggregate
    agg <- aggHicMatrices(
        x=mat,
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

#' Pileup Hi-C domains
#' 
#' pileupDomains expands then extracts regions/domains
#' from Hi-C files, regularizes them so they are the
#' same size, then aggregates them into a single
#' matrix. This can be used to perform aggregate
#' TAD analysis.
#' 
#' It may be necessary to adjust the `zrange` in 
#' `plotMatrix()` since the Hi-C diagonal will 
#' dominate the scale.
#' 
#' If interactions are passed to the function, only
#' intrachromosomal ranges are maintained.
#' 
#' Using small `binSize` values with large ranges
#' may lead to pulling very large sections of a Hi-C map
#' that can crash your R session. If this happens try
#' increasing the `binSize` and `nBlocks` parameters,
#' while lower the `blockSize` parameter.
#' 
#' @param x GRanges or GInteractions object containing
#'  the TADs or Loops to be aggregated.
#' @param files Character file paths to `.hic` files.
#' @param binSize Integer (numeric) describing the
#'  resolution (range widths) of the paired data. Note
#'  that small values for this argument may lead to 
#'  R session crashes.
#' @param buffer Fraction (length one numeric vector)
#'  pair-distance to expand around the resulting
#'  range.
#' @param ndim Numeric vector of length two describing
#'  the new dimensions of the output matrices.
#' @param scale Boolean (TRUE/FALSE) indicating whether
#'  the values in the new matrices should be scaled to the
#'  total signal in each matrix.
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
#' pileupDomains(x=loops, files=hicFile, binSize=50e3, buffer=0.25)
#' 
#' @rdname pileupDomains
#' @export
setMethod("pileupDomains",
          signature(x='GRanges_OR_GInteractions',
                    files='character',
                    binSize='numeric'),
          definition=.pileupDomains)
