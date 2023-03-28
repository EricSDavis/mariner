#' Internal for calcLoopEnrichment
#' @inheritParams calcLoopEnrichment
#' @importFrom rlang abort inform
#' @importFrom glue glue
#' @importFrom DelayedArray RegularArrayGrid DelayedArray blockApply
#' @importFrom utils capture.output
#' @importFrom stats median
#' @noRd
.calcLoopEnrichment <- function(x, files, fg, bg, FUN, nBlocks, verbose,
                                BPPARAM, ...) {

    ## Parameter parsing
    if (nBlocks <= 0) abort("`nBlocks` must be > 0.")

    ## Determine resolution from x
    ## and ensure all pixels are same res.
    binSize <- .getBinSize(x)
    if (length(binSize) != 1L) {
        abort(c(glue("All interactions in `x` must be \\
                     the same width."),
                "i"="Check this with `width(x)`.",
                "i"="Set binSize with `binPairs(x, binSize)`."))
    }

    ## Check & show selection
    .checkBuffer(fg$buffer, bg$buffer)
    buffer <- fg$buffer
    fg <- fg$x
    bg <- bg$x
    .showMultiSelection(fg=fg, bg=bg, buffer=buffer)

    ## Pull Hi-C matrices
    mr <- pixelsToMatrices(x, buffer)

    ## Pull count matrices
    ## Use nBlocks to set blockSize if not provided?
    iarr <- pullHicMatrices(mr, files, binSize, ...)
    cnts <- counts(iarr)

    ## Subset the center pixel of matrix
    cp <- cnts[(buffer+1), (buffer+1),,]

    ## Build array grid
    spacings <- dims <- dim(cnts)
    spacings[3] <- ceiling(spacings[3]/nBlocks)
    grid <- RegularArrayGrid(dims, spacings)

    ## Define functions to apply
    newBody <-
        gsub("fg", "x[fg]", deparse1(body(FUN))) |>
        gsub("bg", "x[bg]", x=_)
    newFUN <- eval(parse(text=paste0("function(x, fg, bg) ", newBody)))
    blockApplyFUN <- \(x) apply(x, c(3,4), \(x) newFUN(x, fg, bg))
    combineFUN <- \(x) do.call("rbind", args=x)

    ## Apply in blocks
    blocks <- blockApply(x=cnts,
                         FUN=blockApplyFUN,
                         grid=grid,
                         verbose=verbose,
                         BPPARAM=BPPARAM)

    ans <- DelayedArray(combineFUN(blocks))
    return(ans)
}

#' Calculate loop enrichment over background.
#'
#' Pulls Hi-C pixels and calculates the enrichment of
#' the selected foreground (`fg`) over the selected
#' background (`bg`).
#'
#' @param x GInteractions object.
#' @param files Character file paths to `.hic` files.
#' @param fg Integer vector of matrix indices for the foreground.
#' @param bg Integer vector of matrix indices for the background.
#' @param FUN Function taking two parameters (i.e., `fg`, `bg`)
#'  defining how enrichment should be calculated. Must produce
#'  a single value (numeric of length one).
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
#'  `pullHicMatrices`. See ?[`pullHicMatrices`].
#' @returns A DelayedMatrix of enrichment scores
#'  where rows are interactions (i.e. loops) and
#'  columns are Hi-C files.
#'
#' @examples
#' ## Read .hic file paths
#' hicFiles <-
#'     system.file("extdata/test_hic", package="mariner") |>
#'     list.files(pattern=".hic", full.names=TRUE)
#'
#' ## Read in loops as GInteractions object
#' loops <-
#'     system.file("extdata", package="mariner") |>
#'     list.files(pattern="WT.*Loops.txt", full.names=TRUE) |>
#'     read.table(header=TRUE) |>
#'     as_ginteractions(keep.extra.columns=FALSE)
#'
#' ## Removes the "chr" prefix for compatibility
#' ## with the preprocessed hic files
#' GenomeInfoDb::seqlevelsStyle(loops) <- 'ENSEMBL'
#'
#' ## Expand binSize of loops
#' loops <- binPairs(x=loops, binSize=100e3)
#'
#' ## Calculate loop enrichment
#' calcLoopEnrichment(x=loops[1:10],
#'                    files=hicFiles)
#'
#' ## Customize different foreground/background
#' ## with selection functions
#' buffer <- 10 # choose pixel radius around center
#' fg <- selectCenterPixel(mhDist=seq(0,4), buffer=buffer)
#' bg <- selectCorners(n=6, buffer=buffer) +
#'     selectOuter(n=2, buffer=buffer)
#' calcLoopEnrichment(x=loops[1:10],
#'                    files=hicFiles,
#'                    fg=fg,
#'                    bg=bg)
#'
#' @rdname calcLoopEnrichment
#' @export
setMethod("calcLoopEnrichment",
          signature(x="GInteractions",
                    files="character"),
          definition=.calcLoopEnrichment)
