#' Internal method for aggregating all interactions
#' and files
#' @importFrom DelayedArray RegularArrayGrid blockApply
#' @importFrom BiocParallel bpparam
#' @importFrom abind abind
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
#' @importFrom InteractionSet GInteractions
#' @returns An aggregated DelayedMatrix object.
#' @noRd
.aggAll <- function(x, FUN, nBlocks, verbose, BPPARAM) {

    ## Define aggregation function to apply
    applyFUN <- \(x) apply(x,c(1,2),FUN)

    ## Build array grid
    a <- counts(x)
    spacings <- dims <- dim(a)
    spacings[3] <- ceiling(spacings[3]/nBlocks)
    grid <- RegularArrayGrid(dims, spacings)

    ## Apply in blocks
    blocks <- blockApply(x=a,
                         FUN=applyFUN,
                         grid=grid,
                         verbose=verbose,
                         BPPARAM=BPPARAM)

    ## Bind into array and aggregate blocks
    ans <- applyFUN(abind(blocks, along=3))

    return(DelayedArray(ans))
}

#' Internal method for aggregating all interactions
#' @importFrom DelayedArray RegularArrayGrid blockApply
#' @importFrom BiocParallel bpparam
#' @importFrom abind abind
#' @returns A 3-dimensional `DelayedArray` where
#'  rows/cols represent Hi-C submatrices
#'  and the 3rd dimension represents files.
#' @noRd
.aggByFiles <- function(x, FUN, nBlocks, verbose, BPPARAM) {

    ## Define aggregation functions to apply
    blockApplyFUN <- \(x) apply(x,c(1,2,4),FUN)
    combineFUN <- \(x) apply(abind(x, along=4), c(1,2,3), FUN)

    ## Build array grid
    a <- counts(x)
    spacings <- dims <- dim(a)
    spacings[3] <- ceiling(spacings[3]/nBlocks)
    grid <- RegularArrayGrid(dims, spacings)

    ## Apply in blocks
    blocks <- blockApply(x=a,
                         FUN=blockApplyFUN,
                         grid=grid,
                         verbose=verbose,
                         BPPARAM=BPPARAM)

    ## Bind into array and aggregate blocks
    ans <- DelayedArray(combineFUN(blocks))

    return(ans)
}

#' Aggregate InteractionArray by interactions
#'
#' Aggregates by interactions in chunks. Since
#' these data are read and written to HDF5 files,
#' parallelization is currently not supported.
#'
#' @inheritParams aggregateCounts
#' @importFrom SummarizedExperiment assay
#' @importFrom DelayedArray RegularArrayGrid DelayedArray
#' @importFrom rhdf5 h5createFile
#' @importFrom HDF5Array HDF5RealizationSink
#' @importFrom abind abind
#' @importFrom InteractionSet interactions
#' @returns A 3-dimensional `DelayedArray` where
#'  rows/cols represent Hi-C submatrices
#'  and the 3rd dimension represents
#'  interactions.
#' @seealso [HDF5Array::write_block]
#' @noRd
.aggByInteractions <- function(x, FUN, nBlocks, verbose,
                               compressionLevel) {

    ## Init & define dims
    a <- counts(x)
    spacings <- dims <- dim(a)
    spacings[3] <- ceiling(spacings[3]/nBlocks)
    chunkDims <- storageDims <- dims[c(1,2,3)]
    chunkDims[3] <- spacings[3]

    ## Create grid for sink
    grid <- RegularArrayGrid(dims, spacings)
    sink_grid <- RegularArrayGrid(storageDims, chunkDims)

    ## Create HDF5-backed realization sink
    h5 <- tempfile(fileext = ".h5")
    h5createFile(h5)
    sink <- HDF5RealizationSink(filepath=h5,
                                name="counts",
                                type="double",
                                dim=storageDims,
                                chunkdim=chunkDims,
                                level=compressionLevel)

    ## Wrap function that operates on each block
    applyFUN <- \(block) apply(block, c(1,2,3), FUN)

    ## Read, apply, and write to HDF5
    ans <- hdf5BlockApply(x=a,
                          FUN=applyFUN,
                          sink=sink,
                          grid=grid,
                          sink_grid=sink_grid,
                          verbose=verbose)

    return(as(ans, "DelayedArray"))
}

#' Internal method for aggregating
#' count matrices from InteractionArray
#' objects
#' @importFrom rlang arg_match inform
#' @noRd
.aggInteractionArray <- function(x, by, FUN, nBlocks, verbose,
                                 BPPARAM, compressionLevel) {

    ## Parse shared arguments
    .checkTypes(c(nBlocks="number", verbose="boolean"))
    if (nBlocks <= 0) abort("`nBlocks` must be > 0.")
    FUN <- match.fun(FUN)
    if (nBlocks > length(x)) {
        nBlocks <- length(x)
        inform(c("nBlocks > length(x)",
                 "i"="setting nBlocks to length(x)"))
    }

    ## Dispatch aggregation functions
    ## Default is aggregate all files & interactions
    if (is.null(by)) {
        ans <- .aggAll(x, FUN, nBlocks, verbose, BPPARAM)
    } else {

        ## If "by" is provided check it and dispatch accordingly
        .checkTypes(c(by="string"))
        by <- arg_match(by, c("files", "interactions"))

        if (identical(by, 'files')) {
            ans <- .aggByFiles(x, FUN, nBlocks, verbose, BPPARAM)
        }

        if (identical(by, 'interactions')) {
            .checkTypes(c(compressionLevel="number"))
            ans <- .aggByInteractions(x, FUN, nBlocks, verbose,
                                      compressionLevel)
        }
    }
    return(ans)
}

#' Aggregate count matrices from
#' InteractionArray objects
#'
#' Aggregation of count matrices is done
#' blocks to avoid large memory usage. Use
#' `nBlocks` to control the number of blocks
#' read into memory at once. Blocks are defined
#' as `length(interactions(x))/nBlocks`.
#'
#' Since interactions are typically the largest
#' dimension in an InteractionArray, using
#' `by=interactions` creates an HDF5-backed array
#' to store these large arrays. Currently parallel
#' processing for HDF5-backed arrays are not
#' supported regardless of the value of `BPPARAM`.
#'
#' Both `by=NULL` and `by=files` support parallel
#' processing.
#'
#' @param x InteractionArray object.
#' @param by String (length one character vector)
#'  describing whether to aggregate by interactions,
#'  files, or neither (i.e. NULL as default).
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
#' @param compressionLevel Number (length one numeric vector)
#'  between 0 (Default) and 9 indicating the compression
#'  level used on HDF5 file.
#' @returns An aggregated `DelayedArray` object.
#'  If `by=interactions` or `by=files` then a 3-dimensional
#'  `DelayedArray` is returned. If `by=NULL` (default) then
#'  A 2-dimensional `DelayedMatrix` is returned.
#' @examples
#' ## Load marinerData
#' if (!require("marinerData", quietly = TRUE))
#'     install.packages("marinerData")
#'
#' ## Read .hic file paths
#' hicFiles <- c(
#'     marinerData::LEUK_HEK_PJA27_inter_30.hic(),
#'     marinerData::LEUK_HEK_PJA30_inter_30.hic()
#' )
#' names(hicFiles) <- c("FS", "WT")
#'
#' ## Read in loops as GInteractions object
#' loops <-
#'     WT_5kbLoops.txt() |>
#'     setNames("WT") |>
#'     read.table(header=TRUE) |>
#'     as_ginteractions(keep.extra.columns=FALSE)
#'
#' ## Removes the "chr" prefix for compatibility
#' ## with the preprocessed hic files
#' GenomeInfoDb::seqlevelsStyle(loops) <- 'ENSEMBL'
#'
#' ## Expand pixel ranges with a 5 pixel buffer on either side
#' loops <-
#'     binPairs(loops, binSize=100e3) |>
#'     pixelsToMatrices(buffer=5)
#'
#' ## Extract 10, 11x11 count matrices from 2 hic files
#' iarr <-
#'     loops[1:10] |>
#'     pullHicMatrices(binSize=100e3,
#'                     files=hicFiles)
#'
#' ## Aggregate all, by files, or by interactions
#' aggHicMatrices(x=iarr)
#' aggHicMatrices(x=iarr, by="files")
#' aggHicMatrices(x=iarr, by="interactions")
#'
#' @rdname aggHicMatrices
#' @export
setMethod("aggHicMatrices",
          signature(x="InteractionArray"),
          definition = .aggInteractionArray)
