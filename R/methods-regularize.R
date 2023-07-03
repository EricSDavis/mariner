## .rescale and .resizeMat functions adapted from
## https://stackoverflow.com/questions/11123152/\
## function-for-resizing-matrices-in-r

#' Internal function for rescaling
#' @param x Numeric vector of positions
#' @param newrange Numeric vector of newrange
#' @noRd
.rescale <- function(x, newrange){
    xrange <- range(x)
    mfac <- (newrange[2]-newrange[1])/(xrange[2]-xrange[1])
    if (is.infinite(mfac)) mfac <- 1
    newrange[1]+(x-xrange[1])*mfac
}

#' Internal function for resizing matrices
#' @param x A matrix.
#' @param ndim Numeric vector of the number of rows
#'  and columns of the new matrix.
#' @returns Resized matrix
#' @noRd
.resizeMat <- function(x, ndim){
    if(!requireNamespace("fields")) stop("`fields` package required.")
    
    odim <- dim(x)

    ## Modify for single row/column matrices
    z <- x
    if (odim[1] == 1L) {
        z <- rbind(z,z)
        odim[1] <- 2
    }
    if (odim[2] == 1L) {
        z <- cbind(z,z)
        odim[2] <- 2
    }

    ## initialize input & output objects
    obj <- list(x=seq_len(odim[1]), y=seq_len(odim[2]), z=z)
    ans <- matrix(NA, nrow=ndim[1], ncol=ndim[2])

    ## rescaling
    loc <- ncord <- as.matrix(expand.grid(seq_len(ndim[1]),
                                          seq_len(ndim[2])))
    loc[,1] <- .rescale(ncord[,1], c(1,odim[1]))
    loc[,2] <- .rescale(ncord[,2], c(1,odim[2]))
    
    ## interpolation
    ans[ncord] <- fields::interp.surface(obj, loc)
    ans
}


#' Regularize JaggedArray objects
#' and return as DelayedArray objects
#'
#' @inheritParams regularize
#' @importFrom rlang abort
#' @importFrom rhdf5 h5createFile h5createDataset h5write
#' @importFrom HDF5Array HDF5Array
#' @importFrom DelayedArray DelayedArray
#' @noRd
.regularizeJA <- function(x, ndim, h5File, scale,
                          nBlocks, verbose,
                          chunkSize, compressionLevel){

    ## Check params
    if (nBlocks <= 0) abort("`nBlocks` must be > 0.")

    ## Create HDF5 file for storing data
    h5createFile(h5File)
    h5createDataset(
        file=h5File,
        dataset="counts",
        dims=c(ndim, x@dim),
        chunk=c(ndim, chunkSize, x@dim[2]),
        storage.mode="double",
        fillValue=NA_real_,
        level=compressionLevel
    )

    ## Split indices into blocks
    ni <- x@dim[1]
    ind <- seq_len(ni)
    blocks <- split(ind, ceiling(ind/(ni/nBlocks)))

    ## Loop through each block and Hi-C file
    for (j in seq_len(x@dim[2])) {
        bid <- 1
        for (block in blocks) {

            ## Subset array into list of matrices
            if (verbose) {
                message("/ Reading and realizing block ",
                        bid, "/", length(blocks),
                        " of file ",
                        j, "/", x@dim[2], " ... ",
                        appendLF=FALSE)
            }
            subArray <- x[,,block,j]
            matList <- list()
            if (is(subArray, "DelayedArray")) {
                if (length(block) == 1L) {
                    sad <- dim(subArray)
                    matList <- matrix(
                        data=subArray,
                        nrow=sad[1],
                        ncol=sad[2]
                    ) |> list()
                } else {
                    matList <- lapply(asplit(subArray, 3), drop)
                }
            }
            if (is(subArray, "JaggedArray")) {
                matList <- as.list(subArray)[[1]]
            }
            if (verbose) message("OK")

            ## Resize & normalize matrices
            if (verbose) {
                message("\\ Processing it ... ", appendLF=FALSE)
            }
            ans <- lapply(matList, \(x) {
                m <- .resizeMat(x, ndim=ndim)
                ifelse(scale, m <- m/sum(m), m <- m)
                m[is.na(m)] <- 0
                m
            })

            ## Turn list of matrices into array for storage
            a <- array(unlist(ans), dim=c(ndim, length(matList), 1))

            ## Store array as an HDF5-backed delayed array
            h5write(
                obj=a,
                file=h5File,
                name="counts",
                index=list(
                    seq_len(ndim[1]),
                    seq_len(ndim[2]),
                    block, j
                )
            )
            if (verbose) message("OK")
            bid <- bid + 1
        }
    }

    ## Return as HDF5-backed array
    DelayedArray(HDF5Array(h5File, "counts"))
}

#' Regularize InteractionJaggedArray objects
#' and return as InteractionArray
#'
#' @inheritParams regularize
#' @importFrom rlang abort
#' @importFrom rhdf5 h5createFile h5createDataset h5write
#' @importFrom HDF5Array HDF5Array
#' @importFrom DelayedArray DelayedArray
#' @noRd
.regularizeIJA <- function(x, ndim, h5File, scale, nBlocks,
                           verbose, chunkSize, compressionLevel) {

    ## Regularize the JaggedArray
    ja <- counts(x)
    da <- regularize(ja, ndim, h5File, scale, nBlocks,
                     verbose, chunkSize, compressionLevel)
    da <- aperm(da, c(3,4,1,2))

    ## Construct and return InteractionArray
    InteractionArray(
        assays = list(
            counts=da
        ),
        interactions=interactions(x),
        colData=colData(x),
        metadata=metadata(x)
    )

}

#' Regularize JaggedArray or
#' InteractionJaggedArray objects
#'
#' InteractionJaggedArray objects and their
#' count matrices (JaggedArray objects) contain
#' variable dimension matrices. The `regularize`
#' function resizes these matrices to the new
#' dimensions supplied in `ndim`. The result is
#' a DelayedArray object (for JaggedArray) or
#' an InteractionArray object (for
#' InteractionJaggedArray).
#'
#' Note that the interaction/binSize/count matrices
#' relationship will be inconsistent in the resulting
#' InteractionArray object and the row/col names will
#' not be available.
#'
#' @param x A JaggedArray or InteractionJaggedArray object.
#' @param ndim Numeric vector of length two describing
#'  the new dimensions of the output matrices.
#' @param h5File Character file path to save `.h5` file.
#' @param scale Boolean (TRUE/FALSE) indicating whether
#'  the values in the new matrices should be scaled to the
#'  total signal in each matrix.
#' @param nBlocks Number of blocks for block-processing JaggedArrays.
#'  Default is 5. Increase this for large datasets. To read and process
#'  all data at once, set this value to 1.
#' @param verbose Boolean (TRUE or FALSE) describing
#'  whether to report block-processing progress.
#' @param chunkSize Number (length one numeric vector)
#'  indicating how many values of `x` to chunk for each
#'  write to HDF5 stored data. This has downstream
#'  implications for accessing subsets later. For small
#'  `compressionLevel` values use smaller `chunkSize`
#'  values and for large `compressionLevel` values use large
#'  (i.e. `length(x)`) values to improve performance.
#' @param compressionLevel Number (length one numeric vector)
#'  between 0 (Default) and 9 indicating the compression
#'  level used on HDF5 file.
#' @param ... Additional arguments.
#'
#' @returns If `x` is a JaggedArray then `regularize` returns
#'  an HDF5-backed 4-dimensional DelayedArray object
#'  where the first and second dimensions are the rows and columns
#'  of the count matrices (`ndim`), the third dimension is the
#'  number of interactions and the fourth dimension is the
#'  number of files. If `x` is an InteractionJaggedArray then
#'  an InteractionArray object is returned where counts returns
#'  the object described above.
#'
#' @examples
#' ## Load marinerData
#' if (!require("marinerData", quietly = TRUE))
#'     BiocManager::install("marinerData")
#'
#' ## Read .hic file paths
#' hicFiles <- c(
#'     LEUK_HEK_PJA27_inter_30.hic(),
#'     LEUK_HEK_PJA30_inter_30.hic()
#' )
#' names(hicFiles) <- c("FS", "WT")
#'
#' ## Create test interactions
#' gi <- read.table(text="
#'             1 51000000 51300000 1 51000000 51500000
#'             2 52000000 52300000 3 52000000 52500000
#'             1 150000000 150500000 1 150000000 150300000
#'             2 52000000 52300000 2 52000000 52800000") |>
#'     as_ginteractions()
#' gi <- c(gi,gi) # make more interactions
#'
#' ## InteractionJaggedArray object
#' ija <- pullHicMatrices(gi, hicFiles, 100e03, half="both")
#'
#' ## Regularize InteractionJaggedArray
#' ia <- regularize(ija, ndim=c(5,5), nBlocks=1)
#' aggHicMatrices(ia, nBlocks=1)
#'
#' ## Regularize JaggedArray
#' ja <- counts(ija)
#' regularize(ja, ndim=c(5,5), nBlocks=1)
#'
#' @rdname regularize
#' @export
setMethod("regularize", "JaggedArray", .regularizeJA)

#' @rdname regularize
#' @export
setMethod("regularize", "InteractionJaggedArray", .regularizeIJA)

