## Show method -----------------------------------------------------------------

#' show for JaggedArray
#' @param object JaggedArray object.
#' @importFrom DelayedArray DelayedArray
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
#' ## Create test interactions
#' gi <- read.table(text="
#'             1 51000000 51300000 1 51000000 51500000
#'             2 52000000 52300000 3 52000000 52500000
#'             1 150000000 150500000 1 150000000 150300000
#'             2 52000000 52300000 2 52000000 52800000") |>
#'     as_ginteractions()
#'
#' ## InteractionJaggedArray object
#' iarr <- pullHicMatrices(gi, hicFiles, 100e03, half="both")
#' arr <- counts(iarr)
#' arr
#'
#' @rdname JaggedArray-class
#' @export
setMethod("show", "JaggedArray", function(object) {
    dims <- object@dim

    ## Show first in set
    first <- object[1,1]
    fd <- dim(first)
    ans <- sprintf(
        "<%i x %i x %i x %i> %s:",
        dims[1], dims[2], fd[1], fd[2], class(object)
    )
    ans <- c(ans, "1, 1,")
    ans <- c(ans, sprintf("<%i x %i> matrix", fd[1], fd[2]))
    ans <- c(ans, capture.output(first)[-c(1,2)])

    ## Show last if not the same as first
    if (any(dims != 1L)) {
        ans[1] <- sprintf(
            "<%i x %i x n x m> %s:",
            dims[1], dims[2], class(object)
        )
        last <- object[dims[1], dims[2]]
        ld <- dim(last)
        ans <- c(ans, "...\n")
        ans <- c(ans, sprintf("%i, %i,", dims[1], dims[2]))
        ans <- c(ans, sprintf("<%i x %i> matrix", ld[1], ld[2]))
        ans <- c(ans, capture.output(last)[-c(1,2)])
    }
    cat(ans, sep='\n')
})


## Subsetting ------------------------------------------------------------------

#' Internal indexing for JaggedArray
#' without modifying underlying HDF5 data
#' @inheritParams [
#' @importFrom rhdf5 h5read
#' @noRd
.accessJaggedArray <- function(x, i, j, ...) {
    # i=interactions, j=files
    dims <- dim(x)
    if (missing(i) | is.null(i)) i <- seq_len(dims[[1]])
    if (missing(j) | is.null(j)) j <- seq_len(dims[[2]])

    ## Apply delayed subset operation
    if (!is.null(x@subsetTree[[1]])) i <- x@subsetTree[[1]][i]
    if (!is.null(x@subsetTree[[2]])) j <- x@subsetTree[[2]][j]

    slices <- h5read(x@h5File, 'slices', index=list(i, seq_len(4)))
    ans <-
        lapply(j, \(m) {
            lapply(seq_len(nrow(slices)), \(s) {
                slice <- seq(slices[s,3], slices[s,4])
                cnts <- h5read(x@h5File, 'counts', index=list(slice, m))
                matrix(
                    data=cnts,
                    nrow=slices[s,1],
                    ncol=slices[s,2],
                    byrow=TRUE
                )
            })
        })

    ans
}

#' Coerce JaggedArray to DelayedArray
#' @param x JaggedArray
#' @returns A DelayedArray object if the conversion
#'  is possible (i.e. all dimensions are the same),
#'  otherwise return the original JaggedArray object.
#' @noRd
.JaggedArrayToDelayedArray <- function(x) {
    ## Read matrix dimensions from HDF5
    slices <- h5read(
        file=x@h5File,
        name='slices',
        index=list(x@subsetTree[[1]], seq_len(4))
    )

    ## Return jagged array if dimensions
    ## are not all the same
    if (!(length(unique(slices[,1])) == 1L &
          length(unique(slices[,2])) == 1L)) {
        return(x)
    }

    ## Otherwise build a DelayedArray.
    ## Since arrays fill from first dimension in R
    ## we need fill columns first then permute.
    cnts <- .accessJaggedArray(x, NULL, NULL) |> unlist()
    a <- array(data=as.vector(cnts),
               dim=c(unique(slices[,1]),
                     unique(slices[,2]),
                     x@dim))
    a <- DelayedArray(a)
    return(a)
}

#' Internal subsetting for JaggedArray
#' that records the subset operation
#' @inheritParams [
#' @importFrom rhdf5 h5read h5write
#' @noRd
.subsetJaggedArray <- function(x, i=NULL, j=NULL) {
    # i=interactions, j=files
    ## Handle missing dims
    dims <- x@dim
    if (missing(i) | is.null(i)) i <- seq_len(dims[[1]])
    if (missing(j) | is.null(j)) j <- seq_len(dims[[2]])

    ## Record subset operation
    if (is.null(x@subsetTree[[1]])) {
        x@subsetTree[[1]] <- i
    } else {
        x@subsetTree[[1]] <- x@subsetTree[[1]][i]
    }
    if (is.null(x@subsetTree[[2]])) {
        x@subsetTree[[2]] <- j
    } else {
        x@subsetTree[[2]] <- x@subsetTree[[2]][j]
    }

    ## Update outer dims
    x@dim[1] <- length(i)
    x@dim[2] <- length(j)

    ## Coerce to DelayedArray if possible
    x <- .JaggedArrayToDelayedArray(x)
    return(x)
}

#' Indexing for JaggedArray
#'
#' Subset a JaggedArray by its interactions
#' ([i,]) or its Hi-C files ([,j]).
#'
#' The object returned will be a JaggedArray
#' if the submatrices contain different dimensions.
#' However, the returned object will automatically
#' be coerced into a DelayedArray if possible (i.e.
#' the dimensions of the rows and columns are the same.)
#'
#' The JaggedArray data is still stored on-disk in
#' an HDF5 file until it is coerced into a DelayedArray
#' or realized as a list of matrices.
#'
#' @param x A JaggedArray object.
#' @param i Numeric vector indicating the indices
#'  of interactions to extract.
#' @param j Numeric vector indicating the indices
#'  of files to extract.
#' @returns Subsetting returns a JaggedArray or DelayedArray object
#'  (see Details).
#' @examples
#' ## Subsetting
#' arr[1,] # DelayedArray
#' arr[,1] # JaggedArray
#'
#' @rdname JaggedArray-class
#' @export
setMethod("[", "JaggedArray", function(x, i, j) {
    if (missing(i)) i <- NULL
    if (missing(j)) j <- NULL
    .subsetJaggedArray(x, i, j)
})

## Accessors -------------------------------------------------------------------

#' Realize as list for JaggedArray
#'
#' `as.list` reads the on-disk data and returns it
#' as an in-memory list of matrices.
#'
#' @param x JaggedArray object.
#' @returns `as.list()` returns a list of matrices.
#' @examples
#' ## Realize as list
#' as.list(arr)
#'
#' @rdname JaggedArray-class
#' @export
setMethod("as.list", "JaggedArray", function(x) {
    .accessJaggedArray(x, NULL, NULL)
})

#' Access HDF5 path for JaggedArray
#' @param object JaggedArray object.
#' @returns `path()` returns a character vector with the path to
#'  the HDF5 file with the JaggedArray data.
#'
#' @examples
#' ## Find the data path
#' path(arr)
#'
#' @rdname JaggedArray-class
#' @export
setMethod("path", "JaggedArray", function(object) {
    object@h5File
})


