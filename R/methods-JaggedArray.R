#' Internal show function for JaggedArray
#' @inheritParams show,JaggedArray
#' @importFrom DelayedArray DelayedArray
#' @noRd
.showJaggedArray <- function(object) {
    dims <- object@dim

    ## Show first in set
    first <- .accessJaggedArray(object, 1, 1)
    fd <- dim(first)
    ans <- sprintf(
        "<%i x %i x %i x %i> %s:",
        dims[1], dims[2], fd[1], fd[2], class(object)
    )
    ans <- c(ans, "1, 1,")
    ans <- c(ans, sprintf("<%i x %i> matrix", fd[1], fd[2]))
    ans <- c(ans, capture.output(DelayedArray(first))[-1])

    ## Show last if not the same as first
    if (any(dims != 1L)) {
        ans[1] <- sprintf(
            "<%i x %i x n x m> %s:",
            dims[1], dims[2], class(object)
        )
        last <- .accessJaggedArray(object, dims[1], dims[2])
        ld <- dim(last)
        ans <- c(ans, "\n...\n")
        ans <- c(ans, sprintf("%i, %i,", dims[1], dims[2]))
        ans <- c(ans, sprintf("<%i x %i> matrix", ld[1], ld[2]))
        ans <- c(ans, capture.output(DelayedArray(last))[-1])
    }
    cat(ans, sep='\n')
}

#' show for JaggedArray
#' @param object JaggedArray object.
#' @rdname JaggedArray-class
#' @export
setMethod("show", "JaggedArray", .showJaggedArray)


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

    ## Assign delayed subsetting operation
    if (!is.null(x@subsetTree[[1]])) i <- x@subsetTree[[1]][i]
    if (!is.null(x@subsetTree[[2]])) j <- x@subsetTree[[2]][j]

    slices <- h5read(x@h5File, 'slices', index=list(i, seq_len(4)))
    ans <-
        lapply(j, \(m) {
            lapply(seq_len(nrow(slices)), \(s) {
                slice <- seq(slices[s,3], slices[s,4])
                cnts <- h5read(x@h5File, 'counts', index=list(slice, m))
                matrix(data=cnts, nrow=slices[s,1], ncol=slices[s,2])
            })
        })

    if (length(i) == 1 & length(j) == 1) {
        ans <- unlist(ans, recursive=FALSE)[[1]]
    }
    ans
}

#' Internal subsetting for JaggedArray
#' that records the subsetting operation
#' @inheritParams [
#' @importFrom rhdf5 h5read h5write
#' @noRd
.subsetJaggedArray <- function(x, i=NULL, j=NULL) {
    # i=interactions, j=files
    ## Handle missing dims
    dims <- x@dim
    if (missing(i) | is.null(i)) i <- seq_len(dims[[1]])
    if (missing(j) | is.null(j)) j <- seq_len(dims[[2]])

    ## Assign to delayed subsetting operation
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

    return(x)
}


#' Indexing for JaggedArray
#' @param object JaggedArray object.
#' @rdname JaggedArray-class
#' @export
setMethod("[", "JaggedArray", function(x, i, j) {
    if (missing(i)) i <- NULL
    if (missing(j)) j <- NULL
    .subsetJaggedArray(x, i, j)
})

#' Realize as list for JaggedArray
#' @param object JaggedArray object.
#' @rdname JaggedArray-class
#' @export
setMethod("as.list", "JaggedArray", function(x) {
    .accessJaggedArray(x, NULL, NULL)
})
