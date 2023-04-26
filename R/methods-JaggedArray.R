#' Internal show function for JaggedArray
#' @inheritParams show,JaggedArray
#' @noRd
.showJaggedArray <- function(object) {
    dims <- object@dim
    first <- .accessJaggedArray(object, 1, 1)
    last <- .accessJaggedArray(object, dims[1], dims[2])
    ans <- ", , 1, 1"
    ans <- c(ans, capture.output(first))
    ans <- c(ans, "\n...\n")
    ans <- c(ans, sprintf(", , %i, %i", dims[1], dims[2]))
    ans <- c(ans, capture.output(last))
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
.accessJaggedArray <- function(x, i, j) {
    # i=interactions, j=files
    dims <- dim(x)
    if (missing(i)) i <- seq_len(dims[[1]])
    if (missing(j)) j <- seq_len(dims[[2]])
    slices <- h5read(x@h5File, 'slices', index=list(i, seq_len(4)))
    ans <-
        lapply(j, \(k) {
            lapply(seq_len(nrow(slices)), \(s) {
                slice <- seq(slices[s,3], slices[s,4])
                cnts <- h5read(x@h5File, 'counts', index=list(slice, k))
                matrix(data=cnts, nrow=slices[s,1], ncol=slices[s,2])
            })
        })

    if (length(i) == 1 & length(j) == 1) {
        ans <- unlist(ans, recursive=FALSE)[[1]]
    }
    ans
}

#' Internal subsetting for JaggedArray
#' that modifies the underlying HDF5 data
#' @inheritParams [
#' @importFrom rhdf5 h5read h5write
#' @noRd
.subsetJaggedArray <- function(x, i, j, h5File=NULL) {
    # i=interactions, j=files
}


#' Indexing for JaggedArray
#' @param object JaggedArray object.
#' @rdname JaggedArray-class
#' @export
setMethod("[", "JaggedArray", function(x, i, j) {
    .accessJaggedArray(x, i, j)
})
