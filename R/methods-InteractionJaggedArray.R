## Accessors -------------------------------------------------------------------

#' dim for InteractionJaggedArray
#' @param x InteractionJaggedArray object.
#' @importFrom InteractionSet interactions
#' @importFrom rhdf5 h5read
#' @rdname InteractionJaggedArray-class
#' @export
setMethod("dim", "InteractionJaggedArray", function(x) {
    ## Create list to hold dimensions
    dims <- vector("list", 4)
    names(dims) <- c("interactions", "files", "rows", "cols")

    ## Extract dimensions of each matrix
    slices <- h5read(x@counts@h5File, 'slices')

    ## Assign dimensions
    dims[[1]] <- length(interactions(x))
    dims[[2]] <- nrow(colData(x))
    dims[[3]] <- slices[,1]
    dims[[4]] <- slices[,2]
    dims
})

#' interactions
#' @param x InteractionJaggedArray object.
#' @rdname InteractionJaggedArray-class
#' @export
setMethod("interactions", "InteractionJaggedArray", function(x) {
    x@interactions
})

#' metadata
#' @param x InteractionJaggedArray object.
#' @rdname InteractionJaggedArray-class
#' @export
setMethod("metadata", "InteractionJaggedArray", function(x) {
    x@metadata
})

#' colData
#' @param x InteractionJaggedArray object.
#' @rdname InteractionJaggedArray-class
#' @export
setMethod("colData", "InteractionJaggedArray", function(x) {
    x@colData
})

#' counts
#' @param object InteractionJaggedArray object.
#' @rdname InteractionJaggedArray-class
#' @export
setMethod("counts", "InteractionJaggedArray", function(object) {
    object@counts
})

## Show ------------------------------------------------------------------------

#' show for InteractionJaggedArray
#' @param object InteractionJaggedArray object.
#' @importFrom InteractionSet interactions
#' @importFrom SummarizedExperiment colData
#' @rdname InteractionJaggedArray-class
#' @export
setMethod("show", "InteractionJaggedArray", function(object) {
    output <- c()

    ## class & dims
    output[1] <- sprintf("class: %s", class(object))
    dims <- dim(object)
    output[2] <- sprintf(
        "dim: %s interaction(s), %s file(s), variable count matrix(es)",
        dims[[1]], dims[[2]]
    )

    ## metadata
    expt <- names(metadata(object))
    output[3] <- sprintf(
        "metadata(%i): %s",
        length(expt), paste0(expt, collapse=", ")
    )

    ## colData
    output[4] <- sprintf(
        "colData: %s",
        paste0(rownames(colData(object)), collapse=", ")
    )
    cdNames <- names(colData(object))
    output[5] <- sprintf(
        "colData names(%d): %s",
        length(cdNames), paste0(cdNames, collapse=", ")
    )

    ## hdf5
    output[6] <- sprintf("HDF5: %s", object@counts@h5File)

    cat(output, sep = "\n")
})


