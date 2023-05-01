## Show ------------------------------------------------------------------------

#' show for InteractionJaggedArray
#' @param object InteractionJaggedArray object.
#' @importFrom InteractionSet interactions
#' @importFrom SummarizedExperiment colData
#'
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
#' iarr
#'
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
    output[6] <- sprintf("HDF5: %s", path(object))

    cat(output, sep = "\n")
})

## Accessors -------------------------------------------------------------------

#' dim for InteractionJaggedArray
#' @param x InteractionJaggedArray object.
#' @importFrom InteractionSet interactions
#' @importFrom rhdf5 h5read
#' @returns `dim()` returns a list of the dimensions
#'  of the interactions, files, and count matrices.
#' @examples
#' ## Show dimensions
#' dim(iarr)
#'
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
#' @returns `interactions()` returns the interactions.
#' @examples
#' ## Access interactions
#' interactions(iarr)
#'
#' @rdname InteractionJaggedArray-class
#' @export
setMethod("interactions", "InteractionJaggedArray", function(x) {
    x@interactions
})

#' metadata
#' @param x InteractionJaggedArray object.
#' @returns `metadata()` returns the metadata.
#' @examples
#' ## Access metadata
#' metadata(iarr)
#'
#' @rdname InteractionJaggedArray-class
#' @export
setMethod("metadata", "InteractionJaggedArray", function(x) {
    x@metadata
})

#' colData
#' @param x InteractionJaggedArray object.
#' @returns `colData()` returns the column data.
#' @examples
#' ## Access colData
#' colData(iarr)
#'
#' @rdname InteractionJaggedArray-class
#' @export
setMethod("colData", "InteractionJaggedArray", function(x) {
    x@colData
})

#' counts
#' @param object InteractionJaggedArray object.
#' @returns `counts()` returns the JaggedArray object
#'  containing count matrix information.
#' @examples
#' ## Access count matrices
#' counts(iarr)
#'
#' @rdname InteractionJaggedArray-class
#' @export
setMethod("counts", "InteractionJaggedArray", function(object) {
    object@counts
})

#' Access HDF5 path for InteractionJaggedArray
#' @param object JaggedArray object.
#' @returns `path()` returns a character vector with
#'  the path to the HDF5 file with the JaggedArray data.
#' @examples
#' ## Access path to HDF5 data
#' path(iarr)
#'
#' @rdname InteractionJaggedArray-class
#' @export
setMethod("path", "InteractionJaggedArray", function(object) {
    object@counts@h5File
})




