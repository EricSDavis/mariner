## Constructors ----------------------------------------------------------------

#' Internal constructor
#' @inheritParams InteractionMatrix-class
#' @noRd
.newInteractionMatrix <- function(assays, interactions, ...) {
    iset <- InteractionSet(assays, interactions, ...)
    new("InteractionMatrix", iset)
}

#' Constructor
#' @rdname InteractionMatrix-class
#' @export
setMethod("InteractionMatrix", c("ANY", "GInteractions"),
          function(assays, interactions, ...) {
              .newInteractionMatrix(assays, interactions, ...)
          })

#' Constructor
#' @param assays,interactions See
#'   \code{?\link[InteractionSet]{InteractionSet}}
#' @rdname InteractionMatrix-class
#' @export
setMethod("InteractionMatrix", c("missing", "missing"),
          function(assays, interactions, ...) {
              .newInteractionMatrix(list(), GInteractions(), ...)
          })

## Accessors -------------------------------------------------------------------

#' Internal accessor for count matrices
#' @inheritParams counts
#' @importFrom SummarizedExperiment assays
#' @noRd
.counts <- function(object) {
    validObject(object)
    ## Check for valid count data
    if (length(assays(object)) == 0) {
        abort("`object` has no counts.")
    }
    if (length(colData(object)) == 0) {
        abort("`object` has no Hi-C files.")
    }
    assay(object, 'counts')
}

#' Access count matrices from
#' InteractionArray or InteractionMatrix
#'
#' @param object InteractionArray or InteractionMatrix object.
#'
#' @importFrom BiocGenerics counts
#'
#' @returns For InteractionMatrix, a 2-dimensional
#'  DelayedArray is returned with rows representing
#'  interactions in `object` and columns for each Hi-C
#'  file in `files`.
#'
#' @examples
#' #################################
#' ## Accessing Hi-C count matrix ##
#' #################################
#'
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
#' ## Create example interactions
#' x <- read.table(text="
#'         9 14000000 14500000 9 14500000 15000000
#'         9 89500000 90000000 9 89500000 90000000
#'         9 23500000 24000000 9 23500000 24000000")
#' x <- as_ginteractions(x)
#'
#' ## Extract 3 pixels from 2 hic files
#' iarr <- pullHicPixels(x, hicFiles, 500e03)
#'
#' ## Access count matrix
#' counts(iarr)
#'
#' @rdname counts
#' @export
setMethod("counts",
          signature(object="InteractionMatrix"),
          definition = .counts)

#' Replace method for counts
#' @param object InteractionMatrix object
#' @param value Value for replacement
#' @importFrom BiocGenerics "counts<-"
#' @importFrom SummarizedExperiment "assays<-"
#' @returns For InteractionMatrix, the replace matrix
#'  replaces the counts assay with matrix-like
#'  objects supplied in `value`.
#' @examples
#' #################################
#' ## Replacing Hi-C count matrix ##
#' #################################
#'
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
#' ## Create example interactions
#' x <- read.table(text="
#'         9 14000000 14500000 9 14500000 15000000
#'         9 89500000 90000000 9 89500000 90000000
#'         9 23500000 24000000 9 23500000 24000000")
#' x <- as_ginteractions(x)
#'
#' ## Extract 3 pixels from 2 hic files
#' imat <- pullHicPixels(x, hicFiles, 500e03)
#'
#' ## Realize as in-memory matrix
#' counts(imat) <- as.matrix(counts(iarr))
#' counts(imat)
#' imat
#'
#' @rdname counts
#' @export
setMethod("counts<-",
          signature(object="InteractionMatrix"),
          definition = function(object, value) {
              validObject(object)
              assays(object)$counts <- value
              object
          })

#' Internal h5File path accessor for InteractionMatrix
#' @inheritParams path
#' @importFrom rlang abort warn
#' @importFrom glue glue
#' @noRd
.path <- function(object) {
    ## If there are assays...
    if (length(object@assays) == 0) {
        abort(glue("`object` has no counts. \\
                             No HDF5 file path to return."))
    }

    ## And they are HDF5...
    if (!is(object@assays@data@listData$counts, "HDF5Array")) {
        abort(glue("Data not stored on disk. \\
                             No HDF5 file path to return."))
    }

    ## Warn if file doesn't exist
    h5File <- object@assays@data@listData$counts@seed@filepath
    if (!file.exists(h5File)) {
        warn("HDF5 file no longer exists.")
    }

    ## Return path
    return(h5File)
}

#' Accessor for h5File path from an InteractionMatrix
#'
#' Returns the file path describing where the on-disk
#' HDF5 data associated with the InteractionMatrix
#' object is stored.
#'
#' If the file no longer exists, the path is returned
#' along with a warning.
#'
#' @param object InteractionMatrix object
#' @returns The path to the HDF5 file associated with
#'  the InteractionMatrix object.
#' @examples
#' #################################
#' ## Accessing path to HDF5 data ##
#' #################################
#'
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
#' ## Create example interactions
#' x <- read.table(text="
#'         9 14000000 14500000 9 14500000 15000000
#'         9 89500000 90000000 9 89500000 90000000
#'         9 23500000 24000000 9 23500000 24000000")
#' x <- as_ginteractions(x)
#'
#' ## Extract 3 pixels from 2 hic files
#' imat <- pullHicPixels(x, hicFiles, 500e03)
#'
#' ## Access path
#' path(imat)
#'
#' @rdname path
#' @export
setMethod("path",
          signature(object="InteractionMatrix"),
          definition = .path)

#' Internal function to update path to HDF5 file
#' @inheritParams `path<-`
#' @importFrom rlang abort warn
#' @importFrom glue glue
#' @noRd
`.path<-` <- function(object, value) {
    ## If there are assays...
    if (length(object@assays) == 0) {
        abort(glue("`object` has no counts. \\
                   No HDF5 file path to replace."))
    }

    ## And they are HDF5...
    if (!is(object@assays@data@listData$counts, "HDF5Array")) {
        abort(glue("Data not stored on disk. \\
                   No HDF5 file path to replace."))
    }

    ## Update path
    object@assays@data@listData$counts@seed@filepath <- value
    object
}

#' Update path to HDF5 file
#'
#' This method circumvents the `assays<-`
#' and `path<-` methods for updating the HDF5
#' path because they are not accessible when
#' the file path is broken.
#'
#' This allows the file path to be updated even
#' if the original linked data no longer exists.
#'
#' @param object InteractionMatrix object
#' @param value String (length-one character vector)
#'  to use for path replacement.
#'
#' @returns Updates path to HDF5 file for the
#'  InteractionMatrix object.
#'
#' @examples
#' #################################
#' ## Updating path to HDF5 data ##
#' ################################
#'
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
#' ## Create example interactions
#' x <- read.table(text="
#'         9 14000000 14500000 9 14500000 15000000
#'         9 89500000 90000000 9 89500000 90000000
#'         9 23500000 24000000 9 23500000 24000000")
#' x <- as_ginteractions(x)
#'
#' ## Extract 3 pixels from 2 hic files
#' h5File <- tempfile(fileext=".h5")
#' imat <- pullHicPixels(x, hicFiles, 500e03, h5File=h5File)
#'
#' ## Move file to new location
#' newFile <- tempfile(fileext="_new.h5")
#' file.rename(from=h5File, to=newFile)
#'
#' ## Update path
#' path(imat) <- newFile
#' path(imat)
#'
#' @rdname path
#' @export
setMethod("path<-",
          signature(object="InteractionMatrix"),
          definition=`.path<-`)

## Show ------------------------------------------------------------------------

#' Internal show method
#' @noRd
.showInteractionMatrix <- function(object) {

}

#' show for InteractionMatrix
#' @param object InteractionMatrix object.
#' @rdname InteractionMatrix-class
#' @export
setMethod("show", "InteractionMatrix", function(object) {
    validObject(object)
    output <- utils::capture.output(callNextMethod(object))
    if (length(assays(object)) == 0) {
        dims <- rep(0, 2)
    } else {
        dims <- dim(assay(object, 'counts'))
    }
    msg <- "dim: count matrix with %s interactions and %s file(s)"
    output[2] <- do.call(sprintf, c(msg, as.list(dims)))
    cat(output, sep = "\n")
})

## Concatenation ---------------------------------------------------------------

#' rbind InteractionMatrix
#' @param ... InteractionMatrix objects to be combined row-wise.
#'  All objects must be the same class.
#' @param deparse.level An integer scalar; see `?base::rbind` for
#'  a description of this argument.
#' @importFrom BiocGenerics rbind
#' @include utils.R
#' @rdname InteractionMatrix-class
#' @export
setMethod("rbind", "InteractionMatrix", .rbindIsetDerived)

#' cbind InteractionMatrix
#' @param ... InteractionMatrix objects to be combined column-wise.
#'  All objects must be the same class.
#' @param deparse.level An integer scalar; see `?base::cbind` for
#'  a description of this argument.
#' @importFrom BiocGenerics cbind
#' @include utils.R
#' @rdname InteractionMatrix-class
#' @export
setMethod("cbind", "InteractionMatrix", .cbindIsetDerived)
