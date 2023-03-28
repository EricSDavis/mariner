## Constructors ----------------------------------------------------------------

#' Internal constructor
#' @importFrom SummarizedExperiment assays
#' @importFrom InteractionSet InteractionSet
#' @noRd
.newInteractionArray <- function(assays, interactions, ...) {
    iset <- InteractionSet(assays, interactions, ...)
    new("InteractionArray", iset)
}

#' Constructor
#' @rdname InteractionArray-class
#' @export
setMethod("InteractionArray", c("ANY", "GInteractions"),
          function(assays, interactions, ...) {
              .newInteractionArray(assays, interactions, ...)
          })

#' Constructor
#' @param assays,interactions See
#'   \code{?\link[InteractionSet]{InteractionSet}}
#' @param ... Additional arguments.
#' @rdname InteractionArray-class
#' @export
setMethod("InteractionArray", c("missing", "missing"),
          function(assays, interactions, ...) {
              .newInteractionArray(list(), GInteractions(), ...)
          })

## Accessors -------------------------------------------------------------------

#' Internal accessor for count matrices
#' from InteractionArray object
#' @inheritParams counts
#' @importFrom SummarizedExperiment assays
#' @importFrom DelayedArray DelayedArray
#' @importFrom rlang abort
#' @importFrom utils capture.output
#' @noRd
.counts <- function(object, showDimnames=FALSE) {

    ## Check for valid count data
    if (length(assays(object)) == 0) {
        abort("`object` has no count matrices.")
    }
    if (length(colData(object)) == 0) {
        abort("`object` has no Hi-C files.")
    }

    ## Get the full object
    cnts <- aperm(assay(object, 'counts'), c(3,4,1,2))

    if (showDimnames) {
        return(CountMatrix(seed=cnts, object=object))
    } else {
        return(cnts)
    }

}


#' Access count matrices from
#' InteractionArray or InteractionMatrix
#'
#' @param object InteractionArray or InteractionMatrix object.
#' @param showDimnames Logical vector of length-one
#'  indicating whether to show dimensions of
#'  count matrices (default FALSE). Only applicable for
#'  InteractionArray objects.
#'
#' @importFrom BiocGenerics counts
#'
#' @returns For InteractionArray, a 4-dimensional
#'  DelayedArray of Hi-C submatrices is returned with
#'  the following dimensions: rows of count matrix,
#'  columns of count matrix, Interactions in `object`,
#'  Hi-C `files`.
#'
#' @examples
#' ######################################
#' ## Accessing Hi-C count submatrices ##
#' ######################################
#'
#' ## Read .hic file paths
#' hicFiles <-
#'     system.file("extdata/test_hic", package="mariner") |>
#'     list.files(pattern=".hic", full.names=TRUE)
#'
#' ## Create example interactions
#' x <- read.table(text="
#'         9 14435000 14490000 9 14740000 14795000
#'         9 89540000 89595000 9 89785000 89840000
#'         9 23700000 23755000 9 23760000 23815000")
#' x <- as_ginteractions(x)
#'
#' ## Extract 3, 11x11 count matrices from 2 hic files
#' iarr <- pullHicMatrices(x, hicFiles, 5e03)
#'
#' ## Access count matrices
#' counts(iarr)
#' counts(iarr, FALSE)
#'
#' @rdname counts
#' @export
setMethod("counts",
          signature(object="InteractionArray"),
          definition=.counts)

## Show ------------------------------------------------------------------------

#' show for InteractionArray
#' @param object InteractionArray object.
#' @rdname InteractionArray-class
#' @export
setMethod("show", "InteractionArray", function(object) {
    output <- utils::capture.output(callNextMethod(object))
    if (length(assays(object)) == 0) {
        dims <- rep(0, 4)
    } else {
        dims <- dim(assay(object, 'counts'))
    }
    msg <- "dim: %s interaction(s), %s file(s), %sx%s count matrix(es)"
    output[2] <- do.call(sprintf, c(msg, as.list(dims)))
    cat(output, sep = "\n")
})


## Concatenation ---------------------------------------------------------------

#' rbind InteractionArray
#' @param ... InteractionArray objects to be combined row-wise.
#'  All objects must be the same class.
#' @param deparse.level An integer scalar; see `?base::rbind` for
#'  a description of this argument.
#' @importFrom BiocGenerics rbind
#' @include utils.R
#' @rdname InteractionArray-class
#' @export
setMethod("rbind", "InteractionArray", .rbindIsetDerived)

#' cbind InteractionArray
#' @param ... InteractionArray objects to be combined column-wise.
#'  All objects must be the same class.
#' @param deparse.level An integer scalar; see `?base::cbind` for
#'  a description of this argument.
#' @importFrom BiocGenerics cbind
#' @include utils.R
#' @rdname InteractionArray-class
#' @export
setMethod("cbind", "InteractionArray", .cbindIsetDerived)

