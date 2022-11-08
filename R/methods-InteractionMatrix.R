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
.counts <- function(x) {
    ## Check for valid count data
    if (length(assays(x)) == 0) {
        abort("`x` has no counts.")
    }
    if (length(colData(x)) == 0) {
        abort("`x` has no Hi-C files.")
    }
    assay(x, 'counts')
}

#' Access count matrices from
#' InteractionArray or InteractionMatrix
#'
#' @param x InteractionArray or InteractionMatrix object.
#' @returns For InteractionMatrix, a 2-dimensional
#'  DelayedArray is returned with rows representing
#'  interactions in `x` and columns for each Hi-C
#'  file in `files`.
#'
#' @examples
#' #################################
#' ## Accessing Hi-C count matrix ##
#' #################################
#'
#' ## Read .hic file paths
#' hicFiles <-
#'     system.file("extdata/test_hic", package="mariner") |>
#'     list.files(pattern=".hic", full.names=TRUE)
#'
#' ## Create example interactions
#' x <- read.table(text="
#'         9 14000000 14500000 9 14500000 15000000
#'         9 89500000 90000000 9 89500000 90000000
#'         9 23500000 24000000 9 23500000 24000000")
#' x <- as_ginteractions(x)
#'
#' ## Extract 3 pixels from 2 hic files
#' iarr <- pullHicPixels(x, 500e03, hicFiles)
#'
#' ## Access count matrix
#' counts(iarr)
#'
#' @rdname counts
#' @export
setMethod("counts",
          signature(x="InteractionMatrix"),
          definition = .counts)

## Show ------------------------------------------------------------------------

#' Internal show method
#' @noRd
.showInteractionMatrix <- function(object) {
    output <- utils::capture.output(callNextMethod(object))
    if (length(assays(object)) == 0) {
        dims <- rep(0, 2)
    } else {
        dims <- dim(assay(object, 'counts'))
    }
    msg <- "dim: count matrix with %s interactions and %s file(s)"
    output[2] <- do.call(sprintf, c(msg, as.list(dims)))
    cat(output, sep = "\n")
}

#' show for InteractionArray
#' @rdname InteractionArray-class
#' @export
setMethod("show",
          signature(object="InteractionMatrix"),
          definition = .showInteractionMatrix)

## Concatenation ---------------------------------------------------------------

#' rbind InteractionMatrix
#' @include utils.R
#' @rdname InteractionMatrix-class
#' @export
setMethod("rbind", "InteractionMatrix", .rbindIsetDerived)

#' cbind InteractionMatrix
#' @include utils.R
#' @rdname InteractionMatrix-class
#' @export
setMethod("cbind", "InteractionMatrix", .cbindIsetDerived)
