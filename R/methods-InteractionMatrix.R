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
    ## Get the full object
    if (length(assays(x)) == 0) {
        abort("`x` has no counts.")
    }
    assay(x, 'counts')
}

#' Access count matrices from InteractionMatrix
#' @param x InteractionMatrix object.
#' @returns 2-dimensional DelayedMatrix where rows
#'  are interactions and columns are Hi-C files.
#'
#' @examples
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
#' ## Extract 3, 11x11 count matrices from 2 hic files
#' iarr <- pullHicPixels(x, 500e03, hicFiles)
#'
#' ## Access count matrices
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
