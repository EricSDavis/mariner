## Constructors ----------------------------------------------------------------

#' Internal constructor
#' @inheritParams InteractionMatrix
#' @noRd
.newInteractionMatrix <- function(assays, interactions, ...) {
    iset <- InteractionSet(assays, interactions, ...)
    new("InteractionMatrix", iset)
}

#' InteractionMatrix-Constructor
#' @param assays,interaction See
#'   \code{? \link[InteractionSet]{InteractionSet}}
setMethod("InteractionMatrix", c("ANY", "GInteractions"),
          function(assays, interactions, ...) {
              .newInteractionMatrix(assays, interactions, ...)
          })

#' Method
setMethod("InteractionMatrix", c("missing", "missing"),
          function(assays, interactions, ...) {
              .newInteractionMatrix(list(), GInteractions(), ...)
          })

## Accessors -------------------------------------------------------------------

#' Accessors
#' @param x InteractionMatrix
#' @noRd
.counts <- function(x) {
    assay(x, 'counts')
}

#'
#'
#' @rdname InteractionMatrix
#' @export
setMethod("counts",
          signature(x="InteractionMatrix"),
          definition = .counts)

#' Show
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

#' Show
setMethod("show",
          signature(object="InteractionMatrix"),
          definition = .showInteractionMatrix)
