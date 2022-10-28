## Constructors ----------------------------------------------------------------

#' Internal constructor
#' @importFrom SummarizedExperiment assays
#' @importFrom InteractionSet InteractionSet
#' @noRd
.newInteractionArray <- function(assays, interactions, ...) {
    iset <- InteractionSet(assays, interactions, ...)
    new("InteractionArray", iset)
}

#' Method
setMethod("InteractionArray", c("ANY", "GInteractions"),
          function(assays, interactions, ...) {
              .newInteractionArray(assays, interactions, ...)
          })

#' Method
setMethod("InteractionArray", c("missing", "missing"),
          function(assays, interactions, ...) {
              .newInteractionArray(list(), GInteractions(), ...)
          })

#' Accessors
#' @param x InteractionArray
#' @noRd
.counts <- function(x) {
    aperm(assay(x, 'counts'), c(3,4,1,2))
}


#' Method
setMethod("counts",
          signature(x="InteractionArray"),
          definition = .counts)

#' Show
.showInteractionArray <- function(object) {
    output <- utils::capture.output(callNextMethod(object))
    if (length(assays(object)) == 0) {
        dims <- rep(0, 4)
    } else {
        dims <- dim(assay(object, 'counts'))
    }
    msg <- "dim: %s interaction(s), %s file(s), %sx%s count matrix(es)"
    output[2] <- do.call(sprintf, c(msg, as.list(dims)))
    cat(output, sep = "\n")
}

#' Show
setMethod("show",
          signature(object="InteractionArray"),
          definition = .showInteractionArray)
