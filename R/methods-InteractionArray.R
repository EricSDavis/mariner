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
#' @importFrom DelayedArray DelayedArray
#' @noRd
.counts <- function(x, showDimnames) {

    ## Get the full object
    cnts <- aperm(assay(x, 'counts'), c(3,4,1,2))

    if (showDimnames) {
        ## Row/colnames
        rows <- assay(x, 'rownames')
        cols <- assay(x, 'colnames')

        ## Get the header and dimensions
        fullHeader <- capture.output(show(cnts))[1]
        dims <- dim(cnts)
        dn <- dimnames(cnts)

        ## Fxn to add row/colnames to matrix
        mat_with_dimnames <- function(x, rows, cols) {
            if (is.vector(x)) {
                x <- matrix(x, length(rows), length(cols))
                x <- DelayedArray(x)
            }
            dimnames(x) <- list(rows, cols)
            capture.output(show(x))[-1]
        }

        ## Construct first & last matrix header
        firstHeader <-
            paste0(c(",,"),
                   ifelse(is.null(dn[[3]]), 1, dn[[1]][1]), ",",
                   ifelse(is.null(dn[[4]]), 1, dn[[4]][1]))

        lastHeader <-
            paste0(c(",,"),
                   ifelse(is.null(dn[[3]]),
                          dims[3], dn[[1]][dims[3]]), ",",
                   ifelse(is.null(dn[[4]]),
                          dims[4], dn[[4]][dims[4]]))

        ## Construct new object to print
        ans <- fullHeader
        ans <- c(ans, firstHeader)
        ans <- c(ans, mat_with_dimnames(cnts[,,1,1],
                                        rows[1,1,], cols[1,1,]))
        if (dims[3] != 1L | dims[4] != 1L) {
            ans <- c(ans, "\n...\n")
            ans <- c(ans, lastHeader)
            ans <- c(ans, mat_with_dimnames(cnts[,,dims[3],dims[4]],
                                            rows[dims[3],dims[4],],
                                            cols[dims[3],dims[4],]))
        } else {
            dimnames(cnts)[c(1,2)] <- list(rows[1,1,], cols[1,1,])
        }
        cat(ans, sep="\n")
        return(invisible(cnts))
    } else {
        return(cnts)
    }

}


#' Method
setMethod("counts",
          signature(x="InteractionArray"),
          definition = .counts)

#' Show
#' @noRd
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
