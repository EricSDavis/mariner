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
#' @noRd
.counts <- function(x, showDimnames) {

    ## Get the full object
    if (length(assays(x)) == 0) {
        abort("`x` has no count matrices.")
    }
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


#' Access count matrices from
#' InteractionArray or InteractionMatrix
#'
#' @param x InteractionArray or InteractionMatrix object.
#' @param showDimnames Logical vector of length-one
#'  indicating whether to show dimensions of
#'  count matrices (default TRUE). Only applicable for
#'  InteractionArray objects.
#'
#' @returns For InteractionArray, a 4-dimensional
#'  DelayedArray of Hi-C submatrices is returned with
#'  the following dimensions: rows of count matrix,
#'  columns of count matrix, Interactions in `x`,
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
#' iarr <- pullHicMatrices(x, 5e03, hicFiles)
#'
#' ## Access count matrices
#' counts(iarr)
#' counts(iarr, FALSE)
#'
#' @rdname counts
#' @export
setMethod("counts",
          signature(x="InteractionArray"),
          definition=.counts)

## Show ------------------------------------------------------------------------

#' Internal show method
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

#' show for InteractionArray
#' @rdname InteractionArray-class
#' @export
setMethod("show",
          signature(object="InteractionArray"),
          definition=.showInteractionArray)


## Concatenation ---------------------------------------------------------------

#' rbind InteractionArray
#' @include utils.R
#' @rdname InteractionArray-class
#' @export
setMethod("rbind", "InteractionArray", .rbindIsetDerived)

#' cbind InteractionArray
#' @include utils.R
#' @rdname InteractionArray-class
#' @export
setMethod("cbind", "InteractionArray", .cbindIsetDerived)

