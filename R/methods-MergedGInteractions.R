## Accessors -------------------------------------------------------------------

#' Sources
#'

#' Internal accessor function for allPairs
#' @inheritParams allPairs
#' @importFrom S4Vectors isEmpty
#' @importFrom glue glue
#' @importFrom rlang abort
.allPairs <- function(x, source) {

    ## Separate solo loops as negative integers
    ap <- x@allPairs
    ap[clst == 0, clst := seq(1, length(clst))*-1, by = grp]

    ## Create a key using ids from merged pairs
    keys <- ap[id %in% x@ids, .(id, grp, clst)]

    ## Join keys by group and cluster
    res <- ap[keys, on = .(grp, clst)]

    ## Save the split vector & drop extra columns
    s <- res$i.id
    res <- res[, -c("id", "grp", "clst", "i.id")]

    ## Split and order result by input vector (slowest step)
    res <-
        split(res, s)[rank(x@ids)] |>
        `names<-`(values = names(x))

    return(res)
}

#' Get all pairs from MergedGInteractions object
#' @inheritParams selectionMethod
#' @export
setMethod("allPairs", signature(x = "MergedGInteractions"),
          definition = .allPairs)

#' Get selectionMethod from MergedGInteractions object
#'
#' @param x MergedGInteractions object.
#' @param ... Additional arguments.
#'
#' @return A character vector describing which selection
#'  method was used for merging.
#'
#' @examples
#' ## Reference BEDPE files (loops called with SIP)
#' bedpeFiles <-
#'     system.file("extdata", package = "mariner") |>
#'     list.files(pattern = "Loops.txt", full.names = TRUE)
#'
#' x <- mergePairs(x = bedpeFiles,
#'                 binSize = 5e03,
#'                 radius = 2,
#'                 column = "APScoreAvg")
#' selectionMethod(x)
#'
#' @rdname selectionMethod
#' @export
setMethod("selectionMethod", "MergedGInteractions", function(x, ...) {
    x@selectionMethod
})
