## Accessors -------------------------------------------------------------------

#' Internal accessor function for allPairs
#' @inheritParams allPairs
#' @importFrom S4Vectors isEmpty
#' @importFrom glue glue
#' @importFrom rlang abort
.allPairs <- function(x, source) {

    ## Pull out all pairs
    res <- x@allPairs

    if (!missing(source)) {

        ## Filter by source
        res <- res[src == source]

        ## Catch empty error
        if (isEmpty(res)) {
            msg <- c(glue("Source '{source}' not found."),
                     'i' = glue("Did you mean one of these: ",
                                "{paste0(x = unique(x@allPairs$src),
                                         collapse=', ')}?"))
            abort(msg)
        }
    }

    ## Remove other columns
    res <- res[,-c("src", "id", "grp", "clst")]

    return(as_ginteractions(res))
}

#' Get all pairs from MergedGInteractions object
#' @inheritParams selectionMethod
#' @export
setMethod("allPairs", signature(x = "MergedGInteractions",
                                source = "character_OR_numeric_OR_missing"),
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
