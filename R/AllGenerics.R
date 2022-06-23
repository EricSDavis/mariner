#' @rdname as_ginteractions
#' @export
setGeneric("as_ginteractions",
           function(df,
                    keep.extra.columns = TRUE,
                    starts.in.df.are.0based = FALSE,
                    ...)
               standardGeneric("as_ginteractions"))

#' @rdname as_ginteractions
#' @export
setGeneric("makeGInteractionsFromDataFrame",
           function(df,
                    keep.extra.columns = TRUE,
                    starts.in.df.are.0based = FALSE,
                    ...)
               standardGeneric("makeGInteractionsFromDataFrame"))

#' @rdname binPairs
#' @export
setGeneric("binPairs",
           function(x,
                    binSize,
                    pos1 = 'center',
                    pos2 = 'center',
                    ...)
               standardGeneric("binPairs"))

#' @rdname binRanges
#' @export
setGeneric("binRanges",
          function(x,
                   binSize,
                   pos = 'center')
              standardGeneric("binRanges"))

#' @rdname shiftRanges
#' @export
setGeneric("shiftRanges",
           function(x, pos)
              standardGeneric("shiftRanges"))

#' @rdname mergePairs
#' @export
setGeneric("mergePairs",
           function(x,
                    binSize = 5000,
                    radius = 2,
                    column)
               standardGeneric("mergePairs"))

#' @rdname selectionMethod
#' @export
setGeneric("selectionMethod", function(x, ...)
    standardGeneric("selectionMethod"))

#' @rdname allPairs
#' @export
setGeneric("allPairs", function(x, ...)
    standardGeneric("allPairs"))

#' @rdname deNovo
#' @export
setGeneric("deNovo", function(x, ...)
    standardGeneric("deNovo"))

#' @rdname aggPairMcols
#' @export
setGeneric("aggPairMcols",
           function(x, columns, funs)
               standardGeneric("aggPairMcols"))
