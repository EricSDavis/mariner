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
                    radius,
                    method = "manhattan",
                    column,
                    selectMax = TRUE,
                    pos = "center")
               standardGeneric("mergePairs"))

#' @rdname selectionMethod
#' @export
setGeneric("selectionMethod", function(x, ...)
    standardGeneric("selectionMethod"))

#' @rdname getPairClusters
#' @export
setGeneric("getPairClusters", function(x, ...)
    standardGeneric("getPairClusters"))

#' @rdname sources
#' @export
setGeneric("sources", function(x, ...)
    standardGeneric("sources"))

#' @rdname subsetBySource
#' @export
setGeneric("subsetBySource", function(x,
                                      include,
                                      exclude)
    standardGeneric("subsetBySource"))

#' @rdname aggPairMcols
#' @export
setGeneric("aggPairMcols",
           function(x, columns, funs)
               standardGeneric("aggPairMcols"))

#' @rdname BinnedGInteractions-class
#' @export
setGeneric("firstBinSize", function(x) standardGeneric("firstBinSize"))

#' @rdname BinnedGInteractions-class
#' @export
setGeneric("secondBinSize", function(x) standardGeneric("secondBinSize"))

#' @rdname BinnedGInteractions-class
#' @export
setGeneric("pairBinsEqual", function(x) standardGeneric("pairBinsEqual"))

#' @rdname pullHicPixels
#' @export
setGeneric("pullHicPixels",
           function(x,
                    hic,
                    binSize,
                    norm = "NONE",
                    ...)
               standardGeneric("pullHicPixels"))
