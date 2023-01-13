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

#' @rdname snapToBins
#' @export
setGeneric("snapToBins",
           function(x,
                    binSize)
               standardGeneric("snapToBins"))

#' @rdname mergePairs
#' @export
setGeneric("mergePairs",
           function(x,
                    radius,
                    method = "manhattan",
                    column = NULL,
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
setGeneric("sources", function(x)
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

#' @rdname pullHicMatrices
#' @export
setGeneric("pullHicMatrices",
           function(x,
                    binSize,
                    files,
                    ...,
                    half="both",
                    h5File = tempfile(fileext = ".h5"),
                    norm = 'NONE',
                    matrix = 'observed',
                    blockSize = 248956422,
                    onDisk = TRUE,
                    compressionLevel = 0,
                    chunkSize = 1)
               standardGeneric("pullHicMatrices"))

#' @rdname pullHicPixels
#' @export
setGeneric("pullHicPixels",
           function(x,
                    binSize,
                    files,
                    ...,
                    half="both",
                    h5File = tempfile(fileext = ".h5"),
                    norm = 'NONE',
                    matrix = 'observed',
                    blockSize = 248956422,
                    onDisk = TRUE,
                    compressionLevel = 0,
                    chunkSize = 1)
               standardGeneric("pullHicPixels"))

#' @rdname pixelsToMatrices
#' @export
setGeneric("pixelsToMatrices", function(x, buffer)
    standardGeneric("pixelsToMatrices"))

#' @rdname InteractionArray-class
#' @export
setGeneric("InteractionArray", function(assays, interactions, ...)
    standardGeneric("InteractionArray"))

setGeneric("counts", package="BiocGenerics")
setGeneric("rbind", package="BiocGenerics")
setGeneric("cbind", package="BiocGenerics")


#' @rdname InteractionMatrix-class
#' @export
setGeneric("InteractionMatrix", function(assays, interactions, ...)
    standardGeneric("InteractionMatrix"))

#' @rdname hdf5BlockApply
#' @export
setGeneric("hdf5BlockApply", function(x, FUN, sink, grid, sink_grid,
                                      verbose=TRUE)
    standardGeneric("hdf5BlockApply"))

#' @rdname aggHicMatrices
#' @export
setGeneric("aggHicMatrices", function(x,
                                      by=NULL,
                                      FUN=sum,
                                      nBlocks=5,
                                      verbose=TRUE,
                                      BPPARAM=bpparam(),
                                      compressionLevel=0)
    standardGeneric("aggHicMatrices"))

#' @rdname selectPixel
#' @export
setGeneric("selectPixel", function(x,
                                   aggFUN=sum,
                                   selectFUN="which.max",
                                   nBlocks=5,
                                   verbose=TRUE)
    standardGeneric("selectPixel"))

#' @rdname changePixelRes
#' @export
setGeneric("changePixelRes", function(x, files, from, to,
                                      aggFUN=sum,
                                      selectFUN="which.max",
                                      nBlocks=5,
                                      verbose=TRUE,
                                      norm="KR",
                                      half="upper",
                                      ...)
    standardGeneric("changePixelRes"))

#' @rdname GInteractions-accessors
#' @export
setGeneric("seqnames1", function(x, ...) standardGeneric("seqnames1"))

#' @rdname GInteractions-accessors
#' @export
setGeneric("seqnames2", function(x, ...) standardGeneric("seqnames2"))

#' @rdname GInteractions-accessors
#' @export
setGeneric("start1", function(x, ...) standardGeneric("start1"))

#' @rdname GInteractions-accessors
#' @export
setGeneric("end1", function(x, ...) standardGeneric("end1"))

#' @rdname GInteractions-accessors
#' @export
setGeneric("start2", function(x, ...) standardGeneric("start2"))

#' @rdname GInteractions-accessors
#' @export
setGeneric("end2", function(x, ...) standardGeneric("end2"))

#' @rdname calcLoopEnrichment
#' @export
setGeneric("calcLoopEnrichment",
           function(x, files,
                    mhDist=c(4,5,6),
                    nBlocks=5,
                    verbose=TRUE,
                    BPPARAM=bpparam(),
                    ...)
    standardGeneric("calcLoopEnrichment"))
