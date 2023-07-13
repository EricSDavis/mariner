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

#' @rdname assignToBins
#' @export
setGeneric("assignToBins",
           function(x,
                    binSize,
                    pos1 = 'center',
                    pos2 = 'center',
                    ...)
               standardGeneric("assignToBins"))

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

#' @rdname clusters
#' @export
setGeneric("clusters", function(x, ...)
    standardGeneric("clusters"))

#' @rdname sources
#' @export
setGeneric("sources", function(x)
    standardGeneric("sources"))

#' @rdname sets
#' @export
setGeneric("sets", function(x,
                                      include,
                                      exclude)
    standardGeneric("sets"))

#' @rdname aggMetadata
#' @export
setGeneric("aggMetadata",
           function(x, columns, funs)
               standardGeneric("aggMetadata"))

#' @rdname pullHicMatrices
#' @export
setGeneric("pullHicMatrices",
           function(x,
                    files,
                    binSize,
                    ...,
                    h5File = tempfile(fileext = ".h5"),
                    half="both",
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
                    files,
                    binSize,
                    ...,
                    h5File = tempfile(fileext = ".h5"),
                    half="both",
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
setGeneric("counts<-", package="BiocGenerics")
setGeneric("rbind", package="BiocGenerics")
setGeneric("cbind", package="BiocGenerics")
setGeneric("path", package="BiocGenerics")
setGeneric("path<-", package="BiocGenerics")


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

#' @rdname plotMatrix
#' @export
setGeneric("plotMatrix",
           function(data,
                    params=NULL,
                    x=NULL,
                    y=NULL,
                    width=NULL,
                    height=NULL,
                    just=c("left", "top"),
                    default.units="inches",
                    draw=TRUE,
                    palette=colorRampPalette(
                        RColorBrewer::brewer.pal(9, 'YlGnBu')
                    ),
                    zrange=NULL,
                    na.color="grey")
           standardGeneric("plotMatrix"))

#' @rdname selection-functions
#' @export
setGeneric("selectRadius", function(x, buffer, invert=FALSE)
    standardGeneric("selectRadius"))

#' @rdname selection-functions
#' @export
setGeneric("selectCenterPixel", function(mhDist, buffer, invert=FALSE)
    standardGeneric("selectCenterPixel"))

#' @rdname selection-functions
#' @export
setGeneric("selectSubmatrix", function(m, invert=FALSE)
    standardGeneric("selectSubmatrix"))

#' @rdname selection-functions
#' @export
setGeneric("selectCoordinates",
           function(rowInd, colInd, buffer, invert=FALSE)
    standardGeneric("selectCoordinates"))

#' @rdname selection-functions
#' @export
setGeneric("selectBlock",
           function(rowInd, colInd, buffer, invert=FALSE)
               standardGeneric("selectBlock"))

#' @rdname selection-functions
#' @export
setGeneric("selectTopLeft",
           function(n, buffer, inset=0, invert=FALSE)
               standardGeneric("selectTopLeft"))

#' @rdname selection-functions
#' @export
setGeneric("selectTopRight",
           function(n, buffer, inset=0, invert=FALSE)
               standardGeneric("selectTopRight"))

#' @rdname selection-functions
#' @export
setGeneric("selectBottomRight",
           function(n, buffer, inset=0, invert=FALSE)
               standardGeneric("selectBottomRight"))

#' @rdname selection-functions
#' @export
setGeneric("selectBottomLeft",
           function(n, buffer, inset=0, invert=FALSE)
               standardGeneric("selectBottomLeft"))

#' @rdname selection-functions
#' @export
setGeneric("selectCorners",
           function(n, buffer, inset=0, invert=FALSE)
               standardGeneric("selectCorners"))

#' @rdname selection-functions
#' @export
setGeneric("selectRows",
           function(rows, buffer, invert=FALSE)
               standardGeneric("selectRows"))

#' @rdname selection-functions
#' @export
setGeneric("selectCols",
           function(cols, buffer, invert=FALSE)
               standardGeneric("selectCols"))

#' @rdname selection-functions
#' @export
setGeneric("selectInner",
           function(n, buffer, invert=FALSE)
               standardGeneric("selectInner"))

#' @rdname selection-functions
#' @export
setGeneric("selectOuter",
           function(n, buffer, invert=FALSE)
               standardGeneric("selectOuter"))

#' @rdname calcLoopEnrichment
#' @export
setGeneric("calcLoopEnrichment",
           function(x, files,
                    fg=selectCenterPixel(mhDist=1, buffer=defaultBuffer()),
                    bg=selectTopLeft(n=4, buffer=defaultBuffer()) +
                        selectBottomRight(n=4, buffer=defaultBuffer()),
                    FUN=\(fg, bg) median(fg+1) / median(bg+1),
                    nBlocks=5,
                    verbose=TRUE,
                    BPPARAM=bpparam(),
                    ...)
               standardGeneric("calcLoopEnrichment"))

#' @rdname adjustEnrichment
#' @export
setGeneric("plotEnrichment",
           function(scores, interactions, k=25, nknots=10, plot=TRUE)
               standardGeneric("plotEnrichment"))

#' @rdname adjustEnrichment
#' @export
setGeneric("adjustEnrichment",
           function(x, interactions, k=25, nknots=10)
               standardGeneric("adjustEnrichment"))

#' @rdname removeShortPairs
#' @export
setGeneric("removeShortPairs",
           function(x, padding=0)
               standardGeneric("removeShortPairs"))

#' @rdname makeRandomGInteractions
#' @export
setGeneric("makeRandomGRanges",
           function(seqinfo,
                    n=100,
                    ...)
               standardGeneric("makeRandomGRanges"))

#' @rdname makeRandomGInteractions
#' @export
setGeneric("makeRandomGInteractions",
           function(seqinfo,
                    n=100,
                    interchromosomal=TRUE,
                    ...)
               standardGeneric("makeRandomGInteractions"))

#' @rdname regularize
#' @export
setGeneric("regularize",
           function(x,
                    ndim=c(10,10),
                    h5File=tempfile(fileext=".h5"),
                    scale=TRUE,
                    nBlocks=5,
                    verbose=TRUE,
                    chunkSize=1,
                    compressionLevel=0,
                    ...)
               standardGeneric("regularize"))

#' @rdname pileupPixels
#' @export
setGeneric("pileupPixels",
           function(x,
                    files,
                    binSize,
                    buffer=5,
                    removeShort=TRUE,
                    minPairDist=0,
                    normalize=TRUE,
                    FUN=sum,
                    nBlocks=5,
                    verbose=TRUE,
                    BPPARAM=bpparam(),
                    ...)
               standardGeneric("pileupPixels"))

#' @rdname pileupDomains
#' @export
setGeneric("pileupDomains",
           function(x,
                    files,
                    binSize,
                    buffer=0.5,
                    ndim=c(100, 100),
                    scale=TRUE,
                    normalize=TRUE,
                    FUN=sum,
                    nBlocks=50,
                    verbose=TRUE,
                    BPPARAM=bpparam(),
                    blockSize=1e6,
                    ...)
           standardGeneric("pileupDomains"))
