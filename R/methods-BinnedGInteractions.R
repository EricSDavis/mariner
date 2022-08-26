## Accessors -------------------------------------------------------------------

#' Accessor for firstBinSize
#' @param x BinnedGInteractions object
#' @rdname BinnedGInteractions-class
#' @export
setMethod("firstBinSize", "BinnedGInteractions",
          function(x) x@firstBinSize)

#' Accessor for secondBinSize
#' @param x BinnedGInteractions object
#' @rdname BinnedGInteractions-class
#' @export
setMethod("secondBinSize", "BinnedGInteractions",
          function(x) x@secondBinSize)

#' Accessor for pairBinsEqual
#' @param x BinnedGInteractions object
#' @rdname BinnedGInteractions-class
#' @export
setMethod("pairBinsEqual", "BinnedGInteractions",
          function(x) x@pairBinsEqual)

## Coersion --------------------------------------------------------------------

#' Coerce GInteractions to BinnedGInteractions
#' @name BinnedGInteractions
#' @rdname BinnedGInteractions-class
setAs(from = "GInteractions", to = "BinnedGInteractions",
      def = function(from) new("BinnedGInteractions", delegate = from))

## Set intra-range methods that adjust binSize ---------------------------------

#' trim method for BinnedGInteractions
#' @param x BinnedGInteractions object
#' @param use.names See \code{? \link[IRanges]{intra-range-methods}}
#' @importFrom GenomicRanges trim
#' @rdname BinnedGInteractions-trim
#' @export
setMethod("trim", "BinnedGInteractions",
          function(x, use.names=TRUE) {
              x <- callNextMethod(x, use.names=use.names)
              new("BinnedGInteractions", delegate = x)
          })

#' resize method for BinnedGInteractions
#' @param x BinnedGInteractions object
#' @param width,fix,use.names See
#'   \code{? \link[IRanges]{intra-range-methods}}
#' @param ... Additional arguments
#' @importFrom GenomicRanges resize
#' @rdname BinnedGInteractions-resize
#' @export
setMethod("resize", "BinnedGInteractions",
          function(x, width, fix="start", use.names=TRUE, ...) {
              x <- callNextMethod(x, width, fix=fix, use.names=use.names, ...)
              new("BinnedGInteractions", delegate = x)
          })

#' narrow method for BinnedGInteractions
#' @param x BinnedGInteractions object
#' @param start,end,width,use.names See
#'   \code{? \link[IRanges]{intra-range-methods}}
#' @importFrom GenomicRanges narrow
#' @rdname BinnedGInteractions-narrow
#' @export
setMethod("narrow", "BinnedGInteractions",
          function(x, start=NA, end=NA, width=NA, use.names=TRUE) {
              x <- callNextMethod(x, start=start, end=end,
                                  width=width, use.names=use.names)
              new("BinnedGInteractions", delegate = x)
          })

#' shift method for BinnedGInteractions
#' @param x BinnedGInteractions object
#' @param shift,use.names See
#'   \code{? \link[IRanges]{intra-range-methods}}
#' @importFrom GenomicRanges shift
#' @rdname BinnedGInteractions-shift
#' @export
setMethod("shift", "BinnedGInteractions",
          function(x, shift=0L, use.names=TRUE) {
              x <- callNextMethod(x, shift=shift, use.names=use.names)
              new("BinnedGInteractions", delegate = x)
          })

#' flank method for BinnedGInteractions
#' @param x BinnedGInteractions object
#' @param width,start,both,use.names See
#'   \code{? \link[IRanges]{intra-range-methods}}
#' @param ignore.strand See
#'   \code{? \link[GenomicRanges]{flank}}
#' @importFrom GenomicRanges flank
#' @rdname BinnedGInteractions-flank
#' @export
setMethod("flank", "BinnedGInteractions",
          function(x, width, start=TRUE, both=FALSE,
                   use.names=TRUE, ignore.strand=FALSE) {
              x <- callNextMethod(x, width, start=start, both=both,
                                  use.names=use.names,
                                  ignore.strand=ignore.strand)
              new("BinnedGInteractions", delegate = x)
          })

#' `regions<-` method for BinnedGInteractions
#' @param x BinnedGInteractions object
#' @param value See
#'   \code{? \link[InteractionSet]{interaction-access}}
#' @importFrom InteractionSet `regions<-`
#' @rdname BinnedGInteractions-regions
#' @export
setReplaceMethod("regions",
                 "BinnedGInteractions",
                 function(x, value) {
                     x <- callNextMethod(x, value=value)
                     new("BinnedGInteractions", delegate = x)
                 })

#' swapAnchors method for BinnedGInteractions
#' @param x BinnedGInteractions object
#' @param ... See
#'   \code{? \link[InteractionSet]{swapAnchors}}
#' @importFrom InteractionSet swapAnchors
#' @rdname BinnedGInteractions-swapAnchors
#' @export
setMethod("swapAnchors",
          "BinnedGInteractions",
          function(x, ...) {
              x <- callNextMethod(x, ...)
              new("BinnedGInteractions", delegate = x)
          })
