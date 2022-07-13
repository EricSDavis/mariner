## Accessors -------------------------------------------------------------------

#' Accessor for firstBinSize
#' @rdname BinnedGInteractions-class
#' @export
setMethod("firstBinSize", "BinnedGInteractions",
          function(x) x@firstBinSize)

#' Accessor for secondBinSize
#' @rdname BinnedGInteractions-class
#' @export
setMethod("secondBinSize", "BinnedGInteractions",
          function(x) x@secondBinSize)

#' Accessor for pairBinsEqual
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
#' @importFrom GenomicRanges trim
#' @rdname BinnedGInteractions-class
#' @export
setMethod("trim", "BinnedGInteractions",
          function(x, use.names=TRUE) {
              x <- callNextMethod(x, use.names=use.names)
              new("BinnedGInteractions", delegate = x)
          })

#' resize method for BinnedGInteractions
#' @importFrom GenomicRanges resize
#' @rdname BinnedGInteractions-class
#' @export
setMethod("resize", "BinnedGInteractions",
          function(x, width, fix="start", use.names=TRUE, ...) {
              x <- callNextMethod(x, width, fix=fix, use.names=use.names, ...)
              new("BinnedGInteractions", delegate = x)
          })

#' narrow method for BinnedGInteractions
#' @importFrom GenomicRanges narrow
#' @rdname BinnedGInteractions-class
#' @export
setMethod("narrow", "BinnedGInteractions",
          function(x, start=NA, end=NA, width=NA, use.names=TRUE) {
              x <- callNextMethod(x, start=start, end=end,
                                  width=width, use.names=use.names)
              new("BinnedGInteractions", delegate = x)
          })

#' shift method for BinnedGInteractions
#' @importFrom GenomicRanges shift
#' @rdname BinnedGInteractions-class
#' @export
setMethod("shift", "BinnedGInteractions",
          function(x, shift=0L, use.names=TRUE) {
              x <- callNextMethod(x, shift=shift, use.names=use.names)
              new("BinnedGInteractions", delegate = x)
          })

#' flank method for BinnedGInteractions
#' @importFrom GenomicRanges flank
#' @rdname BinnedGInteractions-class
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
#' @importFrom InteractionSet `regions<-`
#' @rdname BinnedGInteractions-class
#' @export
setReplaceMethod("regions",
                 "BinnedGInteractions",
                 function(x, value) {
                     x <- callNextMethod(x, value=value)
                     new("BinnedGInteractions", delegate = x)
                 })

#' swapAnchors method for BinnedGInteractions
#' @importFrom InteractionSet swapAnchors
#' @rdname BinnedGInteractions-class
#' @export
setMethod("swapAnchors",
          "BinnedGInteractions",
          function(x, ...) {
              x <- callNextMethod(x, ...)
              new("BinnedGInteractions", delegate = x)
          })
