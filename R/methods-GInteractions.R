## Accessors for GInteractions & InteractionSet
## Classes and their derivatives

#' Access each portion of a GInteractions-like object
#' @param x GInteractions object.
#' @param ... Additional arguments.
#' @importFrom S4Vectors first
#' @importFrom GenomicRanges seqnames
#' @returns A vector of values corresponding to the
#'  requested component of a GInteractions-like object.
#'  For seqnames1 and seqnames2 the RLE is coerced to a
#'  character vector.
#'
#' @examples
#' library(InteractionSet)
#' ## Create example reference interactions objects
#' gi <- read.table(text="
#'     chr1 10 20 chr1 50 60
#'     chr2 30 40 chr2 60 70
#'     chr1 50 60 chr3 10 20") |>
#'     as_ginteractions()
#'
#' iset <- InteractionSet(assays=matrix(nrow=3),
#'                        interactions=gi)
#'
#' ## Access vectors of values
#' seqnames1(gi)
#' start1(gi)
#' end1(gi)
#' seqnames2(gi)
#' start2(gi)
#' end2(gi)
#'
#' ## Also works for InteractionSet-like objects
#' seqnames1(iset)
#' start1(iset)
#' end1(iset)
#' seqnames2(iset)
#' start2(iset)
#' end2(iset)
#'
#' @rdname GInteractions-accessors
#' @export
setMethod("seqnames1", "GInteractions_OR_InteractionSet",
          function(x) as.character(seqnames(first(x))))

#' @importFrom S4Vectors second
#' @importFrom GenomicRanges seqnames
#' @rdname GInteractions-accessors
#' @export
setMethod("seqnames2", "GInteractions_OR_InteractionSet",
          function(x) as.character(seqnames(second(x))))

#' @importFrom S4Vectors first
#' @importFrom GenomicRanges start
#' @rdname GInteractions-accessors
#' @export
setMethod("start1", "GInteractions_OR_InteractionSet",
          function(x) start(first(x)))

#' @importFrom S4Vectors first
#' @importFrom GenomicRanges end
#' @rdname GInteractions-accessors
#' @export
setMethod("end1", "GInteractions_OR_InteractionSet",
          function(x) end(first(x)))

#' @importFrom S4Vectors second
#' @importFrom GenomicRanges start
#' @rdname GInteractions-accessors
#' @export
setMethod("start2", "GInteractions_OR_InteractionSet",
          function(x) start(second(x)))

#' @importFrom S4Vectors second
#' @importFrom GenomicRanges end
#' @rdname GInteractions-accessors
#' @export
setMethod("end2", "GInteractions_OR_InteractionSet",
          function(x) end(second(x)))
