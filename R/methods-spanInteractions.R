#' Internal function taken from GenomicRanges package
#' Returns index of out-of-bound ranges located on non-circular sequences
#' whose length is not NA. Works on a GenomicRanges or GAlignments object.
#' @param x GenomicRanges or GAlignments object.
#' @importFrom GenomeInfoDb isCircular
#' @noRd
.get_out_of_bound_index <- function(x)
{
    if (length(x) == 0L)
        return(integer(0))
    x_seqnames_id <- as.integer(seqnames(x))
    x_seqlengths <- unname(seqlengths(x))
    seqlevel_is_circ <- unname(isCircular(x)) %in% TRUE
    seqlength_is_na <- is.na(x_seqlengths)
    seqlevel_has_bounds <- !(seqlevel_is_circ | seqlength_is_na)
    which(seqlevel_has_bounds[x_seqnames_id] &
              (start(x) < 1L | end(x) > x_seqlengths[x_seqnames_id]))
}

#' Internal spanInteractions
#' @inheritParams spanInteractions
#' @importFrom InteractionSet pairdist intrachr
#' @importFrom rlang warn
#' @importFrom GenomeInfoDb seqinfo Seqinfo `seqinfo<-`
#' @noRd
.spanInteractions <- function(x, buffer) {
    
    ## Check inputs
    .checkTypes(c(buffer="number"))
    if (buffer < 0) abort("Buffer must be >= 0")
    
    ## Keep track of original x indices
    xIndices <- seq_len(length(x))
    
    ## Remove interchromosomal & warn
    if (any(!intrachr(x))) {
        interchr <- which(!intrachr(x))
        warn(c("Interchromosomal interactions are not allowed.",
               "i"="The following indices in `x` have been removed:",
               paste0(interchr, collapse=",")))
        xIndices <- xIndices[intrachr(x)]
        x <- x[intrachr(x)]
    }
    
    ## Get expansion factor
    ifelse(buffer == 0, d <- 0, d <- pairdist(x, 'span')/(1/buffer))
    
    ## Construct ranges
    df <- data.frame(
        seqnames=seqnames1(x),
        start=start1(x)-d,
        end=end2(x)+d
    )
    suppressWarnings({
        gr <- makeGRangesFromDataFrame(
            df=df,
            keep.extra.columns=TRUE,
            seqinf=seqinfo(x)
        )
    })
    
    ## Remove out of bound ranges & warn
    oob <- .get_out_of_bound_index(gr)
    if (length(oob) > 0) {
        oob_x <- xIndices[oob]
        warn(c("Some ranges are out of bounds.",
               "i"="The following indices in `x` have been removed:",
               paste0(oob_x, collapse=",")))
        gr <- gr[-oob]
    }
    
    return(gr)
    
}

#' Span the distance of an interaction with an
#' expansion buffer.
#' 
#' Takes the span (start1 and end2) of the interaction
#' and adds an expansion buffer that is a fraction
#' of the interaction distance.
#'
#' Note, this function removes interchromosomal and
#' ranges that would fall out of bounds of the
#' sequence lengths.
#'
#' @param x GInteractions object.
#' @param buffer Fraction (length one numeric vector)
#'  pair-distance to expand around the resulting
#'  range.
#'
#' @returns A GRanges object.
#'
#' @examples
#' ## Example GInteractions
#' gi <- as_ginteractions(read.table(text="
#'         1 10 40 1 50 80
#'         2 100 110 2 150 160
#'     "))
#' 
#' ## Add seqinfo
#' GenomeInfoDb::seqinfo(gi) <- 
#'     GenomeInfoDb::Seqinfo(
#'         seqnames=c("1", "2"),
#'         seqlengths=c(110, 200)
#'     )
#' 
#' ## Get span with an expansion buffer
#' spanInteractions(gi, buffer=0.1)
#' 
#' @rdname spanInteractions
#' @export
setMethod("spanInteractions",
          signature(x='GInteractions'),
          definition = .spanInteractions)
