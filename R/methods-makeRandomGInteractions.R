#' Creating random GRanges
#' @inheritParams makeRandomGInteractions
#' @param .rows Vector of row positions to seqinfo
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom SummarizedExperiment seqnames
#' @importFrom data.table data.table .SD
#' @returns List of GRanges and vector of row positions.
#' @noRd
.makeRandomGRanges <- function(seqinfo, n, .rows=NULL) {

    ## Suppress NSE notes in R CMD check
    chrom <- NULL

    ## Convert Seqinfo to a chromSizes
    ## data.table to avoid rowname issues
    cs <- data.table(
        chrom=seqnames(seqinfo),
        length=seqlengths(seqinfo)
    )

    ## Sample rows of seqinfo or use provided rows
    if (is.null(.rows)) {
        rows <- sample(seq_len(nrow(cs)), n, replace=TRUE)
    } else {
        stopifnot(length(.rows) == n) #internal catch
        rows <- .rows
    }

    ## Sample ranges per chromosome
    ans <- cs[rows][, {
        l <- .SD$length[1]
        s1 <- sample(l, .N, replace=TRUE)
        s2 <- sample(l, .N, replace=TRUE)
        s3 <- pmin(s1, s2)
        e <- s3 + abs(s1-s2)
        list(start=s3, end=e)
    }, by=chrom]

    ## Convert to GRanges object & return as list
    ## with the sampled rows
    gr <- makeGRangesFromDataFrame(ans, seqinfo=seqinfo)
    return(list(gr=gr, rows=rows))

}

#' Creating random GInteractions
#' @inheritParams makeRandomGInteractions
#' @returns A GInteractions object.
#' @noRd
.makeRandomGInteractions <- function(seqinfo, n, interchromosomal) {

    ## Returns list of GRanges & rows
    anchor1 <- .makeRandomGRanges(seqinfo, n)

    ## Allow interchromosomal by sampling again
    if (interchromosomal) {
        gi <- GInteractions(
            anchor1$gr,
            makeRandomGRanges(seqinfo, n) #returns GRanges
        )
        return(gi)
    }

    ## Otherwise use the rows from anchor1 to sample
    ## on the same chromosome
    anchor2 <- .makeRandomGRanges(seqinfo, n, .rows=anchor1$rows)
    gi <- GInteractions(
        anchor1$gr,
        anchor2$gr
    )
    return(gi)
}

#' @examples
#' ## Define Seqinfo containing chromosome info
#' if (require(TxDb.Hsapiens.UCSC.hg38.knownGene)) {
#'     txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#'     si <- seqinfo(txdb)
#'     si <- keepStandardChromosomes(si)
#' } else {
#'     si <- Seqinfo(
#'         seqnames=c("chr1", "chr2"),
#'         seqlengths=rep(200e6, 2),
#'         genome="hg38"
#'     )
#' }
#'
#' ## Make some GRanges
#' set.seed(123)
#' makeRandomGRanges(si, 100)
#'
#' @rdname makeRandomGInteractions
#' @export
setMethod("makeRandomGRanges", signature(seqinfo="Seqinfo"),
          function(seqinfo, n, .rows=NULL) {
              .makeRandomGRanges(seqinfo, n, .rows)$gr
          })

#' Creating random GRanges & GInteractions
#' @param seqinfo A Seqinfo object containing the
#'  chromosome names, lengths, and genome build.
#' @param n Integer describing the number of
#'  random sequences to generate
#' @param interchromosomal Boolean (TRUE/FALSE) indicating
#'  whether interchromosomal interactions should be allowed.
#'  Default is TRUE.
#' @param ... Additional arguments.
#' @param .rows (internal use only) vector of row positions
#'  to sample from seqinfo.
#' @returns A GRanges or GInteractions object with ranges
#'  selected randomly with replacement on the provided
#'  seqinfo.
#'
#' @examples
#' ## Make some GInteractions
#' set.seed(123)
#' makeRandomGInteractions(si, n=100)
#'
#' ## Make some GInteractions only on same chromosome
#' set.seed(123)
#' makeRandomGInteractions(si, n=100, interchromosomal=FALSE)
#'
#' ## Use specific binSizes
#' n <- 100
#' binOptions <- seq(5e3, 200e3, by=5e3)
#' si <- Seqinfo(seqnames="chr1", seqlengths=200e6, genome="hg38")
#' set.seed(123)
#' bins <- sample(binOptions, n, replace=TRUE)
#' makeRandomGInteractions(si, n) |>
#'     resize(bins) |>
#'     trim()
#'
#' @rdname makeRandomGInteractions
#' @export
setMethod("makeRandomGInteractions", signature(seqinfo="Seqinfo"),
          definition=.makeRandomGInteractions)


