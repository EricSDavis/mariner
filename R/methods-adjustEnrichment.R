#' Calculate rolling windows of loop enrichment
#' @param x A GInteractions Object.
#' @param scores Numeric vector of enrichment scores.
#' @param k Number of observations for rolling window.
#' @param thresh Numeric loop size cutoff (keeps loops)
#'  less than this size.
#' @importFrom InteractionSet pairdist
#' @importFrom data.table frollapply data.table
#' @importFrom stats median
#' @returns A data table with median loop size and enrichment
#'  after calculating the rolling median within a window.
#' @noRd
.rollEnrich <- function(x, scores, k=200, thresh=Inf) {

    ## Suppress NSE notes in R CMD check
    size <- rollMedSize <- rollMedScore <- NULL

    ## Data table of loop size and score
    dat <- data.table(
        size = pairdist(x),
        scores = scores
    )

    ## Order then filter out interactions >thresh
    dat <- dat[order(size)]
    dat <- dat[size <= thresh]

    ## Calculate rollMedians for size and score
    FUN <- \(x) median(x)
    dat[, rollMedSize := frollapply(x=size, n=k, FUN=FUN)]
    dat[, rollMedScore := frollapply(x=scores, n=k, FUN=FUN)]

    ## Take the median per group
    ans <- dat[, .(rollMedScore = median(rollMedScore)), by=rollMedSize]
    ans
}

#' Define function to correct the effect of distance
#' on loop enrichment
#' @param x A GInteractions Object.
#' @param scores Numeric vector of enrichment scores.
#' @param k Number of observations for rolling window.
#' @param nknots integer or function giving the number
#'  of knots to use see `?smooth.spline` for more info.
#' @importFrom stats smooth.spline predict median
#' @importFrom InteractionSet pairdist
#' @returns Vector of corrected loop enrichment scores.
#' @noRd
.adjustEnrich <- function(x, scores, k=25, nknots=10) {
    ## Calculate rolling enrichment
    re <- .rollEnrich(x, scores=scores, k=k)

    ## Fit smoothed spline
    sp <- smooth.spline(na.omit(re), nknots=nknots)
    m <- median(na.omit(re$rollMedScore))

    ## Calculate correction factor
    ratio <- m / predict(sp, pairdist(x))$y

    ## Adjusted scores
    adjusted <- scores*ratio

    ## Return adjusted scores
    return(list(
        rollEnrich = re,
        splineFit = sp,
        corrFactor = ratio,
        adjusted = adjusted
    ))
}

#' Adjust loop enrichment to remove distance-
#' dependent effect.
#' @inheritParams adjustEnrichment
#' @importFrom DelayedArray DelayedArray
#' @noRd
.adjustEnrichment <- function(x, interactions, k, nknots) {
    ## Apply enrichment to each column of a matrix
    ans <- vapply(seq_len(ncol(x)), \(i) {
        ae <- .adjustEnrich(x=interactions, scores=x[,i], k, nknots)
        ae$adjusted
    }, numeric(nrow(x))) |>
        DelayedArray() |>
        `colnames<-`(colnames(x))
    ans
}

#' Adjust loop enrichment to remove distance-
#' dependent effect.
#'
#' @param x A DelayedMatrix or matrix with enrichment
#'  scores.
#' @param interactions A GInteractions Object containing
#'  the interactions used to calculate enrichment scores.
#' @param k Number of observations for rolling window.
#' @param nknots integer or function giving the number
#'  of knots to use see `?smooth.spline` for more info.
#' @returns A DelayedMatrix of enrichment scores where rows
#'  are loops and columns are Hi-C files.
#' @examples
#' ## Load marinerData
#' if (!require("marinerData", quietly = TRUE))
#'     BiocManager::install("marinerData")
#'
#' ## Read .hic file paths
#' hicFiles <- c(
#'     marinerData::LEUK_HEK_PJA27_inter_30.hic(),
#'     marinerData::LEUK_HEK_PJA30_inter_30.hic()
#' )
#' names(hicFiles) <- c("FS", "WT")
#'
#' ## Read in loops as GInteractions object
#' loops <-
#'     WT_5kbLoops.txt() |>
#'     setNames("WT") |>
#'     read.table(header=TRUE, nrows=1000) |>
#'     as_ginteractions(keep.extra.columns=FALSE)
#'
#' ## Removes the "chr" prefix for compatibility
#' ## with the preprocessed hic files
#' GenomeInfoDb::seqlevelsStyle(loops) <- 'ENSEMBL'
#'
#' ## Calculate loop enrichment
#' enrich <- calcLoopEnrichment(
#'     x=assignToBins(loops, 100e03),
#'     files=hicFiles
#' )
#'
#' adjustEnrichment(enrich, loops)
#'
#' @rdname adjustEnrichment
#' @export
setMethod("adjustEnrichment",
          signature(x="DelayedMatrix_OR_matrix",
                    interactions="GInteractions"),
          .adjustEnrichment)

#' Show diagnostic plot of enrichment before
#' and after correction.
#' @inheritParams plotEnrichment
#' @importFrom graphics lines
#' @importFrom graphics abline
#' @importFrom graphics legend
#' @noRd
.plotEnrich <- function(scores, interactions, k, nknots, plot) {

    ## Suppress NSE notes in R CMD check
    size <- rollMedSize <- rollMedScore <- NULL

    ## Adjust enrichment scores & assign results
    ## to variables for plotting
    ae <- .adjustEnrich(interactions, scores, k, nknots)
    re <- ae$rollEnrich
    sp <- ae$splineFit
    m  <- median(na.omit(re$rollMedScore))
    ratio <- ae$corrFactor
    adjusted <- ae$adjusted

    ## Diagnostic plot
    if (plot) {
        plot(re$rollMedSize, re$rollMedScore, type='l',
             main=paste0("k=", k, ", nknots=", nknots), col='lightgrey')
        lines(sp$x, sp$y, col="forestgreen")
        abline(h=m, lty=2)
        cre <- .rollEnrich(interactions, scores=adjusted, k=k)
        lines(cre$rollMedSize, cre$rollMedScore, col="blue")
        legend(x='bottomright',
               legend=c("Before adjustment",
                        "Smooth spline fitted",
                        "After adjustment"),
               fill=c("lightgrey", "forestgreen", "blue"))
    }

    return(ae)
}

#' Show diagnostic plot of loop enrichment before
#' and after distance adjustment.
#' @param interactions A GInteractions Object containing
#'  the interactions used to calculate enrichment scores.
#' @param scores Numeric vector of enrichment scores.
#' @param k Number of observations for rolling window.
#' @param nknots integer or function giving the number
#'  of knots to use see `?smooth.spline` for more info.
#' @param plot Boolean (default=FALSE), of whether to
#'  show diagnostic plot.
#' @returns A plot (and associated data) for visualizing
#'  loop enrichment before and after distance adjustment.
#' @examples
#'
#' plotEnrichment(enrich[,1], loops)
#'
#' @rdname adjustEnrichment
#' @export
setMethod("plotEnrichment",
          signature(scores="numeric",
                    interactions="GInteractions"),
          .plotEnrich)
