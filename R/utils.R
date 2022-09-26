#' Check that a GRanges object has been snapped
#' to bins
#'
#' @param x GRanges object
#' @param binSize Integer (numeric) describing
#'  the new size of each range.
#'
#' @return Logical
#'
#' @importFrom GenomicRanges start end
#'
#' @noRd
.checkSnappedRanges <- function(x, binSize) {
    all((start(x) / binSize) %% 1 == 0) &
    all((end(x) / binSize) %% 1 == 0)
}

#' Check that a GInteractions object has been
#' snapped to bins
#'
#' @param x GInteractions object
#' @param binSize Integer (numeric) describing
#'  the new size of each range.
#'
#' @return Logical
#'
#' @importFrom S4Vectors first second
#'
#' @noRd
.checkSnappedPairs <- function(x, binSize) {
    r1 <- .checkSnappedRanges(x = first(x), binSize = binSize)
    r2 <- .checkSnappedRanges(x = second(x), binSize = binSize)
    return(r1 & r2)
}

#' Check that a GRanges object has been binned
#'
#' Starts are 0-based for interfacing with the
#' `strawr` package. Therefore, widths of correctly
#' binned objects will be `binSize+1`.
#'
#' @param x GRanges object
#' @param binSize Integer (numeric) describing
#'  the new size of each range.
#'
#' @return Logical
#'
#' @importFrom GenomicRanges width
#'
#' @noRd
.checkBinnedRanges <- function(x, binSize) {
    length(which(width(x) != binSize+1)) == 0
}

#' Check that a GInteractions object has been binned
#'
#' Starts are 0-based for interfacing with the
#' `strawr` package. Therefore, widths of correctly
#' binned objects will be `binSize+1`.
#'
#' @param x GInteractions object
#' @param binSize Integer (numeric) describing
#'  the new size of each range.
#'
#' @return Logical
#'
#' @importFrom S4Vectors first second
#'
#' @noRd
.checkBinnedPairs <- function(x, binSize) {

    r1 <- .checkBinnedRanges(x = first(x), binSize = binSize)
    r2 <- .checkBinnedRanges(x = second(x), binSize = binSize)

    return(r1 & r2)
}

#' Return the mode(s)
#' @param x numeric vector
#' @returns A vector of the mode(s)
#' @noRd
.modes <- function(x) {
    ux <- unique(x)
    tab <- tabulate(match(x, ux))
    ux[tab == max(tab)]
}

