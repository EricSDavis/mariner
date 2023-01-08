#' Get binSize or throw error
#' @param x GInteractions object.
#' @importFrom S4Vectors first second
#' @importFrom IRanges width
#' @importFrom rlang abort
#' @noRd
.getBinSize <- function(x) {
    widths <- unique(width(regions(x))) - 1
    if (length(widths) != 1L) {
        abort(c("All ranges in `x` must be equal widths.",
                "i"="Use `binPairs()` to bin into equal widths."))
    }
    return(widths)
}

#' Make a manhattan distance matrix
#' around a central pixel
#'
#' The center pixel is denoted as "0".
#'
#' @param buffer Integer indicatinig the
#'  buffer size, or number of pixels
#'  around a matrix.
#' @importFrom rlang abort
#' @returns A matrix where each position
#'  indicates the manhattan distance from
#'  the center pixel.
#' @noRd
.manhattanMatrix <- function(buffer) {
    if (buffer < 0) abort("`buffer` must be >= 0.")
    if (buffer == 0) return(matrix(0))
    b <- buffer
    a <- b*2
    FUN <- \(a, b) c(seq(a, a-b), seq(a-b+1, a))
    vapply(FUN(a, b), FUN, b=b, integer(a+1))
}

#' Make a matrix indicating the distance
#' from the diagonal
#' @param x A GInteractions object.
#' @param buffer Integer indicatinig the
#'  buffer size, or number of pixels
#'  around a matrix.
#' @importFrom InteractionSet pairdist
#' @noRd
.diagonalMatrix <- function(x, buffer) {

    ## Get binSize
    binSize <- .getBinSize(x)
    if (length(binSize) != 1L) {
        abort(c(glue("All interactions in `x` must be \\
                     the same width."),
                "i"="Check this with `width(x)`.",
                "i"="Set binSize with `binPairs(x, binSize)`."))
    }
    ## How does moving down or right change dist to diag?
    p <- pairdist(x)

    ## TODO: swap these if starts & ends are swapped
    ## Row function
    rows <- \(x) as.integer(seq(x, x-(binSize*(buffer*2)), -binSize))

    ## Apply row function for each column
    FUN <- \(x) {
        cols <- as.integer(seq(x, x+(binSize*(buffer*2)), binSize))
        mat <- vapply(cols, rows, integer(buffer*2+1))
        mat[mat < 0] <- NA
        mat
    }

    ## Apply for each loop
    template <- array(0,c(buffer*2+1,buffer*2+1,1))
    ans <- drop(vapply(p, FUN, template))
    ans

}

#' Create angle matrix
#'
#' @param x Integer indicatinig the
#'  buffer size, or number of pixels
#'  around a matrix.
#' @noRd
.angleMatrix <- function(x) {
    ## Define quadrants
    rotate <- \(x) t(apply(x,2,rev))
    height <- matrix(rep(seq(x,0), x+1), nrow=x+1, ncol=x+1)
    base <- rotate(height)
    q4 <- rotate(atan(height/base)*(180/pi) + 270)
    q3 <- rotate(q4 - 90)
    q2 <- rotate(q3 - 90)
    q1 <- rotate(q2 - 90)

    ## Assign quadrants
    M <- matrix(NA, nrow=x*2+1, ncol=x*2+1)
    upper <- seq(1, x+1)
    lower <- seq(x+1, x*2+1)
    M[lower, lower] <- q4
    M[upper, upper] <- q2
    M[lower, upper] <- q3
    M[upper, lower] <- q1
    M
}


#' Return indices matching provided
#' manhattan distances
#'
#' @param x Integer vector of manhattan distances
#'  to select.
#' @param buffer Integer describing the number of pixels
#'  surrounding the central pixel.
#' @param invert Boolean indicating whether to invert
#'  the selection.
#' @noRd
.selectRadius <- function(x, buffer, invert) {
    mh <- .manhattanMatrix(buffer)
    if (invert) {
        return(which(!mh %in% x))
    } else {
        return(which(mh %in% x))
    }
}

#' Select center with a manhattan distance radius
#'
#' @param mhDist Integer vector of manhattan distances
#'  to select along with center pixel.
#' @param buffer Integer describing the number of pixels
#'  surrounding the central pixel.
#' @param invert Boolean indicating whether to invert
#'  the selection.
#'
#' @noRd
.selectCenterPixel <- function(mhDist, buffer, invert) {
    .selectRadius(c(0, mhDist), buffer, invert)
}

#' TODO: Write function for selecting angles
#' @noRd
.selectAngles <- function() {

}

#' TODO: Write function for setting diagonal
#' @noRd
.selectDiagonal <- function() {

}

#' TODO: Select submatrix
#' Accepts matrix and returns indices
#' with corresponding 1 values.
#' @noRd
.selectSubmatrix <- function(m, invert) {
    if (invert) which(m != 1) else which(m == 1)
}

#' Select coordinates from matrix
#' @noRd
.selectCoordinates <- function(buffer, rowInd, colInd, invert) {
    w <- buffer*2+1
    m <- matrix(1:w^2, nrow=w, ncol=w)
    ind <- m[cbind(rowInd, colInd)]
    if (invert) which(!m %in% ind) else ind
}

#' Select continuous block of matrix
#' @noRd
.selectBlock <- function(buffer, rowInd, colInd, invert) {
    w <- buffer*2+1
    m <- matrix(1:w^2, nrow=w, ncol=w)
    ind <- as.vector(m[rowInd, colInd])
    if (invert) which(!m %in% ind) else ind
}

#' Shortcuts for selecting corners of a matrix
#' TODO: add inset parameter
#' @noRd
.selectTopLeft <- function(buffer, n, invert) {
    w <- buffer*2+1
    rowInd <- seq(1, n)
    colInd <- seq(1, n)
    .selectBlock(buffer, rowInd, colInd, invert)
}

#' Top right corner
#' @noRd
.selectTopRight <- function(buffer, n, invert) {
    w <- buffer*2+1
    rowInd <- seq(1, n)
    colInd <- seq(w-n+1, w)
    .selectBlock(buffer, rowInd, colInd, invert)
}

#' Bottom right corner
#' @noRd
.selectBottomRight <- function(buffer, n, invert) {
    w <- buffer*2+1
    rowInd <- seq(w-n+1, w)
    colInd <- seq(w-n+1, w)
    .selectBlock(buffer, rowInd, colInd, invert)
}

#' Bottom left corner
#' @noRd
.selectBottomLeft <- function(buffer, n, invert) {
    w <- buffer*2+1
    rowInd <- seq(w-n+1, w)
    colInd <- seq(1, n)
    .selectBlock(buffer, rowInd, colInd, invert)
}

#' Select corners
#' @noRd
.selectCorners <- function(buffer, n, invert) {
    ind <- c(
        .selectTopLeft(buffer, n, invert=FALSE),
        .selectTopRight(buffer, n, invert=FALSE),
        .selectBottomRight(buffer, n, invert=FALSE),
        .selectBottomLeft(buffer, n, invert=FALSE)
    )
    if (invert) {
        w <- buffer*2+1
        m <- matrix(1:w^2, nrow=w, ncol=w)
        which(!m %in% ind)
    } else {
        ind
    }
}

#' Select entire rows of matrix
#' @noRd
.selectRows <- function(buffer, rows, invert) {
    w <- buffer*2+1
    m <- matrix(1:w^2, nrow=w, ncol=w)
    ind <- as.vector(m[rows,])
    if (invert) which(!m %in% ind) else ind
}

#' Select entire columns of matrix
#' @noRd
.selectCols <- function(buffer, cols, invert) {
    w <- buffer*2+1
    m <- matrix(1:w^2, nrow=w, ncol=w)
    ind <- as.vector(m[,cols])
    if (invert) which(!m %in% ind) else ind
}

#' TODO: Select offset diagonal
#' Include param for direction 1 (TLtoBR) or 2 (BLtoTR)
#' @noRd
.selectDiagonal <- function() {


}

#' Select inner square of matrix
#' @noRd
.selectInner <- function(buffer, n, invert) {
    cp <- buffer+1
    ind <- seq(cp-n, cp+n)
    .selectBlock(buffer, ind, ind, invert)
}

#' TODO
#' Figure out inset
#' @noRd
.selectOuter <- function(buffer, n, invert) {
    w <- buffer*2+1
    m <- matrix(1:w^2, nrow=w, ncol=w)
    tRim <- seq(1, n)
    bRim <- seq(w, w-n+1)
    tr <- .selectRows(buffer, tRim, invert=FALSE)
    br <- .selectRows(buffer, bRim, invert=FALSE)
    tc <- .selectCols(buffer, tRim, invert=FALSE)
    bc <- .selectCols(buffer, bRim, invert=FALSE)

    ind <- unique(c(tr, br, tc, bc))
    if (invert) which(!m %in% ind) else ind
}


#' Visualize selection for a matrix
#'
#' Note: that buffer must be the same
#' as the selection functions to work appropriately
#'
#' @param x Integer vector of matrix indices.
#' @param buffer Integer describing the number of pixels
#'  surrounding the central pixel.
#' @importFrom rlang inform
#' @importFrom utils capture.output
#' @noRd
.showSelection <- function(x, buffer) {
    m <- matrix(data="- ", nrow=buffer*2+1, ncol=buffer*2+1)
    m[x] <- "X "
    prm <- capture.output(prmatrix(x=m,
                                   rowlab=rep("",nrow(m)),
                                   collab=rep("",ncol(m)),
                                   quote=FALSE))
    for(i in c("'X' = selected; '- ' = unselected", prm)) {
        inform(i, appendLF=TRUE)
    }
}

#' Select a subset of a loop matrix
#' Uses three filters to select a subset of
#' interactions. `mhDist` selects all interactions
#' belonging to provided manhattan distance vector.
#' `diagDist` selects all interactions greater than
#' the provided input
#' @noRd
.select <- function(mhDist, diagDist, angles, invert) {

}


#' Internal for calcLoopEnrichment
#' @inheritParams calcLoopEnrichment
#' @importFrom rlang abort inform
#' @importFrom glue glue
#' @importFrom DelayedArray RegularArrayGrid DelayedArray blockApply
#' @importFrom utils capture.output
#' @importFrom stats median
#' @noRd
.calcLoopEnrichment <- function(x, files, mhDist, nBlocks, verbose,
                                BPPARAM, ...) {

    ## Parameter parsing
    if (nBlocks <= 0) abort("`nBlocks` must be > 0.")

    ## Determine resolution from x
    ## and ensure all pixels are same res.
    binSize <- .getBinSize(x)
    if (length(binSize) != 1L) {
        abort(c(glue("All interactions in `x` must be \\
                     the same width."),
                "i"="Check this with `width(x)`.",
                "i"="Set binSize with `binPairs(x, binSize)`."))
    }

    ## Use mhDist set buffer for matrix ranges
    buffer <- max(mhDist)
    mr <- pixelsToMatrices(x, buffer)

    ## Pull count matrices
    ## Use nBlocks to set blockSize if not provided?
    iarr <- pullHicMatrices(mr, binSize, files, ...)
    cnts <- counts(iarr)

    ## Subset the center pixel of matrix
    cp <- cnts[(buffer+1), (buffer+1),,]

    ## Subset the background
    ## Make manhattan distance matrix
    mh <- .manhattanMatrix(buffer)
    ind <- which(mh %in% mhDist)

    ## Show selection
    .showSelection(ind, buffer)

    ## Build array grid
    spacings <- dims <- dim(cnts)
    spacings[3] <- ceiling(spacings[3]/nBlocks)
    grid <- RegularArrayGrid(dims, spacings)

    blockApplyFUN <- \(x) apply(x, c(3,4), \(x) median(x[ind]+1))
    combineFUN <- \(x) do.call("rbind", args=x)

    ## Apply in blocks
    blocks <- blockApply(x=cnts,
                         FUN=blockApplyFUN,
                         grid=grid,
                         verbose=verbose,
                         BPPARAM=BPPARAM)

    bg <- DelayedArray(combineFUN(blocks))

    return((cp+1)/bg)
}

#' Calculate loop enrichment over background.
#'
#' Pulls Hi-C pixels and calculates the enrichment
#' over background.
#'
#' @param x GInteractions object.
#' @param files Character file paths to `.hic` files.
#' @param mhDist Which manhattan distances to select
#' @param nBlocks Number of blocks for block-processing
#'  arrays. Default is 5. Increase this for large
#'  datasets. To read and process all data at once, set
#'  this value to 1.
#' @param verbose Boolean (TRUE or FALSE) describing
#'  whether to report block-processing progress.
#' @param BPPARAM Parallelization params (passed to
#'  `BiocParallel::bplapply()`). Default is the result
#'  of `BiocParallel::bpparams()`. Parallel processing
#'  is not available when `by=interactions`.
#' @param ... Additional arguments passed to
#'  `pullHicMatrices`. See ?[`pullHicMatrices`].
#' @returns A DelayedMatrix of enrichment scores
#'  where rows are interactions (i.e. loops) and
#'  columns are Hi-C files.
#'
#' @rdname calcLoopEnrichment
#' @export
setMethod("calcLoopEnrichment",
          signature(x="GInteractions",
                    files="character"),
          definition=.calcLoopEnrichment)
