## Generate manhattan, diagonal & angle matrices -------------------------------

#' Make a manhattan distance matrix
#' around a central pixel
#'
#' The center pixel is denoted as "0".
#'
#' @param buffer Integer indicating the
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
#' @param buffer Integer indicating the
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
                "i"="Set binSize with `assignToBins(x, binSize)`."))
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
#' @param x Integer indicating the
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

## Accessor methods ------------------------------------------------------------

#' Names S3 method for autocomplete
#' @param x A `MatrixSelection` object.
#' @rdname MatrixSelection
#' @return different data or metadata concerning a MatrixSelection object
#' @keywords internal
#' @exportS3Method
names.MatrixSelection <- function(x) slotNames(x)

#' Extract `$` operator for MatrixSelection
#' @param name Name of slot.
#' @rdname MatrixSelection
#' @keywords internal
#' @export
setMethod("$", "MatrixSelection", function(x, name) slot(x, name))

#' Extract `[[` operator for MatrixSelection
#' @param x A `MatrixSelection` object.
#' @param i A character or numeric to extract.
#' @rdname MatrixSelection
#' @keywords internal
#' @export
setMethod("[[", "MatrixSelection", function(x, i) {
    ifelse(is.character(i), slot(x,i), slot(x, slotNames(x)[i]))
})

#' Visualize selection for a MatrixSelection object
#'
#' Note: that buffer must be the same
#' as the selection functions to work appropriately
#'
#' @param object A MatrixSelection object.
#' @importFrom rlang inform
#' @importFrom utils capture.output
#' @noRd
.showSelection <- function(object) {
    x <- object@x
    b <- object@buffer
    m <- matrix(data="- ", nrow=b*2+1, ncol=b*2+1)
    m[x] <- "0 "
    prm <- capture.output(
        prmatrix(x=m,
                 rowlab=rep("",nrow(m)),
                 collab=rep("",ncol(m)),
                 quote=FALSE)
    )
    leg <- c("'0' = selected; '-' = unselected", prm)
    for(i in leg) {
        inform(i, appendLF=TRUE)
    }
}


#' Visualize foreground/background selection
#' @param fg Integer vector of foreground
#'  selected matrix indices.
#' @param fg Integer vector of background
#'  selected matrix indices.
#' @param buffer Integer describing the number of pixels
#'  surrounding the central pixel.
#' @noRd
.showMultiSelection <- function(fg, bg, buffer) {
    m <- matrix(data="- ", nrow=buffer*2+1, ncol=buffer*2+1)
    m[fg] <- "0 "
    m[bg] <- "X "
    m[intersect(fg, bg)] <- "* "
    prm <- capture.output(
        prmatrix(x=m,
                 rowlab=rep("",nrow(m)),
                 collab=rep("",ncol(m)),
                 quote=FALSE)
    )
    leg <- c("'0' = foreground;",
             "'X' = background;",
             "'*' = both;",
             "'-' = unselected", prm)
    for(i in leg) {
        inform(i, appendLF=TRUE)
    }
}

#' Visualize selection for a MatrixSelection object
#'
#' Note: that buffer must be the same
#' as the selection functions to work appropriately
#'
#' @param object A MatrixSelection object.
#' @returns A text-based visualization of the select matrix
#'  indices.
#' @examples
#' res <- selectCenterPixel(0, 3)
#' show(res)
#' @rdname selection-functions
#' @export
setMethod("show", "MatrixSelection", .showSelection)

#' Concatenate MatrixSelection objects
#' @rdname MatrixSelection
#' @export
setMethod("+",
          signature(e1="MatrixSelection",
                    e2="MatrixSelection"),
          function(e1, e2) {
              .checkBuffer(e1@buffer, e2@buffer)
              MatrixSelection(x=c(e1@x, e2@x),
                              buffer=e1@buffer)
          })

#' Concatenate MatrixSelection objects
#' @rdname MatrixSelection
#' @export
setMethod("-",
          signature(e1="MatrixSelection",
                    e2="MatrixSelection"),
          function(e1, e2) {
              .checkBuffer(e1@buffer, e2@buffer)
              MatrixSelection(x=setdiff(e1@x, e2@x),
                              buffer=e1@buffer)
          })


## Selection functions ---------------------------------------------------------

#' Return indices matching provided
#' manhattan distances
#' @inheritParams selectRadius
#' @noRd
.selectRadius <- function(x, buffer, invert) {
    mh <- .manhattanMatrix(buffer)
    if (invert) {
        return(which(!mh %in% x))
    } else {
        return(which(mh %in% x))
    }
}

#' Return indices matching provided
#' manhattan distances
#' @param x Integer vector of manhattan distances
#'  to select.
#' @param buffer Integer describing the number of pixels
#'  surrounding the central pixel.
#' @param invert Boolean indicating whether to invert
#'  the selection.
#' @return Numeric vector of matrix indices (byRow).
#' @rdname selection-functions
#' @examples
#' selectRadius(x=c(2,3,4), buffer=5, invert=FALSE)
#' @export
setMethod("selectRadius", "numeric", function(x, buffer, invert) {
    ans <- .selectRadius(x, buffer, invert)
    MatrixSelection(x=ans, buffer=buffer)
})

#' Select center with a manhattan distance radius
#' @inheritParams selectCenterPixel
#' @noRd
.selectCenterPixel <- function(mhDist, buffer, invert) {
    .selectRadius(c(0, mhDist), buffer, invert)
}

#' Select center with a manhattan distance radius
#' @param mhDist Integer vector of manhattan distances
#'  to select along with center pixel.
#' @param buffer Integer describing the number of pixels
#'  surrounding the central pixel.
#' @param invert Boolean indicating whether to invert
#'  the selection.
#' @return Numeric vector of matrix indices (byRow).
#' @rdname selection-functions
#' @examples
#' selectCenterPixel(0, 5)
#' @export
setMethod("selectCenterPixel", "numeric", function(mhDist, buffer, invert) {
    ans <- .selectCenterPixel(mhDist, buffer, invert)
    MatrixSelection(x=ans, buffer=buffer)
})


#' Accepts matrix and returns indices
#' with corresponding 1 values.
#' @inheritParams selectSubmatrix
#' @noRd
.selectSubmatrix <- function(m, invert) {
    if (invert) which(m != 1) else which(m == 1)
}

#' Select center with a manhattan distance radius
#' TODO: Enforce dimensions
#' @param m matrix with 1's indicating selected positions and
#'  0's indicated unselected positions.
#' @param invert Boolean indicating whether to invert
#'  the selection.
#' @return Numeric vector of matrix indices (byRow).
#' @rdname selection-functions
#' @examples
#' selectSubmatrix(m = matrix(rep(c(1,0,1), 3), nrow=3, ncol=3))
#' @export
setMethod("selectSubmatrix", "matrix", function(m, invert) {
    ans <- .selectSubmatrix(m, invert)
    dims <- dim(m)
    m <- matrix(data="- ", nrow=dims[1], ncol=dims[2])
    m[ans] <- "0 "
    prm <- capture.output(prmatrix(x=m,
                                   rowlab=rep("",nrow(m)),
                                   collab=rep("",ncol(m)),
                                   quote=FALSE))
    for(i in c("'0' = selected; '- ' = unselected", prm)) {
        inform(i, appendLF=TRUE)
    }
    ans
})

#' Select coordinates from matrix
#' @inheritParams selectCoordinates
#' @noRd
.selectCoordinates <- function(rowInd, colInd, buffer, invert) {
    w <- buffer*2+1
    m <- matrix(seq(1,w^2), nrow=w, ncol=w)
    ind <- m[cbind(rowInd, colInd)]
    if (invert) which(!m %in% ind) else ind
}

#' Select coordinates from matrix
#'
#' For `selectCoordinates`, `rowInd` and `colInd` are paired
#' such that the selected position in the matrix is
#' `c(rowInd[1:i], colInd[1:j])` for `i` rows and `j` columns.
#'
#' @param rowInd Integer describing the row index coordinate.
#' @param colInd Integer describing the column index coordinate.
#' @param buffer Integer describing the number of pixels
#'  surrounding the central pixel.
#' @param invert Boolean indicating whether to invert
#'  the selection.
#' @return Numeric vector of matrix indices (byRow).
#' @examples
#' selectCoordinates(rowInd=1:3, colInd=1:3, buffer=5)
#' @rdname selection-functions
#' @export
setMethod("selectCoordinates", "numeric",
          function(rowInd, colInd, buffer, invert) {
              ans <- .selectCoordinates(rowInd, colInd, buffer, invert)
              MatrixSelection(x=ans, buffer=buffer)
          })

#' Select continuous block of matrix
#' @inheritParams selectBlock
#' @noRd
.selectBlock <- function(rowInd, colInd, buffer, invert) {
    w <- buffer*2+1
    m <- matrix(seq(1,w^2), nrow=w, ncol=w)
    ind <- as.vector(m[rowInd, colInd])
    if (invert) which(!m %in% ind) else ind
}

#' Select continuous block of matrix
#' @param rowInd Integer describing the row indices.
#' @param colInd Integer describing the column indices.
#' @param buffer Integer describing the number of pixels
#'  surrounding the central pixel.
#' @param invert Boolean indicating whether to invert
#'  the selection.
#' @return Numeric vector of matrix indices (byRow).
#' @examples
#' selectBlock(rowInd=1:3, colInd=1:3, buffer=5)
#' @rdname selection-functions
#' @export
setMethod("selectBlock", "numeric",
          function(rowInd, colInd, buffer, invert) {
              ans <- .selectBlock(rowInd, colInd, buffer, invert)
              MatrixSelection(x=ans, buffer=buffer)
          })

#' Select top left corner of a matrix
#' @inheritParams selectTopLeft
#' @noRd
.selectTopLeft <- function(n, buffer, inset, invert) {
    w <- buffer*2+1
    rowInd <- seq(1, n) + inset
    colInd <- seq(1, n) + inset
    .selectBlock(rowInd, colInd, buffer, invert)
}

#' Select top left corner of a matrix
#' @param n Integer describing the number of square pixels
#'  to select.
#' @param buffer Integer describing the number of pixels
#'  surrounding the central pixel.
#' @param inset Integer describing the number of pixels to
#'  inset the selection from the outer edge of the matrix.
#'  Default of 0 uses no inset.
#' @param invert Boolean indicating whether to invert
#'  the selection.
#' @return Numeric vector of matrix indices (byRow).
#' @examples
#' selectTopLeft(n=3, buffer=5, inset=1, invert=FALSE)
#' @rdname selection-functions
#' @export
setMethod("selectTopLeft", "numeric",
          function(n, buffer, inset, invert) {
              ans <- .selectTopLeft(n, buffer, inset, invert)
              MatrixSelection(x=ans, buffer=buffer)
          })



#' Top right corner
#' @inheritParams selectTopRight
#' @noRd
.selectTopRight <- function(n, buffer, inset, invert) {
    w <- buffer*2+1
    rowInd <- seq(1, n) + inset
    colInd <- seq(w-n+1, w) - inset
    .selectBlock(rowInd, colInd, buffer, invert)
}

#' Select top right corner of a matrix
#' @param n Integer describing the number of square pixels
#'  to select.
#' @param buffer Integer describing the number of pixels
#'  surrounding the central pixel.
#' @param inset Integer describing the number of pixels to
#'  inset the selection from the outer edge of the matrix.
#'  Default of 0 uses no inset.
#' @param invert Boolean indicating whether to invert
#'  the selection.
#' @return Numeric vector of matrix indices (byRow).
#' @examples
#' selectTopRight(n=3, buffer=5, inset=1, invert=FALSE)
#' @rdname selection-functions
#' @export
setMethod("selectTopRight", "numeric",
          function(n, buffer, inset, invert) {
              ans <- .selectTopRight(n, buffer, inset, invert)
              MatrixSelection(x=ans, buffer=buffer)
          })


#' Select bottom right corner of a matrix
#' @inheritParams selectBottomRight
#' @noRd
.selectBottomRight <- function(n, buffer, inset, invert) {
    w <- buffer*2+1
    rowInd <- seq(w-n+1, w) - inset
    colInd <- seq(w-n+1, w) - inset
    .selectBlock(rowInd, colInd, buffer, invert)
}

#' Select bottom right corner of a matrix
#' @param n Integer describing the number of square pixels
#'  to select.
#' @param buffer Integer describing the number of pixels
#'  surrounding the central pixel.
#' @param inset Integer describing the number of pixels to
#'  inset the selection from the outer edge of the matrix.
#'  Default of 0 uses no inset.
#' @param invert Boolean indicating whether to invert
#'  the selection.
#' @return Numeric vector of matrix indices (byRow).
#' @examples
#' selectBottomRight(n=3, buffer=5, inset=1, invert=FALSE)
#' @rdname selection-functions
#' @export
setMethod("selectBottomRight", "numeric",
          function(n, buffer, inset, invert) {
              ans <- .selectBottomRight(n, buffer, inset, invert)
              MatrixSelection(x=ans, buffer=buffer)
          })


#' Select bottom left corner of a matrix
#' @inheritParams selectBottomLeft
#' @noRd
.selectBottomLeft <- function(n, buffer, inset, invert) {
    w <- buffer*2+1
    rowInd <- seq(w-n+1, w) - inset
    colInd <- seq(1, n) + inset
    .selectBlock(rowInd, colInd, buffer, invert)
}

#' Select bottom left corner of a matrix
#' @param n Integer describing the number of square pixels
#'  to select.
#' @param buffer Integer describing the number of pixels
#'  surrounding the central pixel.
#' @param inset Integer describing the number of pixels to
#'  inset the selection from the outer edge of the matrix.
#'  Default of 0 uses no inset.
#' @param invert Boolean indicating whether to invert
#'  the selection.
#' @return Numeric vector of matrix indices (byRow).
#' @examples
#' selectBottomLeft(n=3, buffer=5, inset=1, invert=FALSE)
#' @rdname selection-functions
#' @export
setMethod("selectBottomLeft", "numeric",
          function(n, buffer, inset, invert) {
              ans <- .selectBottomLeft(n, buffer, inset, invert)
              MatrixSelection(x=ans, buffer=buffer)
          })


#' Select corners
#' @inheritParams selectCorners
#' @noRd
.selectCorners <- function(n, buffer, inset, invert) {
    ind <- c(
        .selectTopLeft(n, buffer, inset, invert=FALSE),
        .selectTopRight(n, buffer, inset, invert=FALSE),
        .selectBottomRight(n, buffer, inset, invert=FALSE),
        .selectBottomLeft(n, buffer, inset, invert=FALSE)
    )
    if (invert) {
        w <- buffer*2+1
        m <- matrix(seq(1,w^2), nrow=w, ncol=w)
        which(!m %in% ind)
    } else {
        ind
    }
}

#' Select corners of a matrix
#' @param n Integer describing the number of square pixels
#'  to select.
#' @param buffer Integer describing the number of pixels
#'  surrounding the central pixel.
#' @param inset Integer describing the number of pixels to
#'  inset the selection from the outer edge of the matrix.
#'  Default of 0 uses no inset.
#' @param invert Boolean indicating whether to invert
#'  the selection.
#' @return Numeric vector of matrix indices (byRow).
#' @examples
#' selectCorners(n=3, buffer=5, inset=1, invert=FALSE)
#' @rdname selection-functions
#' @export
setMethod("selectCorners", "numeric",
          function(n, buffer, inset, invert) {
              ans <- .selectCorners(n, buffer, inset, invert)
              MatrixSelection(x=ans, buffer=buffer)
          })


#' Select entire rows of a matrix
#' @inheritParams selectRows
#' @noRd
.selectRows <- function(rows, buffer, invert) {
    w <- buffer*2+1
    m <- matrix(seq(1,w^2), nrow=w, ncol=w)
    ind <- as.vector(m[rows,])
    if (invert) which(!m %in% ind) else ind
}

#' Select entire rows of a matrix
#' @param rows Integer describing which rows to select.
#' @param buffer Integer describing the number of pixels
#'  surrounding the central pixel.
#' @param invert Boolean indicating whether to invert
#'  the selection.
#' @return Numeric vector of matrix indices (byRow).
#' @examples
#' selectRows(rows=1:3, buffer=5, invert=FALSE)
#' @rdname selection-functions
#' @export
setMethod("selectRows", "numeric",
          function(rows, buffer, invert) {
              ans <- .selectRows(rows, buffer, invert)
              MatrixSelection(x=ans, buffer=buffer)
          })

#' Select entire columns of matrix
#' @inheritParams selectCols
#' @noRd
.selectCols <- function(cols, buffer, invert) {
    w <- buffer*2+1
    m <- matrix(seq(1,w^2), nrow=w, ncol=w)
    ind <- as.vector(m[,cols])
    if (invert) which(!m %in% ind) else ind
}

#' Select entire columns of matrix
#' @param cols Integer describing which cols to select.
#' @param buffer Integer describing the number of pixels
#'  surrounding the central pixel.
#' @param invert Boolean indicating whether to invert
#'  the selection.
#' @return Numeric vector of matrix indices (byRow).
#' @examples
#' selectCols(cols=1:3, buffer=5, invert=FALSE)
#' @rdname selection-functions
#' @export
setMethod("selectCols", "numeric",
          function(cols, buffer, invert) {
              ans <- .selectCols(cols, buffer, invert)
              MatrixSelection(x=ans, buffer=buffer)
          })

#' Select inner square of matrix
#' @inheritParams selectInner
#' @noRd
.selectInner <- function(n, buffer, invert) {
    cp <- buffer+1
    ind <- seq(cp-n, cp+n)
    .selectBlock(ind, ind, buffer, invert)
}

#' Select inner square of matrix
#' @param n Integer describing the number of square pixels
#'  to select. Must be length of one.
#' @param buffer Integer describing the number of pixels
#'  surrounding the central pixel.
#' @param invert Boolean indicating whether to invert
#'  the selection.
#' @return Numeric vector of matrix indices (byRow).
#' @examples
#' selectInner(n=1, buffer=5, invert=FALSE)
#' @rdname selection-functions
#' @export
setMethod("selectInner", "numeric",
          function(n, buffer, invert) {
              ans <- .selectInner(n, buffer, invert)
              MatrixSelection(x=ans, buffer=buffer)
          })

#' Select outer square of matrix
#' @inheritParams selectOuter
#' @noRd
.selectOuter <- function(n, buffer, invert) {
    w <- buffer*2+1
    m <- matrix(seq(1,w^2), nrow=w, ncol=w)
    tRim <- seq(1, n)
    bRim <- seq(w, w-n+1)
    tr <- .selectRows(tRim, buffer, invert=FALSE)
    br <- .selectRows(bRim, buffer, invert=FALSE)
    tc <- .selectCols(tRim, buffer, invert=FALSE)
    bc <- .selectCols(bRim, buffer, invert=FALSE)

    ind <- unique(c(tr, br, tc, bc))
    if (invert) which(!m %in% ind) else ind
}

#' Select outer square of matrix
#' @param n Integer describing the number of outer pixels
#'  to select. Must be length of one.
#' @param buffer Integer describing the number of pixels
#'  surrounding the central pixel.
#' @param invert Boolean indicating whether to invert
#'  the selection.
#' @return Numeric vector of matrix indices (byRow).
#' @examples
#' selectOuter(n=1, buffer=5, invert=FALSE)
#' @rdname selection-functions
#' @export
setMethod("selectOuter", "numeric",
          function(n, buffer, invert) {
              ans <- .selectOuter(n, buffer, invert)
              MatrixSelection(x=ans, buffer=buffer)
          })

