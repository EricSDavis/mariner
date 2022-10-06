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

#' Check input types
#'
#' Derived types:
#' string - length one character vector
#' number - length one numeric vector
#' boolean - a length one logical vector that is not NA
#'
#' @param types Named vector or list of arguments and their types
#' @importFrom rlang abort
#' @importFrom glue glue
#' @importFrom assertthat is.string is.number is.flag
#' @returns NULL or an error message incorrect types.
#' @noRd
.checkTypes <- function(types, env=parent.frame()) {
    args <- names(types)
    for(i in seq_along(types)) {
        if (types[i] == "string") {
            if(any(!is.string(get(args[i], envir=env)))) {
                abort(glue(
                    "{args[i]} is not a string \\
                        (a length one character vector)."
                ))
            }
        }
        if (types[i] == "number") {
            if(any(!is.number(get(args[i], envir=env)))) {
                abort(glue(
                    "{args[i]} is not a number \\
                        (a length one numeric vector)."
                ))
            }
        }
        if (types[i] == "boolean") {
            arg <- get(args[i], envir=env)
            if(any(!is.flag(arg) | is.na(arg))) {
                abort(glue(
                    "{args[i]} is not a boolean \\
                        (a length one logical vector that is not NA)."
                ))
            }
        }
    }
}
