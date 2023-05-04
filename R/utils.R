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

## TODO: Write method for combining
## HDF5 data into a single file in blocks

#' Check that a list of objects contains
#' the same data in a slot.
#' @param x List of objects.
#' @param FUN Slot accessor function.
#' @returns Logical that all objects contain the same
#'  data or not.
#' @noRd
.checkEqualSlot <- function(x, FUN) {
    d <- lapply(x, FUN)
    all(vapply(seq_along(d), \(i) identical(d[[1L]], d[[i]]), logical(1L)))
}

#' Internal rbind method for InteractionMatrix/InteractionArray
#' @param ... InteractionMatrix objects to be combined row-wise.
#'  All objects must be the same class.
#' @param deparse.level An integer scalar; see `?base::rbind` for
#'  a description of this argument.
#' @importFrom S4Vectors metadata `metadata<-`
#' @importFrom SummarizedExperiment colData `colData<-`
#' @importFrom rlang abort
#' @importFrom glue glue
#' @noRd
.rbindIsetDerived <- function(..., deparse.level=1) {
    args <- unname(list(...))
    type <- class(args[[1L]]) # get class name

    ## Check equivalent metadata before binding
    if (!.checkEqualSlot(args, metadata)) {
        abort(glue("Can't rbind {type} \\
                    objects with different metadata."))
    }

    ## Check equivalent colData before binding
    if (!.checkEqualSlot(args, colData)) {
        abort(glue("Can't rbind {type} \\
                    objects with different colData."))
    }

    ans <- new(type, callNextMethod())
    metadata(ans) <- metadata(args[[1L]])
    colData(ans) <- colData(args[[1L]])
    ans
}

#' Internal cbind method for InteractionMatrix/InteractionArray
#' @param ... InteractionMatrix objects to be combined column-wise.
#'  All objects must be the same class.
#' @param deparse.level An integer scalar; see `?base::cbind` for
#'  a description of this argument.
#' @importFrom S4Vectors metadata `metadata<-`
#' @importFrom rlang abort
#' @importFrom glue glue
#' @noRd
.cbindIsetDerived <- function(..., deparse.level=1) {
    args <- unname(list(...))
    type <- class(args[[1L]]) # get class name

    ## Check equivalent metadata before binding
    if (!.checkEqualSlot(args, metadata)) {
        abort(glue("Can't cbind {type} \\
                    objects with different metadata."))
    }

    tryCatch({
        ans <- new(type, callNextMethod())
    }, error=\(e) {
        abort(e$message, call=parent.frame(4L))
    })

    metadata(ans) <- metadata(args[[1L]])
    ans
}

#' Stop if buffer is not the same
#' @param b1 buffer (numeric) from first object
#' @param b2 buffer (numeric) from second object
#' @importFrom rlang abort
#' @return NULL or error message if not the same.
#' @noRd
.checkBuffer <- function(b1, b2) {
    if (b1 != b2) {
        abort("`buffer` must be the same for both selections.")
    }
}

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

#' Stop if matrices are not odd and square
#' @param x InteractionArray
#' @importFrom rlang abort
#' @importFrom glue glue
#' @return NULL or error message if not odd and square.
#' @noRd
.checkOddSquareMatrices <- function(x){
  dims <- dim(counts(x))
  
  ## Check that input is a square matrix
  if(dims[1] != dims[2]){
    abort(c("`x` must have square count matrices.",
            "i"="Dimensions of count matrices must be equal.",
            "x"=glue("`dim(counts(x))[1] != dim(counts(x))[2]`",
                     ", {dims[1]} != {dims[2]}."),
            "i"="See `?pullHicMatrices` for more information."))
  }
  
  ## Check that buffer for InteractionArray is odd
  if((dims[1] %% 2) == 0){
    abort(c(glue("Enrichment scores can only be calculated for matrices",
                 " with a center pixel."),
            "i"="Dimensions of count matrices must be odd.",
            "i"=glue("Dimensions of count matrices are {dims[1]} x {dims[2]}."),
            "x"= glue("{dims[1]} is not odd."),
            "i"="See `?pixelsToMatrices` for help."))
  }
}

#' Return default buffer
#' If InteractionArray is supplied,
#' it uses the counts matrices to
#' set the buffer dimensions.
#' @param x InteractionArray
#' @return 5 (set default), 
#'  the buffer of the provided InteractionArray,
#'  or an error message if the InteractionArray 
#'  is not odd and square (no buffer)
#' @rdname defaultBuffer
#' @export
defaultBuffer <- function(x) {
  if (missing(x)) {
    return(5)
  }
  .checkOddSquareMatrices(x)
  buffer <- (dim(counts(x))[1] - 1) /2
  buffer
}


