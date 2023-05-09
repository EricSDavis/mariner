## Show ------------------------------------------------------------------------

#' show for InteractionJaggedArray
#' @param object InteractionJaggedArray object.
#' @importFrom InteractionSet interactions
#' @importFrom SummarizedExperiment colData
#'
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
#' ## Create test interactions
#' gi <- read.table(text="
#'             1 51000000 51300000 1 51000000 51500000
#'             2 52000000 52300000 3 52000000 52500000
#'             1 150000000 150500000 1 150000000 150300000
#'             2 52000000 52300000 2 52000000 52800000") |>
#'     as_ginteractions()
#'
#' ## InteractionJaggedArray object
#' iarr <- pullHicMatrices(gi, hicFiles, 100e03, half="both")
#' iarr
#'
#' @rdname InteractionJaggedArray-class
#' @export
setMethod("show", "InteractionJaggedArray", function(object) {
    output <- c()

    ## class & dims
    output[1] <- sprintf("class: %s", class(object))
    dims <- dim(object)
    output[2] <- sprintf(
        "dim: %s interaction(s), %s file(s), variable count matrix(es)",
        dims[[1]], dims[[2]]
    )

    ## metadata
    expt <- names(metadata(object))
    output[3] <- sprintf(
        "metadata(%i): %s",
        length(expt), paste0(expt, collapse=", ")
    )

    ## colData
    output[4] <- sprintf(
        "colData: %s",
        paste0(rownames(colData(object)), collapse=", ")
    )
    cdNames <- names(colData(object))
    output[5] <- sprintf(
        "colData names(%d): %s",
        length(cdNames), paste0(cdNames, collapse=", ")
    )

    ## hdf5
    output[6] <- sprintf("HDF5: %s", path(object))

    cat(output, sep = "\n")
})

## Accessors -------------------------------------------------------------------

#' dim for InteractionJaggedArray
#' @param x InteractionJaggedArray object.
#' @importFrom InteractionSet interactions
#' @importFrom rhdf5 h5read
#' @returns `dim()` returns a list of the dimensions
#'  of the interactions, files, and count matrices.
#' @examples
#' ## Show dimensions
#' dim(iarr)
#'
#' @rdname InteractionJaggedArray-class
#' @export
setMethod("dim", "InteractionJaggedArray", function(x) {
    ## Initialize list
    ans <- vector("list", 4L)
    names(ans) <- c("interactions", "files", "rows", "cols")

    ## Rearrange dimensions of JaggedArray
    dims <- dim(x@counts)
    ans[[1]] <- dims[[3]]
    ans[[2]] <- dims[[4]]
    ans[[3]] <- dims[[1]]
    ans[[4]] <- dims[[2]]
    ans
})

#' interactions
#' @param x InteractionJaggedArray object.
#' @returns `interactions()` returns the interactions.
#' @examples
#' ## Access interactions
#' interactions(iarr)
#'
#' @rdname InteractionJaggedArray-class
#' @export
setMethod("interactions", "InteractionJaggedArray", function(x) {
    x@interactions
})

#' metadata
#' @param x InteractionJaggedArray object.
#' @returns `metadata()` returns the metadata.
#' @examples
#' ## Access metadata
#' metadata(iarr)
#'
#' @rdname InteractionJaggedArray-class
#' @export
setMethod("metadata", "InteractionJaggedArray", function(x) {
    x@metadata
})

#' colData
#' @param x InteractionJaggedArray object.
#' @returns `colData()` returns the column data.
#' @examples
#' ## Access colData
#' colData(iarr)
#'
#' @rdname InteractionJaggedArray-class
#' @export
setMethod("colData", "InteractionJaggedArray", function(x) {
    x@colData
})

#' counts
#' @param object InteractionJaggedArray object.
#' @returns `counts()` returns the JaggedArray object
#'  containing count matrix information.
#' @examples
#' ## Access count matrices
#' counts(iarr)
#'
#' @rdname InteractionJaggedArray-class
#' @export
setMethod("counts", "InteractionJaggedArray", function(object) {
    object@counts
})

#' Access HDF5 path for InteractionJaggedArray
#' @param object InteractionJaggedArray object.
#' @returns `path()` returns a character vector with
#'  the path to the HDF5 file with the JaggedArray data.
#' @examples
#' ## Access path to HDF5 data
#' path(iarr)
#'
#' @rdname InteractionJaggedArray-class
#' @export
setMethod("path", "InteractionJaggedArray", function(object) {
    object@counts@h5File
})

#' length (number of interactions) for
#' InteractionJaggedArray
#' @param x InteractionJaggedArray object.
#' @returns `length()` returns an integer with the
#'  number of interactions in an InteractionJaggedArray object.
#' @examples
#' ## length
#' length(iarr)
#'
#' @rdname InteractionJaggedArray-class
#' @export
setMethod("length", "InteractionJaggedArray", function(x) {
    length(x@interactions)
})

## Subsetting ------------------------------------------------------------------

#' Internal subsetting for InteractionJaggedArray
#' @noRd
.subsetInteractionJaggedArray <- function(x, i, j) {
    dims <- dim(x)

    ## Catch NULL indices
    if (is.null(i)) i <- seq_len(dims$interactions)
    if (is.null(j)) j <- seq_len(dims$files)

    ## Update interactions
    x@interactions <- x@interactions[i,]

    ## Update colData
    x@colData <- x@colData[j,]

    ## Update counts
    cnts <- x@counts[,,i,j]

    ## Convert to InteractionArray if possible
    if (is(cnts, "JaggedArray")) {
        x@counts <- x@counts[,,i,j]
    }
    if (is(cnts, "DelayedArray")) {
        iarr <- InteractionArray(
            interactions = x@interactions,
            assays = list(
                counts = aperm(cnts, c(3,4,1,2))
            ),
            colData = x@colData,
            metadata = x@metadata
        )
        return(iarr)
    }

    return(x)
}

#' Indexing for InteractionJaggedArray
#'
#' Subset an InteractionJaggedArray by its interactions
#' ([i,]) or its Hi-C files ([,j]).
#'
#' The object returned will be a InteractionJaggedArray
#' if the submatrices contain different dimensions.
#' However, the returned object will automatically
#' be coerced into a InteractionArray if possible (i.e.
#' the dimensions of the rows and columns of
#' submatrices are the same.)
#'
#' @param x An InteractionJaggedArray object.
#' @param i Numeric vector indicating the indices
#'  of interactions to extract.
#' @param j Numeric vector indicating the indices
#'  of files to extract.
#' @returns Subsetting returns an InteractionJaggedArray
#'  or InteractionArray object (see Details).
#' @examples
#' ## Subsetting
#' iarr[1:3,1]
#'
#' @rdname InteractionJaggedArray-class
#' @export
setMethod("[", "InteractionJaggedArray", function(x, i, j) {
    if (missing(i)) i <- NULL
    if (missing(j)) j <- NULL
    .subsetInteractionJaggedArray(x, i, j)
})

## Overlaps --------------------------------------------------------------------

#' Overlap methods for InteractionJaggedArray
#'
#' @param query,subject An InteractionJaggedArray, Vector,
#'  GInteractions or InteractionSet object, depending on
#'  the specified method. At least one of these must be a
#'  `subject` can be missing if query is an
#'  InteractionJaggedArray object.
#' @param ... see ?`findOverlaps` in InteractionSet package
#'  for more information
#'
#' @importMethodsFrom IRanges findOverlaps countOverlaps
#'  overlapsAny subsetByOverlaps
#' @importFrom S4Vectors countQueryHits
#' @importFrom IRanges overlapsAny
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
#' ## Create test interactions
#' gi <- read.table(text="
#'             1 51000000 51300000 1 51000000 51500000
#'             2 52000000 52300000 3 52000000 52500000
#'             1 150000000 150500000 1 150000000 150300000
#'             2 52000000 52300000 2 52000000 52800000") |>
#'     as_ginteractions()
#'
#' ## InteractionJaggedArray object
#' iarr <- pullHicMatrices(gi, hicFiles, 100e03, half="both")
#'
#' ## Shift first two ranges out of range
#' gi2 <- c(binPairs(gi[1:2], binSize=100e3, pos1=-200e3), gi[3:4])
#'
#' ## Find overlaps
#' findOverlaps(iarr, gi2)
#' countOverlaps(iarr, gi2)
#' countOverlaps(iarr, gi2, maxgap=100e3)
#' overlapsAny(iarr, gi2)
#' subsetByOverlaps(iarr, gi2)
#' subsetByOverlaps(iarr, gi2, invert=TRUE)
#'
#' @name findOverlaps
#' @rdname InteractionJaggedArray-overlaps
NULL

for (sig in c("InteractionJaggedArray", "Vector", "missing")) {

    #' @exportMethod findOverlaps
    setMethod(
        "findOverlaps",
        signature(query="InteractionJaggedArray", subject=sig),
        function(query, subject, maxgap=-1L, minoverlap=0L,
                 type=c("any", "start", "end", "within", "equal"),
                 select=c("all", "first", "last", "arbitrary"),
                 ignore.strand=TRUE, ..., use.region="both"){
            call <- as.list(match.call(expand.dots=TRUE))[-1]
            call$query <- interactions(query)
            if (missing(subject)){
                call[["subject"]] <- NULL
            }
            else if (is(subject, class(query))) {
                call$subject <- interactions(subject)
            }
            do.call(findOverlaps, call)
        }
    )

    #' @exportMethod countOverlaps
    setMethod(
        "countOverlaps",
        signature(query="InteractionJaggedArray", subject=sig),
        function(query, subject, maxgap=-1L, minoverlap=0L,
                 type=c("any", "start", "end", "within", "equal"),
                 select=c("all", "first", "last", "arbitrary"),
                 ignore.strand=TRUE, ..., use.region="both"){
            call <- as.list(match.call(expand.dots=TRUE))[-1]
            call$query <- interactions(query)
            if (missing(subject)){
                call[["subject"]] <- NULL
            }
            else if (is(subject, class(query))) {
                call$subject <- interactions(subject)
            }
            hits <- do.call(findOverlaps, call)
            ans <- countQueryHits(hits)
            names(ans) <- names(call$query)
            ans
        }
    )

    #' @exportMethod overlapsAny
    setMethod(
        "overlapsAny",
        signature(query="InteractionJaggedArray", subject=sig),
        function(query, subject, maxgap=-1L, minoverlap=0L,
                 type=c("any", "start", "end", "within", "equal"),
                 ..., use.region="both") {
            call <- as.list(match.call(expand.dots=TRUE))[-1]
            call$query <- interactions(query)
            if (missing(subject)){
                call[["subject"]] <- NULL
            }
            else if (is(subject, class(query))) {
                call$subject <- interactions(subject)
            }
            do.call(overlapsAny, call)
        }
    )

    #' @exportMethod subsetByOverlaps
    setMethod(
        "subsetByOverlaps",
        signature(x="InteractionJaggedArray", ranges=sig),
        function(x, ranges, maxgap=-1L, minoverlap=0L,
                 type=c("any", "start", "end", "within", "equal"),
                 invert=FALSE, ..., use.region="both") {
            call <- as.list(match.call(expand.dots=TRUE))[-1]
            call$x <- interactions(x)
            if (missing(ranges)){
                call[["ranges"]] <- NULL
            }
            else if (is(ranges, class(x))) {
                call$ranges <- interactions(ranges)
            }
            names(call)[names(call)=="x"] <- "query"
            names(call)[names(call)=="ranges"] <- "subject"
            call[["invert"]] <- NULL
            ov_any <- do.call(overlapsAny, call)
            if (invert) ov_any <- !ov_any
            x[which(ov_any),]
        }
    )
}


