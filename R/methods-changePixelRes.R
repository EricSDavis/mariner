#' Internal select pixel function
#' @inheritParams selectPixel
#' @importFrom abind abind
#' @importFrom S4Vectors metadata first second
#' @importFrom SummarizedExperiment assay
#' @importFrom data.table as.data.table
#' @importFrom InteractionSet GInteractions interactions
#' @importFrom IRanges IRanges
#' @returns A GInteractions object with the updated
#'  pixel interactions, along with a column with the
#'  aggregated max/min value for that pixel.
#' @noRd
.selectPixel <- function(x, aggFUN, selectFUN, nBlocks, verbose) {

    ## Parse arguments
    aggFUN <- match.fun(aggFUN)
    selectFUN <- match.fun(selectFUN)

    ## Define binSize from input data
    binSize <- metadata(x)$binSize

    ## Sum across Hi-C files (by interactions)
    aggArr <- aggHicMatrices(x, "interactions", aggFUN, nBlocks, verbose)

    ## Find pixel of interest for each interaction
    poi <- apply(aggArr, 3, selectFUN)
    val <- apply(aggArr, 3, \(x) x[selectFUN(x)])

    ## Get indices for rows/cols
    idx <-
        expand.grid(seq_len(dim(counts(x))[1]),
                    seq_len(dim(counts(x))[2])) |>
        `colnames<-`(value = c("row", "col"))

    ## Pull row/cols & realize as matrices
    ## and paste together
    rows <- assay(x, 'rownames')[,1,][,idx$row] |> as.matrix()
    cols <- assay(x, 'colnames')[,1,][,idx$col] |> as.matrix()
    data <- apply(abind(rows, cols, along=3), c(1,2),
                  paste0, collapse="-")

    ## Select pixels of interest
    pixels <- vapply(seq_len(nrow(data)), \(i) {
        data[i, poi[i]]
    }, FUN.VALUE=character(1L))

    ## Split character ranges into
    pixelDF <-
        strsplit(pixels, "-") |>
        do.call(rbind, args=_) |>
        apply(2, as.numeric) |>
        as.data.table() |>
        `colnames<-`(value=c("start1", "start2"))

    ## Convert pixels to GInteractions Object
    gi <-
        GInteractions(
            anchor1 = GRanges(
                seqnames = seqnames(first(interactions(x))),
                ranges = IRanges(start = pixelDF$start1,
                                 width = binSize+1)
            ),
            anchor2 = GRanges(
                seqnames = seqnames(second(interactions(x))),
                ranges = IRanges(start = pixelDF$start2,
                                 width = binSize+1)
            ),
            value = val
        )

    return(gi)
}

#' Get the pixel representing the strongest
#' or weakest interaction in an InteractionArray
#'
#' @param x InteractionArray object
#' @param aggFUN Function to use for aggregating
#'  across Hi-C files. Must be passable to
#'  `which.max` or `which.min`. Default is "sum".
#' @param selectFUN Function to use for selecting
#'  among aggregated interactions. Must be one of
#'  "which.max" or "which.min".
#' @param nBlocks Number of blocks for block-processing
#'  arrays. Default is 5. Increase this for large
#'  datasets. To read and process all data at once, set
#'  this value to 1.
#' @param verbose Boolean (TRUE or FALSE) describing
#'  whether to report block-processing progress. Default
#'  is TRUE.
#'
#' @examples
#' ## Read .hic file paths
#' hicFiles <-
#'     system.file("extdata/test_hic", package="mariner") |>
#'     list.files(pattern=".hic", full.names=TRUE)
#'
#' ## Read in loops as GInteractions object
#' loops <-
#'     system.file("extdata", package="mariner") |>
#'     list.files(pattern="WT.*Loops.txt", full.names=TRUE) |>
#'     read.table(header=TRUE) |>
#'     as_ginteractions(keep.extra.columns=FALSE)
#'
#' ## Removes the "chr" prefix for compatibility
#' ## with the preprocessed hic files
#' GenomeInfoDb::seqlevelsStyle(loops) <- 'ENSEMBL'
#'
#' ## Rebin loops to 2.5e6 resolution
#' loops <- binPairs(x=loops, binSize=2.5e06)
#'
#' ## Pull 5x5 matrices
#' iarr <- pullHicMatrices(x=loops[1:5],
#'                         files=hicFiles,
#'                         binSize=500e3,
#'                         norm="KR",
#'                         half='upper')
#'
#' ## Select pixel
#' selectPixel(iarr)
#'
#' @returns A GInteractions object with the updated
#'  pixel interactions, along with a column with the
#'  aggregated max/min value for that pixel.
#'
#' @rdname selectPixel
#' @export
setMethod("selectPixel",
          signature(x="InteractionArray"),
          definition=.selectPixel)


#' Internal changePixelRes
#' @inheritParams changePixelRes
#' @returns A GInteractions object with the updated
#'  pixel interactions, along with a column with the
#'  aggregated max/min value for that pixel.
#' @noRd
.changePixelRes <- function(x, files, from, to,
                            aggFUN=sum,
                            selectFUN="which.max",
                            nBlocks=5,
                            verbose=TRUE,
                            norm="KR",
                            half="upper",
                            ...) {
    ## Basic type checking
    .checkTypes(types=list(from="number", to="number"))

    ## Execute steps
    x |>
        binPairs(binSize=from) |>
        pullHicMatrices(binSize=to, files=files, norm=norm, ...) |>
        selectPixel(aggFUN, selectFUN, nBlocks, verbose)
}


#' Change pixels from one resolution to another
#' selecting the new pixel using Hi-C data.
#'
#' A GInteractions object containing pixels of
#' interest is resized to the `from` resolution
#' (if its not already), then count matrices are
#' extracted for each interaction and Hi-C file
#' using the new `to` resolution. Count matrices
#' are aggregated by interactions with the
#' supplied `aggFUN`, and a new pixel is selected
#' with the supplied `selectFUN`. For large
#' datasets, increase `nBlocks` to allow for smaller
#' blocks of data to be processed in memory.
#'
#' @param x GInteractions object.
#' @param files Character file paths to `.hic` files.
#' @param from Number (length one numeric vector) describing
#'  the resolution of `x`. Data will be binned to this
#'  value if it is not already binned.
#' @param to Number (length one numeric vector) describing
#'  the new resolution for the pixels.
#' @param aggFUN Function to use for aggregating
#'  across Hi-C files. Must be passable to
#'  `which.max` or `which.min`. Default is "sum".
#' @param selectFUN Function to use for selecting
#'  among aggregated interactions. Must be one of
#'  "which.max" or "which.min".
#' @param nBlocks Number of blocks for block-processing
#'  arrays. Default is 5. Increase this for large
#'  datasets. To read and process all data at once, set
#'  this value to 1.
#' @param verbose Boolean (TRUE or FALSE) describing
#'  whether to report block-processing progress. Default
#'  is TRUE.
#' @param norm String (length one character vector)
#'  describing the Hi-C normalization to apply. Use
#'  `strawr::readHicNormTypes()` to see accepted values
#'  for each file in `files`.
#' @param half String (character vector of length one)
#'  indicating whether to keep values for the upper
#'  triangular (`half="upper"`) where `start1 < start2`,
#'  lower triangular (`half="lower"`) where
#'  `start1 > start2`, or both (`half="both"`, default).
#'  When `half="upper"` all lower triangular values are `NA`.
#'  When `half="lower"` all upper triangular values are `NA`.
#'  When `half="both"` there are no `NA` values.
#'  For interchromosomal interactions there is no inherent
#'  directionality between chromosomes, so data is returned
#'  regardless of specified order.
#' @param ... Additional arguments passed to `pullHicMatrices()`.
#'  See ?[`pullHicMatrices`].
#'
#' @examples
#' ## Read .hic file paths
#' hicFiles <-
#'     system.file("extdata/test_hic", package="mariner") |>
#'     list.files(pattern=".hic", full.names=TRUE)
#'
#' ## Read in loops as GInteractions object
#' loops <-
#'     system.file("extdata", package="mariner") |>
#'     list.files(pattern="WT.*Loops.txt", full.names=TRUE) |>
#'     read.table(header=TRUE) |>
#'     as_ginteractions(keep.extra.columns=FALSE)
#'
#' ## Removes the "chr" prefix for compatibility
#' ## with the preprocessed hic files
#' GenomeInfoDb::seqlevelsStyle(loops) <- 'ENSEMBL'
#'
#' ## Rebin loops to 2.5e6 resolution
#' loops <- binPairs(x=loops, binSize=2.5e06)
#'
#' ## Change pixel resolution from 2.5e6 to 500e3
#' changePixelRes(x=loops[1:5],
#'                files=hicFiles,
#'                from=2.5e6,
#'                to=500e3)
#'
#' @returns A GInteractions object with the updated
#'  pixel interactions, along with a column with the
#'  aggregated max/min value for that pixel.
#'
#' @rdname changePixelRes
#' @export
setMethod("changePixelRes",
          signature(x="GInteractions",
                    files="character"),
          definition=.changePixelRes)
