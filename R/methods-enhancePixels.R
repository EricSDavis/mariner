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
#'
#' @importFrom abind abind
#' @importFrom S4Vectors metadata first second
#' @importFrom SummarizedExperiment assay
#' @importFrom data.table as.data.table
#' @importFrom InteractionSet GInteractions interactions
#' @importFrom IRanges IRanges
#'
#' @returns A GInteractions object with the updated
#'  pixel interactions, along with a column with the
#'  aggregated max/min value for that pixel.
#'
#' @rdname selectPixel
#' @export
selectPixel <- function(x,
                        aggFUN=sum,
                        selectFUN="which.max",
                        nBlocks=5,
                        verbose=TRUE) {

    ## Parse arguments
    aggFUN <- match.fun(aggFUN)
    selectFUN <- match.fun(selectFUN)

    ## Define binSize from input data
    binSize <- metadata(x)$binSize

    ## Sum across Hi-C files (by interactions)
    aggArr <- aggHicMatrices(x, "interactions", aggFUN, nBlocks, verbose)

    ## Find pixel of interest for each interaction
    poi <- apply(aggArr, 3, selectFUN)

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
    pixels <- sapply(seq_len(nrow(data)), \(i) {
        data[i, poi[i]]
    })

    ## Select accompanying values
    # values <- sapply(seq_len(nrow(data)), \(i) {
    #     aggArr
    # })

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
            )
        )

    return(gi)
}




selectMaxPixel <- function(x) {

    ## Set binSize and data.table
    binSize <- metadata(x)$binSize
    dt <- data.table::as.data.table(interactions(x))

    ## Initialize progress bar
    pb <- progress_bar$new(
        format = "  :step [:bar] :percent elapsed: :elapsedfull",
        clear = FALSE,
        total = 4
    )
    step <- list(step="Selecting max pixel")
    pb$tick(0)

    ## Generate bins for start1 and start2
    pb$tick(tokens=step)
    bins <- dt[,{
        .(x = seq(start1, end1-binSize, binSize),
          y = seq(start2, end2-binSize, binSize))
    }, by = .(groupRow=1:nrow(dt))]

    ## Expand bins to long format (fast cross-join)
    pb$tick(tokens=step)
    longMat <- bins[, {
        do.call(CJ, c(.SD, sorted = F))
    }, .SDcols = c('x', 'y'), by = groupRow]

    ## Add counts for each hic file
    pb$tick(tokens=step)
    cn <- c()
    for(i in seq_len(dim(x)[2])) {
        cn[i] <- paste0("counts",i)
        longMat[,cn[i]] <-
            as.vector(aperm(counts(x), c(2,1,3,4))[,,,i])
    }

    ## Sum counts
    longMat[,countSum := rowSums(.SD), .SDcols=cn]

    ## Subset for max pixels
    pixels <- longMat[longMat[,.I[which.max(countSum)],by=groupRow]$V1]

    ## Close progress bar
    pb$tick(tokens = list(step = "Done!"))
    if(pb$finished) pb$terminate()

    ## Convert pixels to GInteractions Object
    gi <-
        GInteractions(
            anchor1 = GRanges(
                seqnames = seqnames(first(interactions(x))),
                ranges = IRanges(start = pixels$x,
                                 width = binSize+1)
            ),
            anchor2 = GRanges(
                seqnames = seqnames(second(interactions(x))),
                ranges = IRanges(start = pixels$y,
                                 width = binSize+1)
            ),
            countSum = pixels$countSum
        )


    return(gi)

}


enhancePixels <- function(x, files, from, to) {
    x |>
        binPairs(binSize = from) |>
        pullHicMatrices(binSize = to, files = files) |>
        selectMaxPixel()
}
