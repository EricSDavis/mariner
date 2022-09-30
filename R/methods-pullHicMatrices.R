#' Ensure pairs are binned and warn if not
#' @inheritParams pullHicMatrices
#' @importFrom rlang inform
#' @importFrom glue glue
#' @return GInteractions object that has been
#'  binned (or original object if already binned).
#' @noRd
.handleBinning <- function(x, binSize) {

    ## Check for binned GInteractions
    binned <- .checkSnappedPairs(x, binSize)

    ## Inform user and snap GInteractions to bins
    if (!binned) {
        x <- snapToBins(x, binSize)
        msg <- c("Pairs are not binned to `binSize`.",
                 'i' = glue("Snapping to binSize={binSize}, "),
                 'i' = glue("Use `binPairs()` or `snapToBins()` ",
                            "for more binning options."))
        inform(msg)
    }

    ## Return appropriately binned object
    return(x)
}

#' Check
#'
#' @inheritParams pullHicMatrices
#' @importFrom GenomeInfoDb seqinfo seqnames
#' @importFrom strawr readHicChroms
#' @importFrom rlang abort
#' @importFrom glue glue
#' @returns Error if there is a chromosome issue
#' @noRd
.checkHicChroms <- function(x, file) {

    ## Extract chromosomes from x and file
    chromsInX <- seqnames(seqinfo(x))
    chromsInFile <- readHicChroms(file)$name

    ## Ensure all chromosomes in x are in file
    if (!all(chromsInX %in% chromsInFile)) {
        abort(c(
            "There's some chr-chr-craziness going on.",
            'x' = glue("seqnames in `x` are not correctly ",
                       "formatted or do not exist in `file`."),
            "Try the following steps:",
            '>' = "Check `x` with `GenomeInfoDb::seqinfo(x)`.",
            '>' = "Check `file` with `strawr::readHicChroms(file)`.",
            '>' = "Edit seqnames in `x` to match chromosomes in `file`.",
            "Hint:",
            '>' = glue("`GenomeInfoDb::seqlevelsStyle(x)",
                       " <- 'UCSC'` for 'chr' prefix."),
            '>' = glue("`GenomeInfoDb::seqlevelsStyle(x)",
                       " <- 'ENSEMBL'` without 'chr' prefix.")
        ))
    }
}


#' Convert GInteractions to sorted short format
#'
#' Extra short format is the minimal information
#' needed to extract Hi-C contacts with `strawr`.
#' See the description for this format here:
#' https://github.com/aidenlab/juicer/wiki/Pre#extra-short-format-dev.
#'
#' It also orders the interactions and chromosomes
#' appropriately. For intrachromosomal interactions, the
#' anchors should be flipped such that start1 < start2.
#' For interchromosomal interactions, columns should be
#' flipped so that chr1 < chr2 (according to the .hic
#' file's internal chromosome map index).
#'
#' @inheritParams pullHicMatrices
#' @importFrom data.table as.data.table `:=`
#' @importFrom strawr readHicChroms
#' @return `data.table` with columns:
#'  "chr1", "start1", "chr2", "start2".
#' @noRd
.GInteractionsToShortFormat <- function(x, file) {

    ## Convert to data.table format
    x <-
        as.data.table(x)[, c("seqnames1", "start1",
                             "seqnames2", "start2")] |>
        `colnames<-`(c("seqnames1", "pos1",
                       "seqnames2", "pos2"))

    ## Get strawr chromosome map index
    chrMap <- readHicChroms(file)

    ## Get indices for correct ordering in .hic file
    x$chromIndex1 <- chrMap$index[match(x$seqnames1, chrMap$name)]
    x$chromIndex2 <- chrMap$index[match(x$seqnames2, chrMap$name)]

    ## Interchromosomal: Flip column order
    ## so that seqnames1 < seqnames2
    x[chromIndex1 > chromIndex2,
      `:=`(chromIndex1=chromIndex2,
           chromIndex2=chromIndex1,
           seqnames1=seqnames2, pos1=pos2,
           seqnames2=seqnames1, pos2=pos1)]

    ## Intrachromosmal: Flip column order
    ## so that pos1 < pos2
    x[chromIndex1 == chromIndex2 & pos1 > pos2,
      `:=`(pos1=pos2, pos2=pos1)]

    ## Remove indices from data.table
    x[,c("chromIndex1", "chromIndex2") := NULL]

    return(x)
}

#' Replacement for "GInteractionsToShortFormat"
#' @noRd
.orderInteractions <- function(x, file) {

    ## Convert to data.table format
    x <-
        as.data.table(x)[, c("seqnames1", "start1", "end1",
                             "seqnames2", "start2", "end2")]

    ## Get strawr chromosome map index
    chrMap <- readHicChroms(file)

    ## Get indices for correct ordering in .hic file
    x$chromIndex1 <- chrMap$index[match(x$seqnames1, chrMap$name)]
    x$chromIndex2 <- chrMap$index[match(x$seqnames2, chrMap$name)]

    ## Interchromosomal: Flip column order
    ## so that seqnames1 < seqnames2
    x[chromIndex1 > chromIndex2,
      `:=`(chromIndex1=chromIndex2,
           chromIndex2=chromIndex1,
           seqnames1=seqnames2, start1=start2, end1=end2,
           seqnames2=seqnames1, start2=start1, end2=end1)]

    ## Intrachromosmal: Flip column order
    ## so that start1 < start2
    x[chromIndex1 == chromIndex2 & start1 > start2,
      `:=`(start1=start2, start2=start1,
           end1=end2, end2=end1)]

    ## Remove indices from data.table
    x[,c("chromIndex1", "chromIndex2") := NULL]

    return(x)
}

#' Put ranges into blocks
#'
#' This is different than binning ranges because
#' it modifies the blocks to fit the ranges, rather
#' than modifiying the ranges to fit "bins".
#'
#' @param start Integer vector of start positions.
#' @param end Integer vector of end positions.
#' @param blockSize Numeric size of blocks.
#'
#' @importFrom rlang abort
#'
#' @returns List of blocked start and end vectors
#'
#' @noRd
.blockRanges <- function(start, end, blockSize) {

    ## Get regions in terms of blockSize
    s <- start / blockSize
    e <- end / blockSize

    ## Throw error if any ranges are greater than
    ## blockSize/2
    widthsFail <- end-start > blockSize/2
    if (any(widthsFail)) {
        abort(c(
            glue("Each `end` - `start` must be ",
                 "less than `blockSize/2`"),
            "i" = glue("The following indices do not ",
                       "meet this criteria: ",
                       glue_collapse(which(widthsFail), ", "))
        ))
    }

    ## Calculate bins covered for each anchor
    bc <- floor(e) - floor(s)

    ## Shift ranges that intersect a boundary by
    ## half a blockSize
    s[bc == 1] <- floor(s[bc == 1])*blockSize + blockSize/2
    e[bc == 1] <- floor(e[bc == 1])*blockSize + blockSize/2

    ## Ranges that fall within blocks are
    ## expanded to blockSize
    s[bc == 0] <- floor(s[bc == 0])*blockSize
    e[bc == 0] <- (floor(e[bc == 0])+1)*blockSize

    return(list(start = s, end = e))
}

#' Create 2D-blocks and map paired ranges to them
#'
#' Creates a unique set of two-dimenional blocks
#' of width `blockSize` and maps GInteractions
#' pairs to each block.
#'
#' @param x GInteractions object.
#' @param blockSize Numeric size of blocks.
#'
#' @importFrom rlang warn
#' @importFrom InteractionSet regions anchors
#' @importFrom GenomicRanges width
#' @importFrom glue glue
#' @importFrom data.table data.table .I .GRP
#'
#' @returns GInteractions object describing
#'  blocks of data with a column, `xIndex`
#'  connecting them to the input ranges, `x`.
#'
#' @noRd
.mapToBlocks <- function(x, blockSize = 10e06) {

    ## Ensure that the chosen blockSize is
    ## twice the size of largest range
    maxRangeWidth <- max(width(regions(x)))
    if (maxRangeWidth > blockSize/2) {
        blockSize <- maxRangeWidth*2
        warn(c(
            "x" = glue("`blockSize` must be twice the ",
                       "longest range."),
            'i' = glue("Setting `blockSize` = {blockSize}.")
        ))
    }

    ## Extract anchors
    a1 <- anchors(x, 'first')
    a2 <- anchors(x, 'second')

    ## Assign ranges to blocks
    r1 <- .blockRanges(start(a1), end(a1), blockSize)
    r2 <- .blockRanges(start(a2), end(a2), blockSize)

    ## Assemble into data.table
    dt <-
        data.table(seqnames1 = as.character(seqnames(a1)),
                   blockStart1 = r1$start,
                   blockEnd1 = r1$end,
                   seqnames2 = as.character(seqnames(a2)),
                   blockStart2 = r2$start,
                   blockEnd2 = r2$end)

    ## Define unique set of blocks
    blocks <-
        dt[, .(block=.GRP, xIndex=.(.I)), by=names(dt)] |>
        as_ginteractions()

    return(blocks)
}

#' Internal pullHicMatrices
#' @inheritParams pullHicMatrices
#' @return Array of Hi-C submatrices.
#' @noRd
.pullHicMatrices <- function(x, binSize, file) {

    ## TODO: Check that binSize is contained in
    ## file

    ## TODO: Change this to "snapToHicBins"
    ## to enforce the starts and ends of ranges
    ## going to acceptable binSizes in the file
    ## Bin GInteractions if necessary
    x <- .handleBinning(x, binSize)

    ## Ensure seqnames are properly formatted
    .checkHicChroms(x, file)

    ## Convert to short format and sort interactions
    x <- .GInteractionsToShortFormat(x, file)

    ## TODO: Develop method for generating clusters
    ## interactions to extract from file.
    .clusterPairs(x=x,
                  radius = binSize,
                  method = 'manhattan',
                  pos = 'center')



}

#' Pull matrices from a `.hic` file
#'
#' @param x GInteractions object.
#' @param binSize Integer (numeric) describing the
#'  resolution (range widths) of the paired data.
#' @param file Character file path to `.hic` file.
#' @return Array of Hi-C submatrices.
#' @noRd
# #' @rdname pullHicMatrices
# #' @export
# setMethod("pullHicMatrices",
#           signature(x = 'GInteractions',
#                     binSize = 'numeric',
#                     file = 'character'),
#           definition = .pullHicMatrices)
