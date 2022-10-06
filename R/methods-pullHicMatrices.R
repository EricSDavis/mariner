#' Check straw arguments are valid against files
#'
#' @inheritParams pullHicMatrices
#' @importFrom strawr readHicNormTypes readHicBpResolutions
#' @importFrom rlang abort arg_match
#' @importFrom glue glue
#' @noRd
.checkStrawArgs <- function(files, norm, binSize, matrix) {

    ## Check norm and binSize against each file
    for(f in files) {
        arg_match(norm, readHicNormTypes(f))
        if (!binSize %in% as.integer(readHicBpResolutions(f))) {
            abort(c(
                glue("binSize={binSize} is not valid."),
                "i"= glue("Use `readHicBpResolutions()` to \\
                          see allowed values.")
            ))
        }
    }

    ## Check matrix argument
    arg_match(matrix, c("observed", "expected", "oe"))
    return(NULL)
}


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

#' Check that chromosome maps in files
#'
#' Ensure they are identical for all files
#' and all chromosomes in `x` are contained
#' in these files.
#'
#' @inheritParams pullHicMatrices
#' @importFrom GenomeInfoDb seqinfo seqnames
#' @importFrom strawr readHicChroms
#' @importFrom rlang abort
#' @importFrom glue glue
#' @returns Error if there is a chromosome issue
#' @noRd
.checkHicChroms <- function(x, files) {

    ## Ensure all chromosome maps are identical
    ## if this fails, user will have to run
    ## the function in separate calls
    chrMaps <- lapply(files, \(f) readHicChroms(f))
    if (!all(sapply(chrMaps, identical, chrMaps[[1]]))) {
        abort(c(
            "Chromosome maps in `files` are not identical",
            "i"="Check this with `strawr::readHicChroms(files[1])`",
            "*"=glue("This is essential to ensure interactions ",
                     "are ordered correctly."),
            "*"=glue("Call this function multiple times ",
                     "or reprocess the Hi-C maps in the ",
                     "same way to proceed.")
        ))
    }

    ## Extract chromosomes from x and files
    chromsInX <- seqnames(seqinfo(x))
    chromsInFile <- chrMaps[[1]]$name

    ## Ensure all chromosomes in x are in file
    if (!all(chromsInX %in% chromsInFile)) {
        abort(c(
            "There's some chr-chr-craziness going on.",
            'x' = glue("seqnames in `x` are not correctly ",
                       "formatted or do not exist in `files`."),
            "Try the following steps:",
            '>' = "Check `x` with `GenomeInfoDb::seqinfo(x)`.",
            '>' = "Check `file` with `strawr::readHicChroms(files)`.",
            '>' = "Edit seqnames in `x` to match chromosomes in `files`.",
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
.GInteractionsToShortFormat <- function(x, files) {

    ## Convert to data.table format
    x <-
        as.data.table(x)[, c("seqnames1", "start1",
                             "seqnames2", "start2")] |>
        `colnames<-`(c("seqnames1", "pos1",
                       "seqnames2", "pos2"))

    ## Get strawr chromosome map index
    chrMap <- readHicChroms(files)

    ## Get indices for correct ordering in .hic files
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

    ## Return as GInteractions
    return(as_ginteractions(x))
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

#' Pull data from files and return as
#' InteractionSet with DelayedArrays
#'
#' @inheritParams pullHicMatrices
#' @param mDim integer - dimensions for matrices (equal)
#' @param blocks
#'
#' @importFrom rhdf5 h5createFile h5createDataset h5write
#' @importFrom progress progress_bar
#' @importFrom strawr straw
#' @importFrom data.table as.data.table setkeyv CJ
#' @importFrom InteractionSet InteractionSet
#' @importFrom DelayedArray DelayedArray
#' @importFrom HDF5Array HDF5Array
#'
#' @returns Interaction set with DelayedArrays of
#'  extracted data.
#'
#' @noRd
.pullArray <- function(x, binSize, files, norm, matrix,
                       blockSize, onDisk, compressionLevel,
                       chunkSize, mDim, blocks) {
    ## Determine dimensions for dataset
    ## Dim order is nInteractions, nFiles, matrix dims
    dims <- c(length(x), length(files), mDim, mDim)

    if (onDisk) {
        ## Create hdf5 for storage
        h5 <- tempfile(fileext = ".h5")
        h5createFile(h5)

        ## Set default chunkSize
        if (missing(chunkSize)) {
            if (compressionLevel >= 5) {
                chunkSize <- 1
            } else {
                chunkSize <- length(x)
            }
        }

        ## Create dataset for holding array of counts
        ## and row/colnames
        h5createDataset(file = h5,
                        dataset = "counts",
                        dims = dims,
                        chunk = c(chunkSize, 1, dims[c(3,4)]),
                        storage.mode = "double",
                        fillValue = NA_real_,
                        level = compressionLevel)
        h5createDataset(file = h5,
                        dataset = "rownames",
                        dims = dims[c(1,2,3)],
                        chunk = c(chunkSize, 1, dims[3]),
                        storage.mode = "integer",
                        fillValue = NA_integer_,
                        level = compressionLevel)
        h5createDataset(file = h5,
                        dataset = "colnames",
                        dims = dims[c(1,2,4)],
                        chunk = c(chunkSize, 1, dims[4]),
                        storage.mode = "integer",
                        fillValue = NA_integer_,
                        level = compressionLevel)
    } else {
        counts <- array(data=NA, dim = dims)
        rownames <- array(data=NA_integer_, dim=dims[c(1,2,3)])
        colnames <- array(data=NA_integer_, dim=dims[c(1,2,4)])
    }

    ## Start progress bar
    pb <- progress_bar$new(
        format = "  :step [:bar] :percent elapsed: :elapsedfull",
        clear = FALSE, total = nrow(blocks)*length(files)+1)
    pb$tick(0)

    ## Begin extraction for each file and block
    for(j in seq_len(length(files))) {
        for(i in seq_len(nrow(blocks))) {

            ## Update progress
            pb$tick(tokens=list(step=sprintf(
                'Pulling block %s of %s from file %s of %s',
                i, nrow(blocks), j, length(files)
            )))

            ## Extract block data from file
            sparseMat <-
                straw(norm = norm,
                      fname = files[j],
                      chr1loc = blocks$chr1loc[i],
                      chr2loc = blocks$chr2loc[i],
                      unit = "BP",
                      binsize = binSize,
                      matrix = matrix) |>
                as.data.table()

            ## Select interactions belonging to block
            xIndices <- blocks$xIndex[[i]]
            g <- as.data.table(x[xIndices])

            ## Generate bins for rows and columns
            bins <- g[,.(x = seq(start1, end1, binSize),
                         y = seq(start2, end2, binSize),
                         counts = 0),
                      by = .(groupRow = 1:nrow(g))]

            ## For pullHicPixels
            # bins <- g[,.(x = start1,
            #              y = start2,
            #              counts = 0),
            #           by = .(groupRow = 1:nrow(g))]

            ## Expand bins to long format (fast cross-join)
            longMat <- bins[, do.call(CJ, c(.SD, sorted = F)),
                            .SDcols = c('x', 'y'), by = groupRow]
            longMat$counts <- 0

            ## Rename columns and rearrange
            longMat <- longMat[,.(x, y, counts, groupRow)]

            ## Set keys
            setkeyv(sparseMat, c('x', 'y'))

            ## Get counts by key
            longMat$counts <- sparseMat[longMat]$counts

            ## Set unmatched counts (NA) to 0
            longMat[is.na(counts), counts := 0]

            ## Collect row/colname info
            aInfo <- longMat[,.(rowNames = .(unique(x)),
                                colNames = .(unique(y))),
                             by=groupRow]

            ## Fill array by col, row, interactions, file
            ## then rearrange for storage: interactions, file, rows, cols
            a <- array(data=longMat$counts,
                       dim=c(dims[4], dims[3], length(xIndices), 1)) |>
                aperm(c(3, 4, 2, 1))

            ## Fill row/colnames by row/col, interactions, file
            ## then rearrange for storage: interactions, file, row/col
            rn <- array(unlist(aInfo$rowNames),
                        dim=c(mDim, length(xIndices), 1)) |>
                aperm(c(2,1,3))
            cn <- array(unlist(aInfo$colNames),
                        dim=c(mDim, length(xIndices), 1)) |>
                aperm(c(2,1,3))

            if (onDisk) {
                # Write data to hdf5 file
                h5write(obj = a,
                        file = h5,
                        name = "counts",
                        index=list(xIndices, j,
                                   seq_len(dims[3]),
                                   seq_len(dims[4])))
                h5write(obj = rn,
                        file = h5,
                        name = "rownames",
                        index=list(xIndices, j, seq_len(dims[3])))
                h5write(obj = cn,
                        file = h5,
                        name = "colnames",
                        index=list(xIndices, j, seq_len(dims[4])))
            } else {
                counts[xIndices, j, seq_len(dims[3]), seq_len(dims[4])] <- a
                rownames[xIndices, j, seq_len(dims[3])] <- rn
                colnames[xIndices, j, seq_len(dims[3])] <- cn
            }

        }
    }

    ## Close progress bar
    pb$tick(tokens = list(step = "Done!"))
    if(pb$finished) pb$terminate()

    ## Construct InteractionSet object
    if (onDisk) {
        iset <-
            InteractionSet(
                interactions = x,
                assays = list(
                    counts = DelayedArray(HDF5Array(h5,"counts")),
                    rownames = DelayedArray(HDF5Array(h5,"rownames")),
                    colnames = DelayedArray(HDF5Array(h5,"colnames"))
                )
            )
        return(iset)

    } else {
        iset <-
            InteractionSet(
                interactions = x,
                assays = list(
                    counts = DelayedArray(counts),
                    rownames = DelayedArray(rownames),
                    colnames = DelayedArray(colnames)
                )
            )
        return(iset)
    }
}

#' Internal pullHicMatrices
#' @inheritParams pullHicMatrices
#' @importFrom data.table as.data.table
#' @importFrom InteractionSet regions
#' @importFrom S4Vectors width
#' @return Array of Hi-C submatrices.
#' @noRd
.pullHicMatrices <- function(x, binSize, files, norm, matrix,
                             blockSize, onDisk, compressionLevel,
                             chunkSize) {

    ## Check input types
    .checkTypes(c(
        binSize="number",
        norm="string",
        matrix="string",
        blockSize="number",
        onDisk="boolean",
        compressionLevel="number",
        chunkSize="number"
    ))

    ## Parse straw parameters
    .checkStrawArgs(files, norm, binSize, matrix)

    ## Ensure seqnames are properly formatted
    .checkHicChroms(x, files)

    ## Assign GInteractions to bins
    x <- .handleBinning(x, binSize)

    ## Order interactions according to chroms
    ## (ensures blocks will be ordered too)
    x <- .orderInteractions(x, files[1])

    ## Assign x to blocks
    blocks <- as.data.table(.mapToBlocks(x, blockSize))

    ## Paste ranges to get chr1 and chr2 locs
    blocks <-
        blocks[,.(chr1loc = paste(seqnames1, start1, end1, sep=":"),
                  chr2loc = paste(seqnames2, start2, end2, sep=":"),
                  block, xIndex)]

    ## TODO: Code will differ here depending on input cases
    ## Use region widths to dispatch code for
    ## extracting equal or variable dimension slices
    widths <- unique(width(regions(x))) - 1
    if (length(widths) == 1L) {

        ## Set matrix dimensions
        mDim <- (widths/binSize) + 1

        ## Dispatch .pullArray for equal dimensions
        iset <-
            .pullArray(x = x,
                       binSize = binSize,
                       files = files,
                       norm = norm,
                       matrix = matrix,
                       blockSize = blockSize,
                       onDisk = onDisk,
                       compressionLevel = compressionLevel,
                       chunkSize = chunkSize,
                       mDim = mDim,
                       blocks = blocks)
        return(iset)

    } else {
        ## TODO: Dispatch pullHicMatrices (varible dimensions)
        stop("Variable dimension matrices not currently supported.")
    }

}

#' Pull submatrices from `.hic` files
#'
#' @param x GInteractions object.
#' @param files Character file paths to `.hic` files.
#' @param binSize Integer (numeric) describing the
#'  resolution (range widths) of the paired data.
#' @param ... Additional arguments.
#' @param norm String (length one character vector)
#'  describing the Hi-C normalization to apply. Use
#'  `strawr::readHicNormTypes()` to see accepted values
#'  for each file in `files`.
#' @param matrix String (length one character vector)
#'  Type of matrix to extract. Must be one of "observed",
#'  "oe", or "expected". "observed" is observed counts,
#'  "oe" is observed/expected counts, "expected" is
#'  expected counts.
#' @param blockSize Number (length one numeric vector)
#'  describing the size in base-pairs to pull from each
#'  `.hic` file. Default is 248956422 (the length of the
#'  longest chromosome in the human hg38 genome). For
#'  large `.hic` files `blockSize` can be reduced to
#'  conserve the amount of data read in at a time. Larger
#'  `blockSize` values speed up performance, but use more
#'  memory.
#' @param onDisk Boolean (length one logical vector that
#'  is not NA) indicating whether extracted data should
#'  be stored on disk in an HDF5 file. Default is TRUE.
#' @param compressionLevel Number (length one numeric vector)
#'  between 0 (Default) and 9 indicating the compression
#'  level used on HDF5 file.
#' @param chunkSize Number (length one numeric vector)
#'  indicating how many values of `x` to chunk for each
#'  write to HDF5 stored data. This has downstream
#'  implications for accessing subsets later. For small
#'  `compressionLevel` values use smaller `chunkSize`
#'  values and for large `compressionLevel` values use large
#'  (i.e. `length(x)`) values to improve performance.
#'
#' @returns InteractionSet object with a 4-dimensional array
#'  of Hi-C submatrices, rownames, and colnames. Array is
#'  stored with the following dimensions: Interactions in `x`,
#'  Hi-C `files`, rows of submatrix, columns of submatrix.
#'
#' @rdname pullHicMatrices
#' @export
setMethod("pullHicMatrices",
          signature(x = 'GInteractions',
                    binSize = 'numeric',
                    files = 'character'),
          definition = .pullHicMatrices)

