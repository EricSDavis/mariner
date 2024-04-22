#' Check straw arguments are valid against files
#'
#' @inheritParams pullHicMatrices
#' @importFrom strawr readHicNormTypes readHicBpResolutions
#' @importFrom rlang abort arg_match
#' @importFrom glue glue
#' @noRd
.checkStrawArgs <- function(files, half, norm, binSize, matrix) {

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

    ## Check half and matrix arguments
    arg_match(half, c("both", "upper", "lower"))
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
                 'i' = glue("Use `assignToBins()` or `snapToBins()` ",
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
    if (!all(vapply(chrMaps, identical, chrMaps[[1]],
                    FUN.VALUE = logical(1L)))) {
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


#' Order interactions for pulling with strawr
#'
#' This function orders the interactions and chromosomes
#' appropriately. For intrachromosomal interactions, the
#' anchors should be flipped such that start1 < start2.
#' For interchromosomal interactions, columns should be
#' flipped so that chr1 < chr2 (according to the .hic
#' file's internal chromosome map index).
#'
#' Extra short format is the minimal information
#' needed to extract Hi-C contacts with `strawr`.
#' See the description for this format here:
#' https://github.com/aidenlab/juicer/wiki/Pre#extra-short-format-dev.
#'
#' @inheritParams pullHicMatrices
#' @importFrom data.table as.data.table `:=`
#' @importFrom strawr readHicChroms
#' @return `data.table` with columns:
#'  "chr1", "start1", "chr2", "start2".
#' @noRd
.orderInteractions <- function(x, file) {

    ## Suppress NSE notes in R CMD check
    chromIndex1 <- chromIndex2 <- NULL

    ## Convert to data.table format
    x <- as.data.table(x)[, c("seqnames1", "start1", "end1",
                              "seqnames2", "start2", "end2")]
    
    ## Get chromosome map index if cool file
    fileEnding <- tryCatch(.checkIfCool(file), 
               error = function(e){
                 return(FALSE)
               })
    
    if(fileEnding %in% c(".cool",".mcool")){
      chrMap <- readCoolChroms(file)
    } else {
      ## Get strawr chromosome map index if hic file
      chrMap <- readHicChroms(file)
    }
    
    ## Get indices for correct ordering in .hic/.cool file
    x$chromIndex1 <- chrMap$index[match(x$seqnames1, chrMap$name)]
    x$chromIndex2 <- chrMap$index[match(x$seqnames2, chrMap$name)]

    ## Interchromosomal: Flip column order
    ## so that seqnames1 < seqnames2
    chrSwapped <- x[chromIndex1>chromIndex2, which=TRUE]
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

    ## Return as GInteractions and
    ## indices of swapped chroms & anchors
    return(list(x=as_ginteractions(x),
                chrSwapped=chrSwapped))
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

    ## Assign ranges to blocks
    r1 <- .blockRanges(start1(x), end1(x), blockSize)
    r2 <- .blockRanges(start2(x), end2(x), blockSize)

    ## Assemble into data.table
    dt <-
        data.table(seqnames1 = seqnames1(x),
                   blockStart1 = r1$start,
                   blockEnd1 = r1$end,
                   seqnames2 = seqnames2(x),
                   blockStart2 = r2$start,
                   blockEnd2 = r2$end)

    ## Define unique set of blocks
    blocks <-
        dt[, .(block=.GRP, xIndex=.(.I)), by=names(dt)] |>
        as_ginteractions() |>
        suppressWarnings() # expected for interchromosomal

    return(blocks)
}

#' Assign counts to short format data.table
#' and flip interactions as requested to add
#' counts to each block data.
#' @param g data.table with group of interactions.
#' @param longMat long-format interactions in block.
#' @param sparseMat straw output for full block region.
#' @param chrSwapped which groups were previously swapped.
#' @param half Which half (both, upper or lower) to return.
#' @noRd
.assignCounts <- function(g, longMat, sparseMat, chrSwapped, half) {

    ## Suppress NSE notes in R CMD check
    x <- y <- grp <- NULL

    ## Set intrachromosomal flag
    intra <- identical(g$seqnames1, g$seqnames2)

    ## Swap where x > y for intrachromosomal (and keep track)
    if (intra)  {
        swapInd <- longMat[x > y, which=TRUE]
        longMat[swapInd, `:=`(x=y, y=x)]
    }
    if (!intra) {
        longMat[grp %in% chrSwapped, `:=`(x=y, y=x)]
    }

    ## Add counts by binary searching bins
    setkeyv(sparseMat, c('x', 'y'))
    longMat$counts <- sparseMat[longMat,counts,on=c('x','y')]

    ## Set unmatched counts to 0 and lower tri to NA
    longMat[is.na(counts), counts := 0]

    ## Swap back to original order &
    ## return requested half
    if (intra) {
        longMat[swapInd, `:=`(x=y,y=x)]
        if (identical(half, "upper")) {
            longMat[x > y, counts := NA_real_]
        }
        if (identical(half, "lower")) {
            longMat[x < y, counts := NA_real_]
        }
    }

    ## Swap back interactions for interchromosomal
    if (!intra) {
        longMat[grp %in% chrSwapped, `:=`(x=y, y=x)]
    }

    return(longMat)
}

#' Pull data from files and return as
#' InteractionSet with DelayedArrays
#'
#' @inheritParams pullHicMatrices
#' @param mDim integer - dimensions for matrices (equal)
#' @param blocks genomic blocks determining how much data to pull out at once
#' @param filesAreCool logical, TRUE if files are all .cool/.mcool
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
.pullArray <- function(x, files, binSize, h5File, half, norm,
                       matrix, blockSize, onDisk, compressionLevel,
                       chunkSize, mDim1, mDim2, blocks, chrSwapped,
                       filesAreCool) {
  
    ## Suppress NSE notes in R CMD check
    grp <- NULL

    ## Determine dimensions for dataset
    ## Dim order is nInteractions, nFiles, matrix dims
    dims <- c(length(x), length(files), mDim1, mDim2)

    if (onDisk) {
        ## Create hdf5 for storage
        h5createFile(h5File)

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
        h5createDataset(file = h5File,
                        dataset = "counts",
                        dims = dims,
                        chunk = c(chunkSize, 1, dims[c(3,4)]),
                        storage.mode = "double",
                        fillValue = NA_real_,
                        level = compressionLevel)
        h5createDataset(file = h5File,
                        dataset = "rownames",
                        dims = dims[c(1,2,3)],
                        chunk = c(chunkSize, 1, dims[3]),
                        storage.mode = "integer",
                        fillValue = NA_integer_,
                        level = compressionLevel)
        h5createDataset(file = h5File,
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
          if(filesAreCool){
            sparseMat <-
              coolStraw(norm = norm,
                    fname = files[j],
                    chr1loc = blocks$chr1loc[i],
                    chr2loc = blocks$chr2loc[i],
                    binsize = binSize) |>
              as.data.table()
          } else {
            sparseMat <-
                straw(norm = norm,
                      fname = files[j],
                      chr1loc = blocks$chr1loc[i],
                      chr2loc = blocks$chr2loc[i],
                      unit = "BP",
                      binsize = binSize,
                      matrix = matrix) |>
                as.data.table()
          }

            ## Select interactions belonging to block
            xIndices <- blocks$xIndex[[i]]
            g <- as.data.table(x[xIndices])
            g$xIndices <- xIndices

            ## Create submatrix bins for each range
            ## with a fast cross-join
            longMat <- g[,{
                x <- seq(start1, end1-binSize, binSize);
                y <- seq(start2, end2-binSize, binSize);
                x <- as.integer(x) #seq sometimes returns double
                y <- as.integer(y)
                CJ(x, y, sorted=FALSE)
            },
            by=.(grp=xIndices)]

            ## Assign counts
            longMat <- .assignCounts(g, longMat,
                                     sparseMat,
                                     chrSwapped,
                                     half)

            ## Collect row/colname info
            aInfo <- longMat[,.(rowNames = .(unique(x)),
                                colNames = .(unique(y))),
                             by=grp]

            ## Fill array by col, row, interactions, file
            ## then rearrange for storage: interactions, file, rows, cols
            a <- array(data=longMat$counts,
                       dim=c(dims[4], dims[3], length(xIndices), 1)) |>
                aperm(c(3, 4, 2, 1))

            ## Fill row/colnames by row/col, interactions, file
            ## then rearrange for storage: interactions, file, row/col
            rn <- array(unlist(aInfo$rowNames),
                        dim=c(dims[3], length(xIndices), 1)) |>
                aperm(c(2,1,3))
            cn <- array(unlist(aInfo$colNames),
                        dim=c(dims[4], length(xIndices), 1)) |>
                aperm(c(2,1,3))

            if (onDisk) {
                # Write data to hdf5 file
                h5write(obj = a,
                        file = h5File,
                        name = "counts",
                        index=list(xIndices, j,
                                   seq_len(dims[3]),
                                   seq_len(dims[4])))
                h5write(obj = rn,
                        file = h5File,
                        name = "rownames",
                        index=list(xIndices, j, seq_len(dims[3])))
                h5write(obj = cn,
                        file = h5File,
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

    ## Build colData and metadata
    colData <- DataFrame(files=files, fileNames=basename(files))
    metadata <- list(
        binSize=binSize,
        norm=norm,
        matrix=matrix
    )

    ## Construct InteractionArray object
    if (onDisk) {
        iset <-
            InteractionArray(
                interactions = x,
                assays = list(
                    counts = DelayedArray(HDF5Array(h5File,"counts")),
                    rownames = DelayedArray(HDF5Array(h5File,"rownames")),
                    colnames = DelayedArray(HDF5Array(h5File,"colnames"))
                ),
                colData = colData,
                metadata = metadata
            )
    } else {
        iset <-
            InteractionArray(
                interactions = x,
                assays = list(
                    counts = DelayedArray(counts),
                    rownames = DelayedArray(rownames),
                    colnames = DelayedArray(colnames)
                ),
                colData = colData,
                metadata = metadata
            )
    }

    ## Assign dimnames to files (default is basename(files))
    if (is.null(names(files))) {
        dimnames(iset)[2] <- list(iset$fileNames)
    } else {
        dimnames(iset)[2] <- list(names(files))
    }

    return(iset)
}

#' Pull data from files and return as
#' InteractionJaggedArray
#'
#' @inheritParams pullHicMatrices
#' @param mDim integer - dimensions for matrices (not equal)
#' @param blocks genomic blocks determining how much data to pull out at once
#' @param filesAreCool logical, TRUE if files are all .cool/.mcool
#'
#' @importFrom rhdf5 h5createFile h5createDataset h5write
#' @importFrom progress progress_bar
#' @importFrom strawr straw
#' @importFrom data.table as.data.table setkeyv CJ
#' @importFrom InteractionSet InteractionSet
#' @importFrom DelayedArray DelayedArray
#' @importFrom HDF5Array HDF5Array
#'
#' @returns InteractionJaggedArray Object.
#'
#' @noRd
.pullJaggedArray <- function(x, files, binSize, h5File, half, norm,
                             matrix, blockSize, onDisk, compressionLevel,
                             chunkSize, mDim1, mDim2, blocks, chrSwapped,
                             filesAreCool) {

    ## TODO ignore onDisk argument
    if (!onDisk) stop("Jagged arrays must be stored onDisk.")

    ## Create hdf5 for storage
    h5createFile(h5File)

    ## Generate & write slice data
    offset2 <- cumsum(mDim1*mDim2)
    offset1 <- c(0L, offset2[-length(offset2)])+1
    slices <- data.table(mDim1, mDim2, offset1, offset2)
    slices <- as.matrix(slices)
    h5write(slices, h5File, "slices")

    ## Storage for count data (column for each file)
    h5createDataset(file = h5File,
                    dataset = "counts",
                    dims = c(sum(mDim1*mDim2), length(files)),
                    chunk = c(1L, 1L),
                    storage.mode = "double",
                    fillValue = NA_real_,
                    level = compressionLevel)

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
          if(filesAreCool){
            sparseMat <-
              coolStraw(norm = norm,
                    fname = files[j],
                    chr1loc = blocks$chr1loc[i],
                    chr2loc = blocks$chr2loc[i],
                    binsize = binSize) |>
              as.data.table()
          } else {
            sparseMat <-
                straw(norm = norm,
                      fname = files[j],
                      chr1loc = blocks$chr1loc[i],
                      chr2loc = blocks$chr2loc[i],
                      unit = "BP",
                      binsize = binSize,
                      matrix = matrix) |>
                as.data.table()
          }

            ## Select interactions belonging to block
            xIndices <- blocks$xIndex[[i]]
            g <- as.data.table(x[xIndices])
            g$xIndices <- xIndices

            ## Create submatrix bins for each range
            ## with a fast cross-join
            longMat <- g[,{
                x <- seq(start1, end1-binSize, binSize);
                y <- seq(start2, end2-binSize, binSize);
                x <- as.integer(x) #seq sometimes returns double
                y <- as.integer(y)
                CJ(x, y, sorted=FALSE)
            },
            by=.(grp=xIndices)]

            ## Assign counts
            longMat <- .assignCounts(g, longMat,
                                     sparseMat,
                                     chrSwapped,
                                     half)

            ## Enumerate indices for writing
            ## from the slices
            slice <- slices[xIndices,,drop=FALSE]
            idx <- mapply(
                FUN=seq,
                from=slice[,"offset1"],
                to=slice[,"offset2"],
                by=1L,
                SIMPLIFY=FALSE
            )
            idx <- unlist(idx)

            ## Write data to hdf5 file
            h5write(obj=longMat$counts,
                    file=h5File,
                    name="counts",
                    index=list(idx, j))
        }
    }

    ## Close progress bar
    pb$tick(tokens = list(step = "Done!"))
    if(pb$finished) pb$terminate()

    ## Build colData and metadata
    colData <- DataFrame(files=files, fileNames=basename(files))
    metadata <- list(
        binSize=binSize,
        norm=norm,
        matrix=matrix
    )

    ## Create JaggedArray
    ja <- JaggedArray(
        h5File=h5File,
        dim=c(length(x), length(files)),
        subList=vector("list", 2L)
    )

    ## Create InteractionJaggedArray
    iset <- InteractionJaggedArray(
        interactions=x,
        colData=colData,
        counts=ja,
        metadata=metadata
    )
    return(iset)

}

#' Preprocessing before pulling Hi-C data
#' @inheritParams pullHicMatrices
#' @returns Updated `x` and `blocks` describing
#'  ranges to pull.
#' @noRd
.prepareInputs <- function(x, files, binSize, h5File, half, norm,
                           matrix, blockSize, onDisk, compressionLevel,
                           chunkSize) {

    ## Suppress NSE notes in R CMD check
    block <- xIndex <- NULL

    ## Check input types
    .checkTypes(c(
        binSize="number",
        half="string",
        h5File="string",
        norm="string",
        matrix="string",
        blockSize="number",
        onDisk="boolean",
        compressionLevel="number",
        chunkSize="number"
    ))
    
    ## Determine if files are .hic or .cool/.mcool
    filesAreCool <- sapply(files,
           function(x){
             ## If an error results from .checkIfCool, set to false
             isCool <- tryCatch(.checkIfCool(x),
                      error = function(e){
                        return(FALSE)
                        }) |>
               as.logical()
             
             ## No error gives NA, so set NAs to TRUE
             isCool[is.na(isCool)] <- TRUE
             
             return(isCool)
           })
    
    if (!(all(filesAreCool) | all(!filesAreCool))){
      abort(c("Mixed file types are not currently supported",
            "i" = "`files` must either be all `.hic` or all `.(m)cool`",
            "*" = glue("Split files by type and run this function",
                       " multiple times to proceed")))
    }
    
    filesAreCool <- all(filesAreCool)

    ## Parse straw parameters
    if(filesAreCool){
      .checkCoolArgs(files, half, norm, binSize, matrix)
    } else {
      .checkStrawArgs(files, half, norm, binSize, matrix)
    }

    ## Ensure seqnames are properly formatted
    if(filesAreCool){
      .checkCoolChroms(x, files, binSize)
    } else {
      .checkHicChroms(x, files)
    }

    ## Assign GInteractions to bins
    x <- .handleBinning(x, binSize)

    ## Order interactions so blocks are
    ## formatted for straw
    oi <- .orderInteractions(x, files[1])

    ## Assign x to blocks
    blocks <- as.data.table(.mapToBlocks(oi$x, blockSize))

    ## Paste ranges to get chr1 and chr2 locs
    blocks <-
        blocks[,.(chr1loc = paste(seqnames1, start1, end1, sep=":"),
                  chr2loc = paste(seqnames2, start2, end2, sep=":"),
                  block, xIndex)]

    ## Return x, blocks, and indices
    ## for swapped chroms & positions
    ## along with tag if files are cool or hic
    return(list(x=x, blocks=blocks,
                chrSwapped=oi$chrSwapped,
                filesAreCool=filesAreCool))
}

#' Internal pullHicMatrices
#' @inheritParams pullHicMatrices
#' @importFrom data.table as.data.table
#' @importFrom InteractionSet regions
#' @importFrom S4Vectors width
#' @importFrom rlang abort
#' @importFrom glue glue
#' @return Array of Hi-C submatrices.
#' @noRd
.pullHicMatrices <- function(x, files, binSize, h5File, half, norm,
                             matrix, blockSize, onDisk, compressionLevel,
                             chunkSize) {

    ## Prepare inputs for Hi-C processing
    dat <- .prepareInputs(x, files, binSize, h5File, half, norm,
                          matrix, blockSize, onDisk, compressionLevel,
                          chunkSize)

    ## Use region widths to dispatch code for
    ## extracting equal or variable dimension slices
    firstWidths <- unique(width(first(dat$x))) - 1
    secondWidths <- unique(width(second(dat$x))) - 1
    if (length(firstWidths) == 1L & length(secondWidths) == 1L) {

        ## Set matrix dimensions
        ## (round up to binSize if widths < binSize)
        mDim1 <- ceiling(firstWidths/binSize)
        mDim2 <- ceiling(secondWidths/binSize)

        ## Dispatch .pullArray for equal dimensions
        iset <- .pullArray(dat$x, files, binSize, h5File, half, norm,
                           matrix, blockSize, onDisk, compressionLevel,
                           chunkSize, mDim1, mDim2, dat$blocks,
                           dat$chrSwapped, dat$filesAreCool)

    } else {

        ## Set matrix dimensions
        mDim1 <- ceiling((width(first(dat$x)) - 1)/binSize)
        mDim2 <- ceiling((width(second(dat$x)) - 1)/binSize)

        ## Dispatch .pullJaggedArray for
        ## irregular dimensions
        iset <- .pullJaggedArray(dat$x, files, binSize, h5File, half, norm,
                                 matrix, blockSize, onDisk, compressionLevel,
                                 chunkSize, mDim1, mDim2, dat$blocks,
                                 dat$chrSwapped, dat$filesAreCool)
    }
    return(iset)
}

#' Pull submatrices from `.hic` files
#'
#' The dimensions of the pulled submatrix is
#' defined by dividing the widths of anchors in
#' `x` by the `binSize`. When the anchor widths
#' are the same for each interaction, an
#' InteractionArray is returned. However, if the
#' anchor widths differ in `x`, an
#' InteractionJaggedArray is returned instead.
#'
#' @param x GInteractions object containing interactions
#'  to extract from Hi-C files.
#' @param files Character file paths to `.hic` files.
#' @param binSize Integer (numeric) describing the
#'  resolution (range widths) of the paired data.
#' @param ... Additional arguments.
#' @param h5File Character file path to save `.h5` file.
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
#' @param norm String (length one character vector)
#'  describing the Hi-C normalization to apply. Use
#'  `strawr::readHicNormTypes()` to see accepted values
#'  for each `.hic` file in `files`. Use `readCoolNormTypes()`
#'  to see accepted values for each `.cool` or `.mcool` file
#'  in `files`.
#' @param matrix String (length one character vector)
#'  Type of matrix to extract. Must be one of "observed",
#'  "oe", or "expected". "observed" is observed counts,
#'  "oe" is observed/expected counts, "expected" is
#'  expected counts. For `.cool` and `.mcool` files, only
#'  extracting the "observed" matrix is supported.
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
#'  The submatrices returned have rows cooresponding to anchor1
#'  of `x` and columns correspond to anchor2 of `x`.
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
#' ## Read in loop pixels as GInteractions object
#' pixels <-
#'   WT_5kbLoops.txt() |>
#'   setNames("WT") |>
#'   read.table(header=TRUE) |>
#'   as_ginteractions(keep.extra.columns=FALSE) |>
#'   assignToBins(binSize=100e3)
#'
#' ## Removes the "chr" prefix for compatibility
#' ## with the preprocessed hic files
#' GenomeInfoDb::seqlevelsStyle(pixels) <- 'ENSEMBL'
#'
#' ## Expand pixels to regions for pulling
#' ## Hi-C submatrices
#' regions <- pixelsToMatrices(x=pixels, buffer=5)
#'
#' ## Extract 11x11 count matrices from the
#' ## first 100 regions and 2 Hi-C files
#' iarr <- pullHicMatrices(x=regions[1:100],
#'                         files=hicFiles,
#'                         binSize=100e3)
#' iarr
#'
#' ## Access count matrices
#' counts(iarr)
#'
#' ## Display the start bin of each
#' ## interaction in the count
#' ## matrices
#' counts(iarr, showDimnames=TRUE)
#'
#' ## InteractionJaggedArray example
#' gi <- read.table(text="
#'             1 51000000 51300000 1 51000000 51500000
#'             2 52000000 52300000 3 52000000 52500000
#'             1 150000000 150500000 1 150000000 150300000
#'             2 52000000 52300000 2 52000000 52800000") |>
#'     as_ginteractions()
#'
#' iarr <- pullHicMatrices(gi, hicFiles, 100e03, half="both")
#' iarr
#'
#' counts(iarr)
#'
#' @rdname pullHicMatrices
#' @export
setMethod("pullHicMatrices",
          signature(x='GInteractions',
                    files='character',
                    binSize='numeric'),
          definition=.pullHicMatrices)


#' Pull data from files and return as
#' InteractionSet with DelayedArrays
#'
#' @inheritParams pullHicMatrices
#' @param mDim integer - dimensions for matrices (equal)
#' @param blocks genomic blocks determining how much data to pull out at once
#' @param filesAreCool logical, TRUE if files are all .cool/.mcool
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
.pullMatrix <- function(x, files, binSize, h5File, half, norm, matrix,
                       blockSize, onDisk, compressionLevel,
                       chunkSize, blocks, chrSwapped, filesAreCool) {
    ## Determine dimensions for dataset
    ## Dim order is nInteractions, nFiles, matrix dims
    dims <- c(length(x), length(files))

    if (onDisk) {

        ## Overwrite existing file
        if (file.exists(h5File)) {
            file.remove(h5File)
        }

        ## Create hdf5 for storage
        h5createFile(h5File)

        ## Create dataset for holding array of counts
        ## and row/colnames
        h5createDataset(file = h5File,
                        dataset = "counts",
                        dims = dims,
                        chunk = c(chunkSize, 1),
                        storage.mode = "double",
                        fillValue = NA_real_,
                        level = compressionLevel)

    } else {
        counts <- array(data=NA, dim = dims)
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
          if(filesAreCool){
            sparseMat <-
              coolStraw(norm = norm,
                    fname = files[j],
                    chr1loc = blocks$chr1loc[i],
                    chr2loc = blocks$chr2loc[i],
                    binsize = binSize) |>
              as.data.table()
          } else {
            sparseMat <-
                straw(norm = norm,
                      fname = files[j],
                      chr1loc = blocks$chr1loc[i],
                      chr2loc = blocks$chr2loc[i],
                      unit = "BP",
                      binsize = binSize,
                      matrix = matrix) |>
                as.data.table()
          }

            ## Select interactions belonging to block
            xIndices <- blocks$xIndex[[i]]
            g <- as.data.table(x[xIndices])
            g$xIndices <- xIndices

            ## Create submatrix bins for each range
            ## with a fast cross-join (rename)
            longMat <- g[,{
                x <- start1;
                y <- start2;
                CJ(x, y, sorted=FALSE)
            },
            by=.(grp=xIndices)]

            ## Assign counts
            longMat <- .assignCounts(g, longMat,
                                     sparseMat,
                                     chrSwapped,
                                     half)

            ## Fill array by col, row, interactions, file
            ## then rearrange for storage: interactions, file, rows, cols
            a <- array(data=longMat$counts,
                       dim=c(length(xIndices), 1))

            if (onDisk) {
                # Write data to hdf5 file
                h5write(a, h5File, "counts", index = list(xIndices, j))
            } else {
                counts[xIndices, j] <- a
            }

        }
    }

    ## Close progress bar
    pb$tick(tokens = list(step = "Done!"))
    if(pb$finished) pb$terminate()

    ## Build colData and metadata
    colData <- DataFrame(files=files, fileNames=basename(files))
    metadata <- list(
        binSize=binSize,
        norm=norm,
        matrix=matrix
    )

    ## Construct InteractionMatrix object
    if (onDisk) {
        iset <-
            InteractionMatrix(
                interactions = x,
                assays = list(
                    counts = DelayedArray(HDF5Array(h5File,"counts"))
                ),
                colData = colData,
                metadata = metadata
            )
    } else {
        iset <-
            InteractionMatrix(
                interactions = x,
                assays = list(
                    counts = DelayedArray(counts)
                ),
                colData = colData,
                metadata = metadata
            )
    }

    ## Assign dimnames to files (default is basename(files))
    if (is.null(names(files))) {
        dimnames(iset)[2] <- list(iset$fileNames)
    } else {
        dimnames(iset)[2] <- list(names(files))
    }

    return(iset)
}


#' Internal pullHicPixels
#' @inheritParams pullHicPixels
#' @importFrom glue glue
#' @return Matrix of Hi-C submatrices.
#' @noRd
.pullHicPixels <- function(x, files, binSize, h5File, half, norm, matrix,
                           blockSize, onDisk, compressionLevel,
                           chunkSize) {

    ## Prepare inputs for Hi-C processing
    dat <- .prepareInputs(x, files, binSize, h5File, half, norm, matrix,
                          blockSize, onDisk, compressionLevel,
                          chunkSize)

    ## Ensure widths are equal to binSize
    if (!.checkBinnedPairs(x, binSize)) {
        abort(c("`x` are not binned to `binSize`.",
                'i' = glue("All ranges in `x` must be equal \\
                           widths for `pullHicPixels()`."),
                'i' = glue("Use `assignToBins()` to bin or use \\
                           `pullHicMatrices()` instead.")))
    }

    ## Dispatch pulling pixels
    iset <- .pullMatrix(dat$x, files, binSize, h5File, half,
                        norm, matrix, blockSize, onDisk,
                        compressionLevel, chunkSize,
                        dat$blocks, dat$chrSwapped, dat$filesAreCool)

    return(iset)

}


#' Pull contact frequency from `.hic` files
#'
#' @param x GInteractions object containing interactions
#'  to extract from Hi-C files.
#' @param files Character file paths to `.hic` files.
#' @param binSize Integer (numeric) describing the
#'  resolution (range widths) of the paired data.
#' @param ... Additional arguments.
#' @param h5File Character file path to save `.h5` file.
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
#' @returns InteractionSet object with a 2-dimensional array
#'  of Hi-C interactions (rows) and Hi-C sample (columns).
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
#' ## Read in loop pixels as GInteractions object
#' pixels <-
#'   WT_5kbLoops.txt() |>
#'   setNames("WT") |>
#'   read.table(header=TRUE) |>
#'   as_ginteractions(keep.extra.columns=FALSE) |>
#'   assignToBins(binSize=100e3)
#'
#' ## Removes the "chr" prefix for compatibility
#' ## with the preprocessed hic files
#' GenomeInfoDb::seqlevelsStyle(pixels) <- 'ENSEMBL'
#'
#' ## Extract the first 100 pixels
#' imat <- pullHicPixels(x=pixels[1:100],
#'                       files=hicFiles,
#'                       binSize=100e3)
#' imat
#'
#' ## Access count matrix
#' counts(imat)
#'
#' @rdname pullHicPixels
#' @export
setMethod("pullHicPixels",
          signature(x='GInteractions',
                    files='character',
                    binSize='numeric'),
          definition=.pullHicPixels)
