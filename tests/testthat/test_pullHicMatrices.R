library(mariner)
library(data.table, include.only = "fread")
library(InteractionSet)
library(glue, include.only = "glue")
library(strawr)
library(dplyr, include.only = "arrange")
library(GenomeInfoDb)
library(GenomicRanges)

## Shared objects --------------------------------------------------------------

## Test .hic files
hicFiles <-
    system.file("extdata/test_hic", package = "mariner") |>
    list.files(pattern = ".hic", full.names = TRUE)

## Reference BEDPE files (loops called with SIP)
bedpeFiles <-
    system.file("extdata", package = "mariner") |>
    list.files(pattern = "Loops.txt", full.names = TRUE)

## Read in bedpeFiles as a list of GInteractions
giList <-
    lapply(bedpeFiles, fread) |>
    lapply(as_ginteractions)

## Create out-of-ordered, interchromosomal ranges
## (for testing)
giList <-
    lapply(giList, \(x) {
        ## Create a sample of regions to permute
        set.seed(123)
        s <- sample(x = seq_len(length(regions(x))),
                    size = length(regions(x)),
                    replace = FALSE)

        ## Scramble regions
        regions(x) <- regions(x)[s]
        x
    })

## Merge pairs
mgi <- mergePairs(x = giList,
                  radius = 50e03)

## Bin MergedGInteractions
# bgi <- binPairs(x = mgi, binSize = 50e03)
bgi <- binPairs(x = mgi, binSize = 250e03)

## Bin with snapping
sgi <- snapToBins(x = mgi, binSize = 50e03)

## Test pullHicMatrices --------------------------------------------------------

## TODO: Change .pullHicMatrices to pullHicMatrices
test_that("pullHicMatrices checks for binning", {

    ## Inform if not binned
    .handleBinning(x = mgi, binSize = 50e03) |>
        expect_message("Pairs are not binned.*")

    ## Still returns binned object
    .handleBinning(x = mgi, binSize = 50e03) |>
        suppressMessages() |>
        expect_s4_class("MergedGInteractions")

    ## Accepts binned pairs
    .handleBinning(x = bgi,
                   binSize = 50e03) |>
        expect_s4_class("MergedGInteractions")

    ## Throws message when using a non-divisible binSize
    .handleBinning(x = bgi,
                   binSize = 250e03) |>
        expect_message(".*Snapping to binSize=250000.*")

    ## No message for divisible binSizes
    .handleBinning(x = sgi,
                   binSize = 10e03) |>
        expect_message(NA)

})

test_that("straw returns data in expected order", {

    ## Define chromosome map & arrange by index
    chrMap <-
        readHicChroms(hicFiles[1]) |>
        arrange(index) |>
        head(n=25)

    ## Check every combination of chromosomes
    for (i in seq(2, nrow(chrMap))) {
        for (j in seq(2, nrow(chrMap))) {
            chr1loc <- glue('{chrMap[i, "name"]}:0:0')
            chr2loc <- chrMap[j, "name"]
            sparseMat <-
                strawr::straw(norm = "NONE",
                              fname = hicFiles[1],
                              chr1loc = chr1loc,
                              chr2loc = chr2loc,
                              unit = "BP",
                              binsize = 2500000,
                              matrix = "observed")

            ## All values in column "x" should be 0
            if (nrow(sparseMat) != 0) {
                if (i <= j) {
                    expect_true(all(sparseMat$x == 0),
                                label = paste(i,j))
                } else {
                    expect_false(all(sparseMat$x == 0),
                                 label = paste(i,j))
                }
            }
        }
    }
})

test_that("Check chromosomes in .hic file", {

    ## Assign to x (to avoid modifying in place)
    x <- bgi

    ## Seqnames mismatch throws error
    .checkHicChroms(x, hicFiles[1]) |>
        expect_error("There's.*craziness.*")

    ## Corrected seqnames don't throw error
    seqlevelsStyle(x) <- "ENSEMBL"
    .checkHicChroms(x, hicFiles[1]) |>
        expect_null()

})

test_that("Blocking works for individual ranges", {

    .blockRanges(start = c(0, 40, 45, 50, 110, 120, 0),
                 end = c(10, 50, 55, 60, 120, 130, 25),
                 blockSize = 50) |>
        expect_identical(
            list(
                start = c(0, 25, 25, 50, 100, 100, 0),
                end = c(50, 75, 75, 100, 150, 150, 50)
            )
        )

    .blockRanges(start = 0, end = 100, blockSize = 50) |>
        expect_error(".*blockSize/2.*")

    .blockRanges(start = 0, end = 50, blockSize = 50) |>
        expect_error(".*blockSize/2.*")

    .blockRanges(start = 0, end = 26, blockSize = 50) |>
        expect_error(".*blockSize/2.*")

})

test_that("GInteractions correctly map to blocks", {

    x <-
        GInteractions(
            GRanges(seqnames = "chr1",
                    ranges = IRanges(start = c(0, 30, 10, 0, 0),
                                     end = c(20, 50, 30, 20, 20))),
            GRanges(seqnames = c(rep("chr1", 4), "chr2"),
                    ranges = IRanges(start = c(0, 30, 40, 10, 10),
                                     end = c(20, 50, 60, 30, 30)))
        )

    exp <-
        GInteractions(
            GRanges(seqnames = "chr1",
                    ranges = IRanges(start = c(0, 25, 0, 0),
                                     end = c(50, 75, 50, 50))),
            GRanges(seqnames = c(rep("chr1", 3), "chr2"),
                    ranges = IRanges(start = c(0, 25, 25, 0),
                                     end = c(50, 75, 75, 50)))
        )
    exp$block <- seq(1L,4L)
    exp$xIndex <- list(c(1L,4L), 2L, 3L, 5L)

    .mapToBlocks(x, 50) |>
        expect_identical(exp)

})

test_that("Straw args are checked correctly", {

    .checkStrawArgs(file = hicFiles,
                    norm = "KR",
                    binSize = 10e03,
                    matrix = "oe") |>
        expect_null()

    .checkStrawArgs(file = hicFiles,
                    norm = "SCALE",
                    binSize = 10e03,
                    matrix = "oe") |>
        expect_error()

    .checkStrawArgs(file = hicFiles,
                    norm = "none",
                    binSize = 10e03,
                    matrix = "oe") |>
        expect_error()

    .checkStrawArgs(file = hicFiles,
                    norm = "KR",
                    binSize = 101e03,
                    matrix = "oe") |>
        expect_error()

    .checkStrawArgs(file = hicFiles,
                    norm = "KR",
                    binSize = 10e03,
                    matrix = "obs") |>
        expect_error()
})

test_that("develop the function (matrix cases)", {

    ## Assign to x (to avoid modifying in place)
    x <- bgi
    seqlevelsStyle(x) <- "ENSEMBL"
    # binSize = 5e03#1e06
    # blockSize = 248956422#100e06
    # file = hicFiles
    # norm = "NONE"
    # matrix = "observed"
    # onDisk = FALSE
    # compressionLevel = 0
    # # chunkSize = length(x)
    # rm(chunkSize)

    iset <-
        .pullHicMatrices(x = resize(x, width = 2.5e06*3, fix = 'center'),
                         binSize = 2.5e06,
                         file = hicFiles,
                         norm = "KR",
                         matrix = "observed",
                         blockSize = 248956422,
                         onDisk = TRUE,
                         compressionLevel = 0,
                         chunkSize = 1)

    interactions(iset)
    aperm(assay(iset), c(3,4,1,2)) |>
        {\(x) apply(x[,,,1], c(1,2), sum)}()

    ## Check for length-1  arguments
    .checkVectorLengths(list(
        binSize=binSize,
        blockSize=blockSize,
        norm=norm,
        matrix=matrix,
        onDisk=onDisk,
        compressionLevel=compressionLevel
    ))

    ## Parse straw parameters
    .checkStrawArgs(file, norm, binSize, matrix)

    ## Ensure seqnames are properly formatted
    .checkHicChroms(x, file)

    ## Assign GInteractions to bins
    x <- .handleBinning(x, binSize)

    ## Order interactions according to chroms
    ## (ensures blocks will be ordered too)
    x <- .orderInteractions(x, file[1])

    ## Assign x to blocks
    blocks <- as.data.table(.mapToBlocks(x, blockSize))
    blockMap <- blocks[,c("block", "xIndex")]

    ## Paste ranges to get chr1 and chr2 locs
    blocks <-
        blocks[,.(chr1loc = paste(seqnames1, start1, end1, sep=":"),
                  chr2loc = paste(seqnames2, start2, end2, sep=":"))]

    ## TODO: Code will differ here depending on input cases


    ## Use region widths to dispatch code for
    ## extracting equal or variable dimension slices
    widths <- unique(width(regions(x))) - 1
    if (length(widths) == 1L) {

        ## Set matrix dimensions
        mDim <- (widths/binSize) + 1

        ## TODO: Dispatch pullHicMatrices (equal dimensions)

    } else {
        ## TODO: Dispatch pullHicMatrices (varible dimensions)
    }

    ## Determine dimensions for dataset
    ## Dim order is nInteractions, nFiles, matrix dims
    dims <- c(length(x), length(file), mDim, mDim)

    if (onDisk) {
        ## Create hdf5 for storage
        library(rhdf5)
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
    pb <- progress::progress_bar$new(
        format = "  :step [:bar] :percent elapsed: :elapsedfull",
        clear = F, total = nrow(blocks)*length(file)+1)
    pb$tick(0)

    for(j in seq_len(length(file))) {
        for(i in seq_len(nrow(blocks))) {

            ## Update progress
            pb$tick(tokens=list(step=sprintf(
                'Pulling block %s of %s from file %s of %s',
                i, nrow(blocks), j, length(file)
            )))

            ## Extract block data from file
            sparseMat <-
                straw(norm = norm,
                      fname = file[j],
                      chr1loc = blocks$chr1loc[i],
                      chr2loc = blocks$chr2loc[i],
                      unit = "BP",
                      binsize = binSize,
                      matrix = matrix) |>
                as.data.table()

            ## Select interactions belonging to block
            xIndices <- blockMap$xIndex[[i]]
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
    ## END LOOP

    ## Close progress bar
    pb$tick(tokens = list(step = "Done!"))
    if(pb$finished) pb$terminate()

    ## Construct InteractionSet object
    if (onDisk) {
        library(DelayedArray)
        library(HDF5Array)
        iset <-
            InteractionSet::InteractionSet(
                interactions = x,
                assays = list(
                    counts = DelayedArray(HDF5Array(h5,"counts")),
                    rownames = DelayedArray(HDF5Array(h5,"rownames")),
                    colnames = DelayedArray(HDF5Array(h5,"colnames"))
                )
            )

    } else {
        iset <-
            InteractionSet::InteractionSet(
                interactions = x,
                assays = list(
                    counts = DelayedArray(counts),
                    rownames = DelayedArray(rownames),
                    colnames = DelayedArray(colnames)
                )
            )
    }


    assay(iset) <- as(assay(iset), "HDF5Array")
    tmp
    path(assay(iset, "counts"))

    iset
    ## TODO: make custom print method for arrays
    library(SummarizedExperiment)
    tmp <-
        assay(iset[1:3]) |>
        aperm(c(3,4,1,2))

    tmp

    # hictoolsr::calcApa(x, file[1], res = 250e03, buffer = 50)

    ## Get result out
    # tmp <-
    #     h5read(h5, "counts", index = list(xIndices, 1, 1:2, 1:2)) |>
    #     aperm(c(3,4,1,2)) |>
    #     {\(x) x[,,1,1]}()
    tmp <-
        h5read(h5, "counts", index = list(xIndices[9:10], 1, 1:2, 1:2)) |>
        aperm(c(3,4,1,2)) # rearrange to form matrix
    rownames(tmp) <-
        h5read(h5,
               "rownames",
               index=list(xIndices[9:10],1,1:2))[,,1:2]
    colnames(tmp) <-
        h5read(h5,
               "colnames",
               index=list(xIndices[9:10],1,1:2))[,,1:2]



    ## Example select one matrix
    # b <- a[,,10,1]
    # rownames(b) <- unlist(aInfo$rowNames[10])
    # colnames(b) <- unlist(aInfo$colNames[10])
    # b


    ## Slower way and doesn't store correct row/colnames
    # system.time({
    #     tmp <-
    #         lapply(unique(longMat$groupRow), \(i){
    #             reshape2::acast(longMat[groupRow==i], x~y, value.var = "counts")
    #         }) |>
    #         simplify2array()
    # })


})
