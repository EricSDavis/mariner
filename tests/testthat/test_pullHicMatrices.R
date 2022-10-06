library(mariner)
library(data.table)
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
# giList <-
#     lapply(giList, \(x) {
#         ## Create a sample of regions to permute
#         set.seed(123)
#         s <- sample(x = seq_len(length(regions(x))),
#                     size = length(regions(x)),
#                     replace = FALSE)
#
#         ## Scramble regions
#         regions(x) <- regions(x)[s]
#         x
#     })

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

    .checkStrawArgs(files = hicFiles,
                    norm = "KR",
                    binSize = 10e03,
                    matrix = "oe") |>
        expect_null()

    .checkStrawArgs(files = hicFiles,
                    norm = "SCALE",
                    binSize = 10e03,
                    matrix = "oe") |>
        expect_error()

    .checkStrawArgs(files = hicFiles,
                    norm = "none",
                    binSize = 10e03,
                    matrix = "oe") |>
        expect_error()

    .checkStrawArgs(files = hicFiles,
                    norm = "KR",
                    binSize = 101e03,
                    matrix = "oe") |>
        expect_error()

    .checkStrawArgs(files = hicFiles,
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
    # files = hicFiles
    # norm = "NONE"
    # matrix = "observed"
    # onDisk = FALSE
    # compressionLevel = 0
    # # chunkSize = length(x)
    # rm(chunkSize)


    iset <-
        pullHicMatrices(x = x,
                        binSize = 2.5e06,
                        files = hicFiles[1],
                        onDisk=TRUE)
    hictoolsr::calcApa(x |> as.data.table() |> as_ginteractions(),
                       hicFiles[1], res = 2.5e06, buffer = 2)

    assay(iset) |>
        aperm(c(3,4, 1,2)) |>
        {\(x) apply(x[,,,1], c(1,2), sum)}()

    pullHicMatrices(x = x[1:10],
                    binSize = 2.5e06,
                    files = hicFiles)

    .pullHicMatrices(x = x[1:10],
                     binSize = 2.5e06,
                     files = hicFiles,
                     norm = "KR",
                     matrix = "observed",
                     blockSize = 248956422,
                     onDisk = TRUE,
                     compressionLevel = 0,
                     chunkSize = 1)

    interactions(iset)
    aperm(assay(iset), c(3,4,1,2)) |>
        {\(x) apply(x[,,,1], c(1,2), sum)}()

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

    # hictoolsr::calcApa(x, files[1], res = 250e03, buffer = 50)

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
