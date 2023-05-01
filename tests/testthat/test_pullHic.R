library(mariner)
library(marinerData)
library(data.table)
library(InteractionSet)
library(glue, include.only = "glue")
library(strawr)
library(dplyr, include.only = "arrange")
library(GenomeInfoDb)
library(GenomicRanges)

## Shared objects --------------------------------------------------------------

## Read .hic file paths
hicFiles <- c(
    LEUK_HEK_PJA27_inter_30.hic(),
    LEUK_HEK_PJA30_inter_30.hic()
)
names(hicFiles) <- c("FS", "WT")

## Reference BEDPE files (loops called with SIP)
bedpeFiles <- c(
    FS_5kbLoops.txt(),
    WT_5kbLoops.txt()
)
names(bedpeFiles) <- c("FS", "WT")

## Read in bedpeFiles as a list of GInteractions
giList <-
    lapply(bedpeFiles, fread) |>
    lapply(as_ginteractions)

## Create out-of-order, interchromosomal ranges
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
bgi <- binPairs(x = mgi, binSize = 50e03)

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
                    half="both",
                    norm = "KR",
                    binSize = 10e03,
                    matrix = "oe") |>
        expect_null()

    .checkStrawArgs(files = hicFiles,
                    half="neither",
                    norm = "KR",
                    binSize = 10e03,
                    matrix = "oe") |>
        expect_error()

    .checkStrawArgs(files = hicFiles,
                    half="both",
                    norm = "SCALE",
                    binSize = 10e03,
                    matrix = "oe") |>
        expect_error()

    .checkStrawArgs(files = hicFiles,
                    half="both",
                    norm = "none",
                    binSize = 10e03,
                    matrix = "oe") |>
        expect_error()

    .checkStrawArgs(files = hicFiles,
                    half="both",
                    norm = "KR",
                    binSize = 101e03,
                    matrix = "oe") |>
        expect_error()

    .checkStrawArgs(files = hicFiles,
                    half="both",
                    norm = "KR",
                    binSize = 10e03,
                    matrix = "obs") |>
        expect_error()
})

test_that("pullHicPixels pulls correct counts", {

    ## Assign to x (to avoid modifying in place)
    x <- binPairs(bgi[1:10], 2.5e06)
    seqlevelsStyle(x) <- "ENSEMBL"

    ## Give names to hicFiles
    namedHicFiles <- hicFiles
    names(namedHicFiles) <- c("Mut", "WT")

    ## Pull pixels
    imat <-
        pullHicPixels(x=x,
                      binSize=2.5e06,
                      files=namedHicFiles,
                      matrix="observed",
                      norm="NONE",
                      half="both")

    ## Get reference counts with straw directly
    chr1loc <- paste(seqnames1(x), start1(x), start1(x), sep=":")
    chr2loc <- paste(seqnames2(x), start2(x), start2(x), sep=":")
    pullCounts <- \(fname,i) {
        straw(norm="NONE",
              fname=fname,
              chr1loc=chr1loc[i],
              chr2loc=chr2loc[i],
              unit="BP",
              binsize=2.5e06,
              matrix="observed")$counts
    }
    Mut <-
        lapply(seq_along(chr1loc), \(i) {
            cnts <- pullCounts(namedHicFiles["Mut"],i)
            if (length(cnts)==0) cnts <- 0
            cnts
        }) |>
        do.call(rbind, args=_)

    WT <-
        lapply(seq_along(chr1loc), \(i) {
            cnts <- pullCounts(namedHicFiles["WT"],i)
            if (length(cnts)==0) cnts <- 0
            cnts
        }) |>
        do.call(rbind, args=_)

    ## Test
    expect_identical(as.matrix(counts(imat)),
                     as.matrix(data.frame(Mut, WT)))
})

test_that("pullHicMatrices for square, regular arrays", {

    ## Intrachromosomal
    ## On-diagonal square
    gi <- read.table(text="
            1 51500000 52000000 1 51500000 52000000
            1 150000000 150500000 1 150000000 150500000") |>
        as_ginteractions()
    res <- pullHicMatrices(gi, hicFiles[1], 100e03, half="upper")

    ## half="both"
    res <- pullHicMatrices(gi, hicFiles[1], 100e03, half="both")
    exp <- array(data=c(67, 22, 13, 9, 8, 22, 59, 20, 10, 7,
                        13, 20, 51, 15, 6, 9, 10, 15, 49, 23,
                        8, 7, 6, 23, 66, 63, 25, 15, 4, 2, 25,
                        68, 28, 7, 3, 15, 28, 87, 45, 11, 4, 7,
                        45, 100, 26, 2, 3, 11, 26, 86),
                 dim=c(5,5,2,1))
    expect_identical(interactions(res), gi)
    expect_identical(unname(as.array(counts(res))), exp)

    ## half="upper"
    res <- pullHicMatrices(gi, hicFiles[1], 100e03, half="upper")
    exp <- array(data=c(67, NA, NA, NA, NA, 22, 59, NA, NA, NA,
                        13, 20, 51, NA, NA, 9, 10, 15, 49, NA,
                        8, 7, 6, 23, 66, 63, NA, NA, NA, NA, 25,
                        68, NA, NA, NA, 15, 28, 87, NA, NA, 4, 7,
                        45, 100, NA, 2, 3, 11, 26, 86),
                 dim=c(5,5,2,1))
    expect_identical(interactions(res), gi)
    expect_identical(unname(as.array(counts(res))), exp)

    ## half="lower"
    res <- pullHicMatrices(gi, hicFiles[1], 100e03, half="lower")
    exp <- array(data=c(67, 22, 13, 9, 8, NA, 59, 20, 10, 7, NA,
                        NA, 51, 15, 6, NA, NA, NA, 49, 23, NA,
                        NA, NA, NA, 66, 63, 25, 15, 4, 2, NA, 68,
                        28, 7, 3, NA, NA, 87, 45, 11, NA, NA, NA,
                        100, 26, NA, NA, NA, NA, 86),
                 dim=c(5,5,2,1))
    expect_identical(interactions(res), gi)
    expect_identical(unname(as.array(counts(res))), exp)


    ## Intrachromosomal
    ## Off-diagonal
    ## Define square where start1 < start2
    ## (i.e. upper triangular)
    gi <- read.table(text="
            1 51500000 51700000 1 51800000 52000000
            1 150000000 150200000 1 150300000 150500000") |>
        as_ginteractions()
    expect_true(all(start1(gi) < start2(gi)))

    ## half="both"
    res <- pullHicMatrices(gi, hicFiles[1], 100e03, half="both")
    exp <- array(data=c(9, 10, 8, 7, 4, 7, 2, 3),
                 dim=c(2,2,2,1))
    expect_identical(interactions(res), gi)
    expect_identical(unname(as.array(counts(res))), exp)

    ## half="upper"
    res <- pullHicMatrices(gi, hicFiles[1], 100e03, half="upper")
    exp <- array(data=c(9, 10, 8, 7, 4, 7, 2, 3),
                 dim=c(2,2,2,1))
    expect_identical(interactions(res), gi)
    expect_identical(unname(as.array(counts(res))), exp)

    ## half="lower"
    res <- pullHicMatrices(gi, hicFiles[1], 100e03, half="lower")
    exp <- array(data=NA_real_, dim=c(2,2,2,1))
    expect_identical(interactions(res), gi)
    expect_identical(unname(as.array(counts(res))), exp)


    ## Intrachromosomal
    ## Off-diagonal
    ## Define square where start1 > start2
    ## (i.e. lower triangular)
    gi <- read.table(text="
            1 51800000 52000000 1 51500000 51700000
            1 150300000 150500000 1 150000000 150200000") |>
        as_ginteractions()
    expect_true(all(start1(gi) > start2(gi)))

    ## half="both"
    res <- pullHicMatrices(gi, hicFiles[1], 100e03, half="both")
    exp <- array(data=c(9, 8, 10, 7, 4, 2, 7, 3),
                 dim=c(2,2,2,1))
    expect_identical(interactions(res), gi)
    expect_identical(unname(as.array(counts(res))), exp)

    ## half="upper"
    res <- pullHicMatrices(gi, hicFiles[1], 100e03, half="upper")
    exp <- array(data=NA_real_, dim=c(2,2,2,1))
    expect_identical(interactions(res), gi)
    expect_identical(unname(as.array(counts(res))), exp)

    ## half="lower"
    res <- pullHicMatrices(gi, hicFiles[1], 100e03, half="lower")
    exp <- array(data=c(9, 8, 10, 7, 4, 2, 7, 3),
                 dim=c(2,2,2,1))
    expect_identical(interactions(res), gi)
    expect_identical(unname(as.array(counts(res))), exp)


    ## Interchromosomal
    ## Define square where seqnames1 < seqnames2
    gi <- read.table(text="
            1 50000000 55000000 2 55000000 60000000
            1 150000000 155000000 2 155000000 160000000") |>
        as_ginteractions() |>
        suppressWarnings()
    expect_true(all(seqnames1(gi) < seqnames2(gi)))

    ## half="both"
    res <- pullHicMatrices(gi, hicFiles[1], 2.5e06, half="both")
    exp <- array(data=c(1, 1, 1, 0, 0, 1, 0, 0),
                 dim=c(2,2,2,1))
    expect_identical(interactions(res), gi)
    expect_identical(unname(as.array(counts(res))), exp)

    ## half="upper"
    res <- pullHicMatrices(gi, hicFiles[1], 2.5e06, half="upper")
    expect_identical(interactions(res), gi)
    expect_identical(unname(as.array(counts(res))), exp)

    ## half="lower"
    res <- pullHicMatrices(gi, hicFiles[1], 2.5e06, half="lower")
    expect_identical(interactions(res), gi)
    expect_identical(unname(as.array(counts(res))), exp)

    ## Interchromosomal
    ## Define square where seqnames1 > seqnames2
    gi <- read.table(text="
            2 50000000 55000000 1 55000000 60000000
            2 150000000 155000000 1 155000000 160000000") |>
        as_ginteractions() |>
        suppressWarnings()
    expect_true(all(seqnames1(gi) > seqnames2(gi)))

    ## half="both"
    res <- pullHicMatrices(gi, hicFiles[1], 2.5e06, half="both")
    exp <- array(data=c(0, 0, 1, 1, 2, 1, 2, 0),
                 dim=c(2,2,2,1))
    expect_identical(interactions(res), gi)
    expect_identical(unname(as.array(counts(res))), exp)

    ## half="upper"
    res <- pullHicMatrices(gi, hicFiles[1], 2.5e06, half="upper")
    expect_identical(interactions(res), gi)
    expect_identical(unname(as.array(counts(res))), exp)

    ## half="lower"
    res <- pullHicMatrices(gi, hicFiles[1], 2.5e06, half="lower")
    expect_identical(interactions(res), gi)
    expect_identical(unname(as.array(counts(res))), exp)

})

test_that("pullHicMatrices for rectangular, regular arrays", {

    ## Intrachromosomal
    ## On-diagonal rectangle
    gi <- read.table(text="
            1 51000000 51300000 1 51000000 51500000
            1 150000000 150300000 1 150000000 150500000") |>
        as_ginteractions()

    ## half="both"
    res <- pullHicMatrices(gi, hicFiles[1], 100e03, half="both")
    exp <- array(data=c(53, 15, 5, 15, 68, 19, 5, 19, 69, 1,
                        8, 12, 4, 5, 2, 63, 25, 15, 25, 68, 28,
                        15, 28, 87, 4, 7, 45, 2, 3, 11),
                 dim=c(3,5,2,1))
    expect_identical(interactions(res), gi)
    expect_identical(unname(as.array(counts(res))), exp)

    ## half="upper"
    res <- pullHicMatrices(gi, hicFiles[1], 100e03, half="upper")
    exp <- array(data=c(53, NA, NA, 15, 68, NA, 5, 19, 69, 1,
                        8, 12, 4, 5, 2, 63, NA, NA, 25, 68, NA,
                        15, 28, 87, 4, 7, 45, 2, 3, 11),
                 dim=c(3,5,2,1))
    expect_identical(interactions(res), gi)
    expect_identical(unname(as.array(counts(res))), exp)

    ## half="lower"
    res <- pullHicMatrices(gi, hicFiles[1], 100e03, half="lower")
    exp <- array(data=c(53, 15, 5, NA, 68, 19, NA, NA, 69, NA,
                        NA, NA, NA, NA, NA, 63, 25, 15, NA, 68,
                        28, NA, NA, 87, NA, NA, NA, NA, NA, NA),
                 dim=c(3,5,2,1))
    expect_identical(interactions(res), gi)
    expect_identical(unname(as.array(counts(res))), exp)


    ## Intrachromosomal
    ## Off-diagonal
    ## Define rectangle where start1 < start2
    ## (i.e. upper triangular)
    gi <- read.table(text="
            1 51000000 51300000 1 51300000 51800000
            1 150000000 150300000 1 150300000 150800000") |>
        as_ginteractions()
    expect_true(all(start1(gi) < start2(gi)))

    ## half="both"
    res <- pullHicMatrices(gi, hicFiles[1], 100e03, half="both")
    exp <- array(data=c(1, 8, 12, 4, 5, 2, 1, 5, 2, 2, 3, 8,
                        1, 4, 2, 4, 7, 45, 2, 3, 11, 6, 7, 5,
                        3, 3, 3, 6, 2, 3),
                 dim=c(3,5,2,1))
    expect_identical(interactions(res), gi)
    expect_identical(unname(as.array(counts(res))), exp)

    ## half="upper"
    res <- pullHicMatrices(gi, hicFiles[1], 100e03, half="upper")
    exp <- array(data=c(1, 8, 12, 4, 5, 2, 1, 5, 2, 2, 3, 8,
                        1, 4, 2, 4, 7, 45, 2, 3, 11, 6, 7, 5,
                        3, 3, 3, 6, 2, 3),
                 dim=c(3,5,2,1))
    expect_identical(interactions(res), gi)
    expect_identical(unname(as.array(counts(res))), exp)

    ## half="lower"
    res <- pullHicMatrices(gi, hicFiles[1], 100e03, half="lower")
    exp <- array(data=NA_real_, dim=c(3,5,2,1))
    expect_identical(interactions(res), gi)
    expect_identical(unname(as.array(counts(res))), exp)


    ## Intrachromosomal
    ## Off-diagonal
    ## Define rectangle where start1 > start2
    ## (i.e. lower triangular)
    gi <- read.table(text="
            1 51300000 51800000 1 51000000 51300000
            1 150300000 150800000 1 150000000 150300000") |>
        as_ginteractions()
    expect_true(all(start1(gi) > start2(gi)))

    ## half="both"
    res <- pullHicMatrices(gi, hicFiles[1], 100e03, half="both")
    exp <- array(data=c(1, 4, 1, 2, 1, 8, 5, 5, 3, 4, 12, 2, 2,
                        8, 2, 4, 2, 6, 3, 6, 7, 3, 7, 3, 2, 45,
                        11, 5, 3, 3),
                 dim=c(5,3,2,1))
    expect_identical(interactions(res), gi)
    expect_identical(unname(as.array(counts(res))), exp)

    ## half="upper"
    res <- pullHicMatrices(gi, hicFiles[1], 100e03, half="upper")
    exp <- array(data=NA_real_, dim=c(5,3,2,1))
    expect_identical(interactions(res), gi)
    expect_identical(unname(as.array(counts(res))), exp)

    ## half="lower"
    res <- pullHicMatrices(gi, hicFiles[1], 100e03, half="lower")
    exp <- array(data=c(1, 4, 1, 2, 1, 8, 5, 5, 3, 4, 12, 2, 2,
                        8, 2, 4, 2, 6, 3, 6, 7, 3, 7, 3, 2, 45,
                        11, 5, 3, 3),
                 dim=c(5,3,2,1))
    expect_identical(interactions(res), gi)
    expect_identical(unname(as.array(counts(res))), exp)


    ## Interchromosomal
    ## Pull rectangle where seqnames1 < seqnames2
    gi <- read.table(text="
            1 50000000 57500000 2 50000000 62500000
            1 150000000 157500000 2 150000000 162500000") |>
        as_ginteractions() |>
        suppressWarnings()
    expect_true(all(seqnames1(gi) < seqnames2(gi)))

    ## half="both"
    res <- pullHicMatrices(gi, hicFiles[1], 2.5e06, half="both")
    exp <- array(data=c(1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1,
                        0, 3, 2, 0, 2, 0, 2, 1, 0, 1, 0, 0, 0,
                        2, 1, 1, 0),
                 dim=c(3,5,2,1))
    expect_identical(interactions(res), gi)
    expect_identical(unname(as.array(counts(res))), exp)

    ## half="upper"
    res <- pullHicMatrices(gi, hicFiles[1], 2.5e06, half="upper")
    expect_identical(interactions(res), gi)
    expect_identical(unname(as.array(counts(res))), exp)

    ## half="lower"
    res <- pullHicMatrices(gi, hicFiles[1], 2.5e06, half="lower")
    expect_identical(interactions(res), gi)
    expect_identical(unname(as.array(counts(res))), exp)


    ## Interchromosomal
    ## Pull rectangle where seqnames1 > seqnames2
    gi <- read.table(text="
            2 50000000 57500000 1 50000000 62500000
            2 150000000 157500000 1 150000000 162500000") |>
        as_ginteractions() |>
        suppressWarnings()
    expect_true(all(seqnames1(gi) > seqnames2(gi)))

    ## half="both"
    res <- pullHicMatrices(gi, hicFiles[1], 2.5e06, half="both")
    exp <- array(data=c(1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0,
                        2, 0, 2, 0, 0, 0, 2, 1, 2, 1, 0, 2, 0,
                        1, 1, 1, 1),
                 dim=c(3,5,2,1))
    expect_identical(interactions(res), gi)
    expect_identical(unname(as.array(counts(res))), exp)

    ## half="upper"
    res <- pullHicMatrices(gi, hicFiles[1], 2.5e06, half="upper")
    expect_identical(interactions(res), gi)
    expect_identical(unname(as.array(counts(res))), exp)

    ## half="lower"
    res <- pullHicMatrices(gi, hicFiles[1], 2.5e06, half="lower")
    expect_identical(interactions(res), gi)
    expect_identical(unname(as.array(counts(res))), exp)

})


test_that("Pull irregular arrays", {

    gi <- read.table(text="
            1 51000000 51300000 1 51000000 51500000
            2 52000000 52300000 3 52000000 52500000
            1 150000000 150500000 1 150000000 150300000
            2 52000000 52300000 2 52000000 52800000") |>
        as_ginteractions()

    iset <- pullHicMatrices(gi, hicFiles, 100e03, half="both")

    ## Are counts pulled & accessed correctly?
    iset1 <- pullHicMatrices(gi[1], hicFiles[1], 50e03, half="both")
    iset2 <- pullHicMatrices(gi, hicFiles[1], 50e03, half="both")
    exp <- counts(iset1) |> as.matrix()
    expect_identical(counts(iset2)[1,1] |> as.matrix(), exp)
    expect_identical(as.list(counts(iset2))[[1]][[1]], exp)

})
