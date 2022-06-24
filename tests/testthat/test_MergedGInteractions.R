library(mariner)
library(GenomicRanges)
library(InteractionSet)
library(data.table)
library(S4Vectors)

## Shared objects --------------------------------------------------------------

## Define example anchor regions
gr1 <-
    GRanges(seqnames = "chr1",
            ranges = IRanges(start = c(30,40,40,70,80),
                             end = c(40,50,50,80,90)))
gr2 <-
    GRanges(seqnames = "chr1",
            ranges = IRanges(start = c(30,30,50,10,30),
                             end = c(40,40,60,20,40)))

## Form GInteractions and convert to data.table
dt <- GInteractions(gr1, gr2) |> as.data.table()

## Split into two files
dts <- split(dt, c(rep(1,3), rep(2, 2)))

## Reference BEDPE files (loops called with SIP)
bedpeFiles <-
    system.file("extdata", package = "mariner") |>
    list.files(pattern = "Loops.txt", full.names = TRUE)


## Test getters and setters ----------------------------------------------------

test_that("selectionMethod accessor works", {
    x <- mergePairs(bedpeFiles)
    selectionMethod(x) |>
        expect_identical("Mean of modes")

    x <- mergePairs(bedpeFiles, column = "APScoreAvg")
    selectionMethod(x) |>
        expect_identical("Select by column 'APScoreAvg'")
})

test_that("allPairs accessor works", {

    ## Merge pairs and add names
    x <- mergePairs(bedpeFiles)
    names(x) <- paste0("loop", 1:length(x))

    x[3:1] |>
        allPairs() |>
        names() |>
        expect_identical(paste0("loop", 3:1))

    x[1:3] |>
        allPairs() |>
        names() |>
        expect_identical(paste0("loop", 1:3))

    expect_equal(length(allPairs(x)), length(x))


    ## Use data.table input
    x <- mergePairs(dts, binSize = 5)
    expect_equal(length(allPairs(x)), length(x))
})


test_that("deNovo method works", {

    ## Merge pairs
    x <- mergePairs(dts, binSize = 10, radius = 2)

    ## deNovo pairs
    lapply(deNovo(x), length) |>
        unname() |>
        unlist() |>
        expect_equal(c(2,1))

})


test_that("Aggregating metadata columns works", {

    ## Merge pairs
    x <- mergePairs(x = bedpeFiles,
                    binSize = 5e03,
                    radius = 0)

    aggPairMcols(x, columns = c("APScoreAvg", "avg"), fun = "mean")

    ## Testing character vs function input
    aggPairMcols(x, columns = "APScoreAvg", funs = "mean") |>
        mcols() |>
        colnames() |>
        expect_identical("mean.APScoreAvg")

    aggPairMcols(x, columns = "APScoreAvg", funs = mean) |>
        mcols() |>
        colnames() |>
        expect_identical("fun1.APScoreAvg")

    ## Testing values
    expect_identical(object = aggPairMcols(x[1:10],
                                           columns = "APScoreAvg",
                                           funs = "mean")$mean.APScoreAvg,
                     expected = allPairs(x[1:10]) |>
                         lapply(`[[`, "APScoreAvg") |>
                         lapply(mean) |>
                         unlist())

    ## Multiple input types
    aggPairMcols(x,
                 columns = "APScoreAvg",
                 fun = c("mean", min, \(x) mean(x))) |>
        mcols() |>
        colnames() |>
        expect_identical(c("mean.APScoreAvg",
                           "fun2.APScoreAvg",
                           "fun3.APScoreAvg"))

    aggPairMcols(x,
                 columns = c("APScoreAvg", "avg"),
                 fun = "mean") |>
        mcols() |>
        colnames() |>
        expect_identical(c("mean.APScoreAvg",
                           "mean.avg"))

    ## Throw errors
    aggPairMcols(x, columns = "APScoreAvg", funs = \(x) x) |>
        expect_error("Improper aggregation function.")

    aggPairMcols(x, columns = c("blah", "foo"), funs = "mean") |>
        expect_error("^Column.*not exist.$")

})
