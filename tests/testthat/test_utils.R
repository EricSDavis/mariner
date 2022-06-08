library(mariner)
## Shared objects --------------------------------------------------------------
library(GenomicRanges)

## Define example GRanges and binSize
gr <- GRanges("chr1:1-10000")
bs <- 10e03

## Perform binning
br <- GRanges("chr1:0-10000")


## Test .checkBinnedRanges() ---------------------------------------------------

test_that("Checking for binned ranges works", {

    .checkBinnedRanges(x = gr, binSize = bs) |>
        expect_false()

    .checkBinnedRanges(x = br, binSize = bs) |>
        expect_true()

    .checkBinnedRanges(x = c(gr, br), binSize = bs) |>
        expect_false()

    .checkBinnedRanges(x = c(br, br), binSize = bs) |>
        expect_true()

})

test_that("Checking for binned pairs works", {

    .checkBinnedPairs(x = GInteractions(gr, gr),
                      binSize = bs) |>
        expect_false()

    .checkBinnedPairs(x = GInteractions(br, br),
                      binSize = bs) |>
        expect_true()

    .checkBinnedPairs(x = GInteractions(br, gr),
                      binSize = bs) |>
        expect_false()

    .checkBinnedPairs(x = GInteractions(c(br, br), c(br, br)),
                      binSize = bs) |>
        expect_true()

    .checkBinnedPairs(x = GInteractions(c(gr, gr), c(gr, gr)),
                      binSize = bs) |>
        expect_false()

    .checkBinnedPairs(x = GInteractions(c(br, gr), c(gr, gr)),
                      binSize = bs) |>
        expect_false()

})
