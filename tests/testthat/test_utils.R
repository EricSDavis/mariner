library(mariner)
## Shared objects --------------------------------------------------------------
library(GenomicRanges)
library(InteractionSet)
library(rlang)
library(glue)

## Define example GRanges and binSize
gr <- GRanges("chr1:1-10000")
bs <- 10e03

## Perform binning
br <- GRanges("chr1:0-10000")

## Test .checkSnapped* ---------------------------------------------------------

test_that("Checking for snapped GRanges works", {

    x <- GRanges(seqnames = "chr1",
                 ranges = IRanges(start = c(0, 2),
                                  end = c(50, 48)))

    .checkSnappedRanges(x, 50) |>
        expect_false()

    .checkSnappedRanges(x, 2) |>
        expect_true()

})


test_that("Checking for snapped GInteractions works", {

    x <- GInteractions(anchor1 = c(GRanges("chr1:0-50"),
                                   GRanges("chr1:2-48")),
                       anchor2 = c(GRanges("chr1:0-50"),
                                   GRanges("chr1:2-48")))

    .checkSnappedPairs(x, 50) |>
        expect_false()

    .checkSnappedPairs(x, 2) |>
        expect_true()

})


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

## Test .modes() ---------------------------------------------------------------

test_that(".modes returns correct results", {
    .modes(c(1,1,2,1)) |>
        expect_equal(1)

    .modes(c(1,2,2,1,3)) |>
        expect_equal(c(1,2))

    .modes(c(1)) |>
        expect_equal(1)
})

## Test .checkVectorLengths ----------------------------------------------------

test_that("Function to check vector lengths works.", {

    .checkVectorLengths(list(arg1=1)) |>
        expect_null()

    .checkVectorLengths(list(arg1=c(1,2))) |>
        expect_error(".*arg1.*not.*supported.*")

    .checkVectorLengths(list(arg1=1, arg2="a")) |>
        expect_null()

    .checkVectorLengths(list(arg1=1, arg2="a", arg3=TRUE)) |>
        expect_null()

    .checkVectorLengths(list(arg1=1, arg2="a",
                             arg3=c(NA, TRUE), arg4=c(NA, TRUE))) |>
        expect_error(".*arg3.*not.*supported.*")

})
