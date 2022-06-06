library(mariner)

## Shared objects --------------------------------------------------------------

## Construct GRanges
library(GenomicRanges)
gr1 <-
    GRanges(seqnames = 'chr1',
            ranges = IRanges(start = 10000,
                             end = 20000))
gr2 <-
    GRanges(seqnames = 'chr1',
            ranges = IRanges(start = 30000,
                             end = 40000))

## Construct GInteractions
library(InteractionSet)
gi1 <-
    GInteractions(anchor1 = gr1,
                  anchor2 = gr2)

## Tests -----------------------------------------------------------------------

test_that("as_ginteractions works with data.frames", {

    ## Construct GInteractions from data.frame
    gi2 <-
        data.frame(chr1 = "chr1", x1 = 10000, x2 = 20000,
                   chr2 = "chr1", y1 = 30000, y2 = 40000) |>
        as_ginteractions()

    expect_identical(gi2, gi1)

})
