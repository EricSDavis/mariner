library(mariner)

## Shared objects --------------------------------------------------------------

library(GenomicRanges)

## Create example GRanges
gr1 <- GRanges(seqnames = "chr1",
               ranges = IRanges(start = rep(5000,3),
                                end = rep(6000,3)),
               strand = c('+', '-', '*'))

gr2 <- gr1 |> promoters(upstream = 2000, downstream = 200)


## Tests -----------------------------------------------------------------------

test_that("Binning works with keyword start", {
    exp <- GRanges(seqnames = "chr1",
                   ranges = IRanges(start = c(5000, 6000, 5000),
                                    end = c(6000, 7000, 6000)),
                   strand = c('+', '-', '*'))
    obj <- binRanges(gr1, 1000, 'start')
    expect_identical(obj, exp)

})

test_that("Binning works with keyword end", {
    exp <- GRanges(seqnames = "chr1",
                   ranges = IRanges(start = c(6000, 5000, 6000),
                                    end = c(7000, 6000, 7000)),
                   strand = c('+', '-', '*'))
    obj <- binRanges(gr1, 1000, 'end')
    expect_identical(obj, exp)

})

test_that("Binning works with keyword center", {
    exp <- GRanges(seqnames = "chr1",
                   ranges = IRanges(start = c(5000, 5000, 5000),
                                    end = c(6000, 6000, 6000)),
                   strand = c('+', '-', '*'))
    obj <- binRanges(gr1, 1000, 'center')
    expect_identical(obj, exp)

})

test_that("Binning to transcription start site", {
    exp <- GRanges(seqnames = "chr1",
                   ranges = IRanges(start = c(4000, 6000, 4000),
                                    end = c(6000, 8000, 6000)),
                   strand = c('+', '-', '*'))
    obj <- binRanges(gr2, 2000, 1000)
    expect_identical(obj, exp)
})
