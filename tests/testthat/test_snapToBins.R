library(mariner)
library(InteractionSet)
library(GenomicRanges)

test_that("GRanges snap correctly to bins", {
    x <- GRanges(seqnames = c("chr1"),
                 ranges = IRanges(start = c(1, 1, 25, 19, 21),
                                  end = c(15, 11, 31, 31, 39)))

    snapToBins(x, binSize = 5) |>
        expect_identical(
            GRanges(seqnames = c("chr1"),
                    ranges = IRanges(start = c(0, 0, 25, 20, 20),
                                     end = c(15, 10, 30, 30, 40)))
        )
    snapToBins(x, binSize = 10) |>
        expect_identical(
            GRanges(seqnames = c("chr1"),
                    ranges = IRanges(start = c(0, 0, 20, 20, 20),
                                     end = c(20, 10, 30, 30, 40)))
        )
    snapToBins(x, binSize = 20) |>
        expect_identical(
            GRanges(seqnames = c("chr1"),
                    ranges = IRanges(start = c(0, 0, 20, 20, 20),
                                     end = c(20, 20, 40, 40, 40)))
        )
})

test_that("GInteractions correctly snap to bins", {

    ## Sample GInteractions object
    x <- GInteractions(anchor1 = c(GRanges("chr1:1-15"),
                                   GRanges("chr1:1-11")),
                       anchor2 = c(GRanges("chr1:25-31"),
                                   GRanges("chr1:19-31")))

    snapToBins(x, binSize = 5) |>
        expect_identical(
            GInteractions(c(GRanges("chr1:0-15"),
                            GRanges("chr1:0-10")),
                          c(GRanges("chr1:25-30"),
                            GRanges("chr1:20-30")))
        )
    snapToBins(x, binSize = 10) |>
        expect_identical(
            GInteractions(c(GRanges("chr1:0-20"),
                            GRanges("chr1:0-10")),
                          c(GRanges("chr1:20-30"),
                            GRanges("chr1:20-30")))
        )
    snapToBins(x, binSize = 20) |>
        expect_identical(
            GInteractions(c(GRanges("chr1:0-20"),
                            GRanges("chr1:0-20")),
                          c(GRanges("chr1:20-40"),
                            GRanges("chr1:20-40")))
        )

    ## binSize must be > 0
    snapToBins(x, binSize = 0) |>
        expect_error("`binSize` must be > 0")

    ## Extended GInteractions object
    mgi <- mergePairs(x, radius = 0)

    snapToBins(mgi, binSize = 5) |>
        expect_s4_class("GInteractions")

    snapToBins(mgi, binSize = 5) |>
        expect_s4_class("MergedGInteractions")

})
