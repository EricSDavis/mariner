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

test_that("shiftRanges works with keywords (start)", {
    exp <-
        GRanges(seqnames = "chr1",
                ranges = IRanges(start = c(5000, 6000, 5000)),
                strand = c('+', '-', '*'))

    obj <- shiftRanges(gr1, 'start')

    expect_identical(obj, exp)

})

test_that("shiftRanges works with keywords (end)", {
    exp <-
        GRanges(seqnames = "chr1",
                ranges = IRanges(start = c(6000, 5000, 6000)),
                strand = c('+', '-', '*'))

    obj <- shiftRanges(gr1, 'end')

    expect_identical(obj, exp)

})

test_that("shiftRanges works with keywords (center)", {
    exp <-
        GRanges(seqnames = "chr1",
                ranges = IRanges(start = c(5500, 5500, 5500)),
                strand = c('+', '-', '*'))

    obj <- shiftRanges(gr1, 'center')

    expect_identical(obj, exp)

})

test_that("shiftRanges only accepts keywords 'start', 'center', 'end'", {
    shiftRanges(GRanges(), pos = 'blah') |>
        expect_error()
    shiftRanges(GRanges(), pos = 's') |>
        expect_error()
})

test_that("shiftRanges works with an integer", {
    exp <-
        GRanges(seqnames = "chr1",
                ranges = IRanges(start = c(5100, 5900, 5100)),
                strand = c('+', '-', '*'))

    obj <- shiftRanges(gr1, 100)

    expect_identical(obj, exp)
})

test_that("shiftRanges works with an integer (promoter example)", {
    exp <-
        GRanges(seqnames = "chr1",
                ranges = IRanges(start = c(5000, 6000, 5000)),
                strand = c('+', '-', '*'))
    obj <- shiftRanges(gr2, c(2000))

    expect_identical(obj, exp)
})

test_that("shiftRanges works with integers", {
    exp <-
        GRanges(seqnames = "chr1",
                ranges = IRanges(start = c(5100, 5800, 5300)),
                strand = c('+', '-', '*'))

    obj <- shiftRanges(gr1, c(100, 200, 300))

    expect_identical(obj, exp)
})

test_that("Incorrect pos argument triggers error", {
    shiftRanges(gr1, c(100, 200)) |>
        expect_error()
    shiftRanges(gr1, c(100, 200, 300, 400)) |>
        expect_error()
})
