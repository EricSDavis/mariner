library(mariner)

test_that("Constructors, accessors, and transformations work", {

    ## Generate example
    library(InteractionSet)
    library(GenomicRanges)
    example(GInteractions, echo=FALSE)

    ## Create GInteractions of different widths
    gi2 <- GInteractions(resize(first(gi), width = 10),
                         second(gi))
    gi3 <- GInteractions(resize(first(gi), width = 10),
                         resize(second(gi), width = 10))
    gi4 <- GInteractions(resize(first(gi), width = 10),
                         resize(second(gi), width = 20))

    ## Test constructor
    BinnedGInteractions(gi) |>
        expect_error("First pair.*equal widths.")
    BinnedGInteractions(gi2) |>
        expect_error("Second pair.*equal widths.")
    validObject(BinnedGInteractions(gi3)) |>
        expect_true()
    validObject(BinnedGInteractions(gi4)) |>
        expect_true()

    ## Test accessors
    bgi <- BinnedGInteractions(gi3)
    expect_equal(firstBinSize(bgi), 10L)
    expect_equal(secondBinSize(bgi), 10L)
    expect_true(pairBinsEqual(bgi))

    ## Test trim
    bgi <- BinnedGInteractions(gi3)
    suppressWarnings({seqlengths(bgi) <- c("chrA"=670, "chrB" = 700)})
    trim(bgi) |>
        expect_error("First pair.*equal widths.")

    ## Test resize
    bgi <- BinnedGInteractions(gi3)
    bgi <- resize(bgi, width = 5, fix = 'start', use.names = TRUE)
    expect_equal(firstBinSize(bgi), 5L)
    expect_equal(secondBinSize(bgi), 5L)
    expect_true(pairBinsEqual(bgi))

    ## Test narrow
    bgi <- BinnedGInteractions(gi3)
    bgi <- narrow(bgi, start = 5, end = 10)
    expect_equal(firstBinSize(bgi), 6L)
    expect_equal(secondBinSize(bgi), 6L)
    expect_true(pairBinsEqual(bgi))

    ## Test flank
    bgi <- BinnedGInteractions(gi3)
    bgi <- flank(bgi, width = 200L)
    expect_equal(firstBinSize(bgi), 200L)
    expect_equal(secondBinSize(bgi), 200L)
    expect_true(pairBinsEqual(bgi))

    ## Test regions<-
    bgi <- BinnedGInteractions(gi3)
    regions(bgi) <- resize(regions(bgi), 100L)
    expect_equal(firstBinSize(bgi), 100L)
    expect_equal(secondBinSize(bgi), 100L)
    expect_true(pairBinsEqual(bgi))

    ## Test swapAnchors error
    bgi <- BinnedGInteractions(gi4)
    expect_equal(firstBinSize(bgi), 10L)
    expect_equal(secondBinSize(bgi), 20L)
    expect_false(pairBinsEqual(bgi))
    swapAnchors(bgi) |>
        expect_error("First pair.*equal widths.")

})
