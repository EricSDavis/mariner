library(mariner)

## Shared objects --------------------------------------------------------------

## Construct GInteractions
library(InteractionSet)
gi1 <-
    data.frame(chr1 = "chr1", x1 = 10000, x2 = 20000,
               chr2 = "chr1", y1 = 30000, y2 = 40000) |>
    as_ginteractions()


## Tests -----------------------------------------------------------------------

test_that("binPairs can expand the resolution (from start)", {

    ## Construct GInteractions from data.frame
    exp <-
        data.frame(chr1 = "chr1", x1 = 0, x2 = 20000,
                   chr2 = "chr1", y1 = 20000, y2 = 40000) |>
        as_ginteractions()

    obj <-
        binPairs(x = gi1,
                 binSize = 20000,
                 pos1 = 'start',
                 pos2 = 'start')

    expect_identical(obj, exp)

})

test_that("Bin to the TSS of one anchor", {

    gi2 <-
        GInteractions(anchor1 = anchors(gi1, 'first'),
                      anchor2 = gi1 |>
                          anchors('second') |>
                          promoters(upstream = 2e03, downstream = 200))

    ## Construct GInteractions from data.frame
    exp <-
        data.frame(chr1 = "chr1", x1 = 10000, x2 = 11000,
                   chr2 = "chr1", y1 = 30000, y2 = 31000) |>
        as_ginteractions()

    obj <-
        binPairs(x = gi2,
                 binSize = 1000,
                 pos1 = 'start',
                 pos2 = 2000)

    expect_identical(obj, exp)

})

test_that("binPairs retains metadata", {

    ## Construct GInteractions from data.frame
    exp <-
        data.frame(chr1 = "chr1", x1 = 0, x2 = 20000,
                   chr2 = "chr1", y1 = 20000, y2 = 40000,
                   var1 = "var1", var2 = 10) |>
        as_ginteractions()

    gi1$var1 <- "var1"
    gi1$var2 <- 10

    obj <-
        binPairs(x = gi1,
                 binSize = 20000,
                 pos1 = 'start',
                 pos2 = 'start')

    expect_identical(obj, exp)

})

test_that("binPairs doesn't modify class", {

    ## Load required packages
    library(data.table)

    ## Reference BEDPE files (loops called with SIP)
    bedpeFiles <-
        system.file("extdata", package = "mariner") |>
        list.files(pattern = "Loops.txt", full.names = TRUE)

    ## Assemble list of GInteractions
    giList <-
        lapply(bedpeFiles, fread) |>
        lapply(as_ginteractions)

    ## Merge pairs
    mgi <- mergePairs(x = giList, radius = 10e03)

    ## Add example metadata
    mgi$metadata1 <- "foo"

    ## Binning
    bmgi <-
        binPairs(x = mgi,
                 binSize = 10e03,
                 pos1 = 'start',
                 pos2 = 'end')

    expect_true(is(bmgi, "MergedGInteractions"))
    expect_equal(ncol(mcols(bmgi)), 1)


    ## For non merged data
    bp <-
        binPairs(x = giList[[1]],
                 binSize = 10e03,
                 pos1 = 'start',
                 pos2 = 'end')

    expect_true(is(bp, "GInteractions"))
    expect_equal(ncol(mcols(bp)), 9)

})

test_that("Binning checks for appropriate binSize", {
    binPairs(gi1, binSize = 0) |>
        expect_error("`binSize` must be > 0")
})

test_that("Anchors can be binned independently", {
    bgi <- binPairs(gi1, binSize=c(5e03, 1e03))
    expect_equal(unique(width(regions(bgi)))-1,
                 c(5e03, 1e03))
})
