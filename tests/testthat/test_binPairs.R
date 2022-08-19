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
