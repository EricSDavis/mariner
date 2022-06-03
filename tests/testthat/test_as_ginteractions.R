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

test_that("as_ginteractions works with data.tables", {

    ## Construct GInteractions from data.table
    library(data.table)
    gi2 <-
        data.table(chr1 = "chr1", x1 = 10000, x2 = 20000,
                   chr2 = "chr1", y1 = 30000, y2 = 40000) |>
        as_ginteractions()

    expect_identical(gi2, gi1)

})

test_that("as_ginteractions works with DataFrames", {

    ## Construct GInteractions from DataFrame
    library(S4Vectors)
    gi2 <-
        DataFrame(chr1 = "chr1", x1 = 10000, x2 = 20000,
                  chr2 = "chr1", y1 = 30000, y2 = 40000) |>
        as_ginteractions()

    expect_identical(gi2, gi1)

})

test_that("alias (makeGInteractionsFromDataFrame) works", {

    ## Construct GInteractions from DataFrame
    library(S4Vectors)
    gi2 <-
        DataFrame(chr1 = "chr1", x1 = 10000, x2 = 20000,
                  chr2 = "chr1", y1 = 30000, y2 = 40000) |>
        makeGInteractionsFromDataFrame()

    expect_identical(gi2, gi1)
})

test_that("Additional metadata can be added in as_ginteractions", {

    ## Construct GInteractions from DataFrame
    gi2 <-
        data.frame(chr1 = "chr1", x1 = 10000, x2 = 20000,
                   chr2 = "chr1", y1 = 30000, y2 = 40000,
                   pval = 0.05, dist = 10000) |>
        as_ginteractions()

    ## Add metadata to gi1 object
    gi3 <- gi1 |> `mcols<-`(value = DataFrame(pval = 0.05, dist = 10000))

    expect_identical(gi2, gi3)

})

test_that("Additional metadata can be ignored in as_ginteractions", {

    ## Construct GInteractions from DataFrame
    gi2 <-
        data.frame(chr1 = "chr1", x1 = 10000, x2 = 20000,
                   chr2 = "chr1", y1 = 30000, y2 = 40000,
                   pval = 0.05, dist = 10000) |>
        as_ginteractions(keep.extra.columns = FALSE)


    expect_identical(gi2, gi1)

})

test_that("Require at least 6 columns for as_ginteractions", {

    data.frame(chr1 = "chr1", x1 = 10000, x2 = 20000,
               chr2 = "chr1", y1 = 30000) |>
        as_ginteractions() |>
        expect_error()

})

test_that("Exported tables can be read with as_ginteractions", {

    ## Make tempdir for testing
    tf <- tempfile(fileext = ".txt")

    ## Write out gi1 object
    write.table(x = gi1, file = tf)

    ## Read back in as GInteractions
    gi2 <-
        read.table(file = tf) |>
        as_ginteractions()

    expect_identical(gi2, gi1)

})

test_that("as_ginteractions supports 10-column format", {

    ## Construct GInteractions from data.frame
    gi2 <-
        data.frame(seqnames1 = "chr1", start1 = 10000, end1 = 20000,
                   width1 = 10001, strand1 = '*',
                   seqnames2 = "chr1", start2 = 30000, end2 = 40000,
                   width2 = 10001, strand2 = '*') |>
        as_ginteractions()

    expect_identical(gi2, gi1)

})

test_that("as_ginteractions supports 10-column format (chr)", {

    ## Construct GInteractions from data.frame
    gi2 <-
        data.frame(chr1 = "chr1", start1 = 10000, end1 = 20000,
                   width1 = 10001, strand1 = '*',
                   chr2 = "chr1", start2 = 30000, end2 = 40000,
                   width2 = 10001, strand2 = '*') |>
        as_ginteractions()

    expect_identical(gi2, gi1)

})

test_that("as_ginteractions supports 10-column format (chrom)", {

    ## Construct GInteractions from data.frame
    gi2 <-
        data.frame(chrom1 = "chr1", start1 = 10000, end1 = 20000,
                   width1 = 10001, strand1 = '*',
                   chrom2 = "chr1", start2 = 30000, end2 = 40000,
                   width2 = 10001, strand2 = '*') |>
        as_ginteractions()

    expect_identical(gi2, gi1)

})

test_that("as_ginteractions supports 10-column format (chrom)", {

    ## Construct GInteractions from data.frame
    gi2 <-
        data.frame(chrom1 = "chr1", start1 = 10000, end1 = 20000,
                   width1 = 10001, strand1 = '*',
                   chrom2 = "chr1", start2 = 30000, end2 = 40000,
                   width2 = 10001, strand2 = '*') |>
        as_ginteractions()

    expect_identical(gi2, gi1)

})
