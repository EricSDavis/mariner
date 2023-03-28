library(mariner)

## Shared objects --------------------------------------------------------------

## Read .hic file paths
hicFiles <-
    system.file("extdata/test_hic", package="mariner") |>
    list.files(pattern=".hic", full.names=TRUE)

## Read in loops as GInteractions object
loops <-
    system.file("extdata", package="mariner") |>
    list.files(pattern="WT.*Loops.txt", full.names=TRUE) |>
    read.table(header=TRUE) |>
    as_ginteractions(keep.extra.columns=FALSE)

## Removes the "chr" prefix for compatibility
## with the preprocessed hic files
GenomeInfoDb::seqlevelsStyle(loops) <- 'ENSEMBL'

## Expand pixel ranges with a 5 pixel buffer on either side
loops <-
    binPairs(loops, binSize=100e3) |>
    pixelsToMatrices(buffer=5)

## Extract 10, 11x11 count matrices from 2 hic files
iarr <-
    loops[1:10] |>
    pullHicMatrices(binSize=100e3,
                    files=hicFiles)

## Tests -----------------------------------------------------------------------

test_that("aggHicMatrices works as expected", {

    ## Aggregate everything (by=NULL)
    ans <- aggHicMatrices(x=iarr, verbose=FALSE)
    expect_s4_class(ans, "DelayedMatrix")
    expect_identical(dim(ans), c(11L, 11L))

    ## Aggregate everything (by="files")
    ans <- aggHicMatrices(x=iarr, by="files", verbose=FALSE)
    expect_s4_class(ans, "DelayedArray")
    expect_identical(dim(ans), c(11L, 11L, 2L))

    ## Aggregate everything (by="interactions")
    ans <- aggHicMatrices(x=iarr, by="interactions", verbose=FALSE)
    expect_s4_class(ans, "DelayedArray")
    expect_identical(dim(ans), c(11L, 11L, 10L))

})
