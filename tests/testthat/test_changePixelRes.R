library(mariner)
library(marinerData)

## Read .hic file paths
hicFiles <- c(
    LEUK_HEK_PJA27_inter_30.hic(),
    LEUK_HEK_PJA30_inter_30.hic()
)
names(hicFiles) <- c("FS", "WT")

## Read in loops as GInteractions object
loops <-
    WT_5kbLoops.txt() |>
    setNames("WT") |>
    read.table(header=TRUE) |>
    as_ginteractions(keep.extra.columns=FALSE)

## Removes the "chr" prefix for compatibility
## with the preprocessed hic files
GenomeInfoDb::seqlevelsStyle(loops) <- 'ENSEMBL'

## Rebin loops to 2.5e6 resolution
loops <- binPairs(x=loops, binSize=2.5e06)

test_that("selectPixel chooses correctly", {

    ## Pull 5x5 matrices
    iarr <- pullHicMatrices(x=loops[1:5],
                            binSize=500e3,
                            files=hicFiles,
                            norm="KR",
                            half='upper')

    ## Select pixel
    expect_equal(round(selectPixel(iarr)$value),
                 c(543, 504, 465, 720, 541))

})

test_that("changePixelRes works as expected", {

    newPixels <-
        changePixelRes(x=loops[1:5],
                       files=hicFiles,
                       from=2.5e6,
                       to=500e3,
                       half="upper")

    expect_equal(round(newPixels$value),
                 c(543, 504, 465, 720, 541))

    expect_identical(unique(unlist(width(newPixels))), 500001L)

})
