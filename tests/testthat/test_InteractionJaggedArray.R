library(mariner)
library(marinerData)

## Read .hic file paths
hicFiles <- c(
    LEUK_HEK_PJA27_inter_30.hic(),
    LEUK_HEK_PJA30_inter_30.hic()
)
names(hicFiles) <- c("FS", "WT")

## Create test interactions
gi <- read.table(text="
            1 51000000 51300000 1 51000000 51500000
            2 52000000 52300000 3 52000000 52500000
            1 150000000 150500000 1 150000000 150300000
            2 52000000 52300000 2 52000000 52800000") |>
    as_ginteractions()

## InteractionJaggedArray object
iarr <- pullHicMatrices(gi, hicFiles, 100e03, half="both")

test_that("InteractionJaggedArray accessors", {

    ## InteractionJaggedArray Accessors
    expect_s4_class(interactions(iarr), "GInteractions")
    expect_identical(metadata(iarr),
                     list(binSize=1e5, norm="NONE", matrix="observed"))
    expect_identical(nrow(colData(iarr)), 2L)
    expect_identical(length(dim(iarr)), 4L)
    expect_snapshot(counts(iarr))
    expect_error(path(iarr), NA)

})
