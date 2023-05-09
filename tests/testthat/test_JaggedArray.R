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

test_that("JaggedArray accessors", {

    ## path accessor
    cnts <- counts(iarr)
    path(cnts)

    ## dim accessor
    cnts <- counts(iarr)
    dim(cnts)

    ## Only show break for more than 2 matrices
    ja <- counts(iarr)
    expect_snapshot(ja[,,c(1,3),1])

    ## Extracting & subsetting counts for JaggedArray ####
    ## Possible to extract counts
    cnts <- counts(iarr)
    expect_snapshot(cnts)
    expect_identical(as.list(cnts)[[1]][[1]][1,1], 53)
    expect_identical(as.list(cnts)[[1]][[4]][1,1], 29)

    ## Reversing by subsetting works
    cnts <- counts(iarr)[,,4:1, 1]
    expect_snapshot(cnts)
    expect_identical(as.list(cnts)[[1]][[1]][1,1], 29)
    expect_identical(as.list(cnts)[[1]][[4]][1,1], 53)

    ## Successive subsetting works
    cnts2 <- cnts[,,4:1, 1]
    expect_snapshot(cnts2)
    expect_identical(as.list(cnts2)[[1]][[1]][1,1], 53)
    expect_identical(as.list(cnts2)[[1]][[4]][1,1], 29)

    ## These should auto convert to DelayedArray
    expect_s4_class(counts(iarr)[,,1:2,], "DelayedArray")
    expect_s4_class(counts(iarr)[,,1,], "DelayedArray")

    ## None of these should error
    expect_s4_class(counts(iarr)[,,1,], "DelayedArray")
    expect_s4_class(counts(iarr)[,,,1], "JaggedArray")
    expect_s4_class(counts(iarr)[], "JaggedArray")
    expect_s4_class(counts(iarr)[,,1:2, 1], "DelayedArray")
    expect_s4_class(counts(iarr)[,,2:1, 1], "DelayedArray")

    ## Test for correct subscripting
    ja <- counts(iarr)
    expect_error(ja[,], "incorrect.*subscripts")
    expect_error(ja[,,], "incorrect.*subscripts")
    expect_error(ja[,,], "incorrect.*subscripts")
    expect_error(ja[,,,], NA)

    ## Error for trying to subset the first two dimensions
    ja <- counts(iarr)
    expect_error(ja[1,1,1,1], "Indexing.*not supported")
    expect_error(ja[1,,1,1], "Indexing.*not supported")
    expect_error(ja[,1,1,1], "Indexing.*not supported")

})
