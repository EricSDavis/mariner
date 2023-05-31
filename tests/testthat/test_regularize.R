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
gi <- c(gi,gi) # make more interactions

## InteractionJaggedArray object
ija <- pullHicMatrices(gi, hicFiles, 100e03, half="both")

test_that("regularize JaggedArray and InteractionJaggedArray", {

    ## JaggedArray
    ja <- counts(ija)
    expect_identical(dim(regularize(ja)), c(10L, 10L, 8L, 2L))
    expect_no_error(regularize(ja, nBlocks=1))
    expect_error(regularize(ja, nBlocks=0), "must be > 0")
    expect_identical(regularize(ja, scale=FALSE)[1,1,1,1], 53)

    ## Works with subset JaggedArray
    expect_equal(regularize(ja[,,3:1,], scale=FALSE)[1,1,1,1], 63)
    expect_equal(regularize(ja[,,,1], scale=FALSE)[1,1,1,1], 53)

    ## Blocks aggregate correctly
    expect_identical(
        apply(regularize(ja, c(3,3)), c(1,2), sum, na.rm=TRUE),
        apply(regularize(ja, c(3,3), nBlocks=1), c(1,2), sum, na.rm=TRUE)
    )

    ## InteractionJaggedArray can be regularized
    ## than aggregated with aggHicMatrices
    ia <- regularize(ija)
    expect_error(print(counts(ia, showDimnames=TRUE)),
                 "Dimnames not available")
    agg <-
        regularize(ija, nBlocks=1) |>
        aggHicMatrices(nBlocks=1, verbose=FALSE)
    expect_s4_class(agg, "DelayedMatrix")
    expect_equal(dim(agg), c(10L, 10L))

    agg <-
        regularize(ija, nBlocks=1) |>
        aggHicMatrices(nBlocks=1, by='files', verbose=FALSE)
    expect_s4_class(agg, "DelayedArray")
    expect_equal(dim(agg), c(10L, 10L, 2L))

    agg <-
        regularize(ija, nBlocks=1) |>
        aggHicMatrices(nBlocks=1, by='interactions', verbose=FALSE)
    expect_s4_class(agg, "DelayedArray")
    expect_equal(dim(agg), c(10L, 10L, 8L))

    ## Works with subset InteractionJaggedArray
    sub <- ija[3:1, 1] |> regularize(scale=FALSE) |> counts()
    expect_s4_class(sub, "DelayedArray")
    expect_equal(dim(sub), c(10L, 10L, 3L, 1L))
    expect_equal(sub[1,1,3,1], c("FS"=53))

})

test_that("Interpolation works with single row/column", {
    
    ## Create regions of dimensions 1x3, 3x1, 1x1
    gi <- read.table(
        text="
    1 51000000 51100000 1 51000000 51300000
    1 51000000 51300000 1 51000000 51100000 
    1 51000000 51100000 1 51000000 51100000
    "
    )
    gi <- as_ginteractions(gi)
    
    ## InteractionJaggedArray
    ija <- pullHicMatrices(x=gi, files=hicFiles[1], binSize=100e3)
    
    ## Internally this should create 2x3, 3x2, 2x2
    ## to avoid interpolation errors.
    expect_s4_class(regularize(ija, ndim=c(3,3)), "InteractionArray")
    
    regularize(ija, ndim=c(3,3), scale=FALSE) |> counts() |>
        expect_snapshot()
    
    ## Not sure what behavior we want, but this is
    ## something...
    regularize(ija, ndim=c(1,3), scale=FALSE) |>
        counts() |>
        expect_snapshot()
    regularize(ija, ndim=c(3,1), scale=FALSE) |>
        counts() |>
        expect_snapshot()
    
})