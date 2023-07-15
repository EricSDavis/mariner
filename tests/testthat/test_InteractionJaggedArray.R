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
    expect_identical(length(iarr), 4L)

})

test_that("InteractionJaggedArray subsetting", {

    ## i
    sub <- iarr[1:3,]
    expect_identical(counts(sub), counts(iarr)[,,1:3,])
    expect_identical(dim(sub), list(
        interactions=3L,
        files=2L,
        rows=c(3, 3, 5),
        cols=c(5, 5, 3)
    ))
    expect_identical(
        interactions(sub),
        interactions(iarr)[1:3,]
    )

    ## i reversed
    sub2 <- iarr[3:1,]
    expect_identical(counts(sub2), counts(iarr)[,,3:1,])
    expect_identical(dim(sub2), list(
        interactions=3L,
        files=2L,
        rows=c(5,3,3),
        cols=c(3,5,5)
    ))
    expect_identical(
        interactions(sub2),
        interactions(iarr)[3:1,]
    )

    ## missing i and j
    expect_identical(
        as.list(counts(iarr)),
        as.list(counts(iarr[]))
    )

    ## both i and j
    sub2 <- iarr[3:1,1]
    expect_identical(counts(sub2), counts(iarr)[,,3:1,1])
    expect_identical(dim(sub2), list(
        interactions=3L,
        files=1L,
        rows=c(5,3,3),
        cols=c(3,5,5)
    ))
    expect_identical(
        interactions(sub2),
        interactions(iarr)[3:1,]
    )

    ## Multiple subsets
    sub2 <- iarr[3:1,]
    sub3 <- sub2[3:1,]
    expect_identical(dim(sub3), dim(iarr[1:3,]))
    expect_identical(counts(sub3), counts(iarr)[,,1:3,])

    ## Returns InteractionArray
    sub <- iarr[1:2,] # return InteractionArray
    expect_s4_class(sub, "InteractionArray")
    expect_identical(
        as.matrix(counts(iarr)[,,1,1]),
        as.matrix(counts(sub)[,,1,1])
    )

    ## TODO build row/col names into
    ## InteractionJaggedArray
    expect_error(counts(iarr[1,], showDimnames=TRUE), NA)

})

test_that("InteractionJaggedArray subsetByOverlaps", {


    ## Shift first two ranges out of range
    gi2 <- c(assignToBins(gi[1:2], binSize=100e3, pos1=-200e3), gi[3:4])

    ## findOverlaps
    expect_identical(length(findOverlaps(iarr, iarr)), 4L)
    expect_identical(length(findOverlaps(iarr, gi2)), 2L)
    expect_identical(length(findOverlaps(iarr)), 4L)

    ## countOverlaps
    expect_identical(countOverlaps(iarr), rep(1L, 4))
    expect_identical(countOverlaps(iarr, gi), rep(1L, 4))
    expect_identical(
        countOverlaps(iarr, gi2),
        c(rep(0L, 2), rep(1L, 2))
    )

    ## overlapsAny
    expect_identical(
        overlapsAny(iarr, gi2),
        c(rep(FALSE, 2), rep(TRUE,2))
    )
    expect_identical(overlapsAny(iarr, iarr), rep(TRUE, 4L))
    expect_identical(overlapsAny(iarr), rep(TRUE, 4L))

    ## subsetByOverlaps
    expect_identical(length(subsetByOverlaps(iarr, gi2)), 2L)
    expect_identical(
        length(subsetByOverlaps(iarr, gi2, maxgap=100e3)),
        4L
    )
    expect_identical(
        subsetByOverlaps(iarr[,1], gi2) |> counts(),
        counts(iarr)[,,3:4,1]
    )
    expect_s4_class(
        subsetByOverlaps(iarr, gi2, invert=TRUE),
        "InteractionArray"
    )
    expect_identical(
        subsetByOverlaps(iarr, gi2, invert=TRUE),
        iarr[1:2]
    )

    ## Works with vector subclass (GRanges)
    expect_error(findOverlaps(iarr, first(gi2)), NA)
    expect_error(subsetByOverlaps(iarr, first(gi2)), NA)
    expect_error(countOverlaps(iarr, first(gi2)), NA)


})
