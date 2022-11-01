

test_that("counts accessor", {
    ## Access count matrices

    x <- iarr[1:3,1:2]

    counts(iarr, showDimnames=TRUE)
    counts(iarr, showDimnames=FALSE)
    counts(iarr, TRUE)
    counts(iarr[1:3, 1:2])
    counts(iarr[3:4, 1])
    counts(iarr[1:7,1:2])
    counts(iarr[1:7,1])
    counts(iarr[1,1:2])
    counts(iarr[1,1])
})
