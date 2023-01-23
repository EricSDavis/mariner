library(mariner)

## Shared objects --------------------------------------------------------------

## Read .hic file paths
hicFiles <-
    system.file("extdata/test_hic", package="mariner") |>
    list.files(pattern=".hic", full.names=TRUE)

## Create example interactions
x <- read.table(text="
        9 14000000 14500000 9 14500000 15000000
        9 89500000 90000000 9 89500000 90000000
        9 23500000 24000000 9 23500000 24000000")
x <- as_ginteractions(x)


test_that("counts accessor for InteractionMatrix", {

    ## Extract 3 pixels from 2 hic files
    imat <- pullHicPixels(x, hicFiles, 500e03)
    ref <- assay(imat, 'counts')

    expect_identical(counts(imat), ref)
    expect_error(counts(InteractionMatrix()),
                 ".*has no counts.$")
})

test_that("counts<- replace method for InteractionMatrix", {

    imat <- pullHicPixels(x, hicFiles, 500e3)
    counts(imat) <- as.matrix(counts(imat))
    expect_equal(class(counts(imat)), c("matrix", "array"))

})

test_that("Concatenating InteractionMatrix objects", {

    ## Reference
    imat <- pullHicPixels(x, hicFiles, 500e03)

    ## Different metadata
    imat2 <- pullHicPixels(x, hicFiles, 500e03, norm="KR")

    ## Different colData
    imat3 <- pullHicPixels(x, hicFiles[1], 500e03)

    ## rbind
    expect_s4_class(rbind(imat, imat), class(imat))
    expect_equal(dim(rbind(imat, imat)), c(6,2))
    expect_error(rbind(imat, imat2),
                 "Can't rbind.*Matrix.*with different metadata.$")
    expect_error(rbind(imat, imat3),
                 "Can't rbind.*Matrix.*with different colData.$")

    ## cbind
    expect_s4_class(cbind(imat, imat), class(imat))
    expect_equal(dim(cbind(imat, imat)), c(3, 4))
    expect_equal(dim(cbind(imat, imat3)), c(3, 3))
    expect_error(cbind(imat, imat2),
                 "Can't cbind.*Matrix.*with different metadata.$")
    expect_error(cbind(imat[1], imat[2]),
                 "interactions must be identical in 'cbind'$")

})

test_that("show method for InteractionArray", {

    ## Extract 3 pixels from 2 hic files
    imat <- pullHicPixels(x, hicFiles, 500e03)

    expect_snapshot(show(imat))
    expect_snapshot(show(imat[1:2,]))
    expect_snapshot(show(imat[1:2,1]))
    expect_snapshot(show(imat[,2]))

})
