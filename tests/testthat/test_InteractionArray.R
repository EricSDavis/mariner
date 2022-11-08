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
loops <- pixelsToMatrices(loops, buffer=5)


test_that("counts accessor for InteractionArray", {

    ## Extract 10, 11x11 count matrices from 2 hic files
    iarr <-
        loops[1:10] |>
        pullHicMatrices(binSize=5e03,
                        files=hicFiles)

    ref <- aperm(assay(iarr, 'counts'), c(3,4,1,2))

    expect_identical(counts(iarr, FALSE), ref)
    expect_snapshot(counts(iarr, TRUE))
    expect_snapshot(counts(iarr[1:3, 1:2]))
    expect_snapshot(counts(iarr[3:4, 1]))
    expect_snapshot(counts(iarr[1:7,1:2]))
    expect_snapshot(counts(iarr[1:7,1]))
    expect_snapshot(counts(iarr[1,1:2]))
    expect_snapshot(counts(iarr[1,1]))

    expect_error(counts(InteractionArray()),
                 ".*has no count matrices.$")
})


test_that("Concatenating InteractionArray objects", {

    ## Reference
    iarr <-
        loops[1:10] |>
        pullHicMatrices(binSize=5e03,
                        files=hicFiles)

    ## Different metadata
    iarr2 <-
        loops[1:10] |>
        pullHicMatrices(binSize=5e03,
                        files=hicFiles,
                        norm = "KR")
    ## Different colData
    iarr3 <-
        loops[1:10] |>
        pullHicMatrices(binSize=5e03,
                        files=hicFiles[1])

    ## rbind
    expect_s4_class(rbind(iarr, iarr), "InteractionArray")
    expect_equal(dim(rbind(iarr, iarr)), c(20, 2))
    expect_error(rbind(iarr, iarr2),
                 "Can't rbind.*Array.*with different metadata.$")
    expect_error(rbind(iarr, iarr3),
                 "Can't rbind.*with different colData.$")

    ## cbind
    expect_s4_class(cbind(iarr, iarr), "InteractionArray")
    expect_equal(dim(cbind(iarr, iarr)), c(10, 4))
    expect_equal(dim(cbind(iarr, iarr3)), c(10, 3))
    expect_error(cbind(iarr, iarr2),
                 "Can't cbind.*Array.*with different metadata.$")
    expect_error(cbind(iarr[1], iarr[2]),
                 "interactions must be identical in 'cbind'$")

})

test_that("show method for InteractionArray", {

    ## Extract 10, 11x11 count matrices from 2 hic files
    iarr <-
        loops[1:10] |>
        pullHicMatrices(binSize=5e03,
                        files=hicFiles)
    expect_snapshot(show(iarr))
    expect_snapshot(show(iarr[1:3,]))
    expect_snapshot(show(iarr[1:3,1]))
    expect_snapshot(show(iarr[,2]))

})
