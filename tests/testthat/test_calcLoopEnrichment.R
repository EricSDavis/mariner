# library(mariner)

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


test_that("Getting binSize", {

    .getBinSize(x=loops) |>
        expect_identical(5000)

    .getBinSize(x=binPairs(loops, 10e03)) |>
        expect_identical(10000)

})

test_that("Able to make manhattan distance matrix", {

    .manhattanMatrix(-1) |>
        expect_error("`buffer` must be >= 0.")

    .manhattanMatrix(0) |>
        expect_identical(matrix(0))

    .manhattanMatrix(1) |>
        expect_identical(
            matrix(data=c(2L, 1L, 2L, 1L, 0L, 1L, 2L, 1L, 2L),
                   nrow=3, ncol=3))

    .manhattanMatrix(5) |>
        expect_identical(
            matrix(data=c(10, 9, 8, 7, 6, 5, 6, 7, 8, 9, 10,
                          9, 8, 7, 6, 5, 4, 5, 6, 7, 8, 9, 8,
                          7, 6, 5, 4, 3, 4, 5, 6, 7, 8, 7, 6,
                          5, 4, 3, 2, 3, 4, 5, 6, 7, 6, 5, 4,
                          3, 2, 1, 2, 3, 4, 5, 6, 5, 4, 3, 2,
                          1, 0, 1, 2, 3, 4, 5, 6, 5, 4, 3, 2,
                          1, 2, 3, 4, 5, 6, 7, 6, 5, 4, 3, 2,
                          3, 4, 5, 6, 7, 8, 7, 6, 5, 4, 3, 4,
                          5, 6, 7, 8, 9, 8, 7, 6, 5, 4, 5, 6,
                          7, 8, 9, 10, 9, 8, 7, 6, 5, 6, 7, 8,
                          9, 10) |> as.integer(),
                   nrow=11, ncol=11))

})

test_that(".selectRadius works", {

    .selectRadius(c(1,3,5), buffer=5, invert=FALSE) |>
        expect_identical(
            c(6L, 16L, 18L, 26L, 28L, 30L, 36L, 38L, 40L,
              42L, 46L, 48L, 50L, 52L, 54L, 56L, 58L, 60L,
              62L, 64L, 66L, 68L, 70L, 72L, 74L, 76L, 80L,
              82L, 84L, 86L, 92L, 94L, 96L, 104L, 106L, 116L)
        )

    .selectCenterPixel(mhDist=1:3, buffer=5, invert=FALSE) |>
        expect_identical(
            c(28L, 38L, 39L, 40L, 48L, 49L, 50L, 51L, 52L,
              58L, 59L, 60L, 61L, 62L, 63L, 64L, 70L, 71L,
              72L, 73L, 74L, 82L, 83L, 84L, 94L)
        )

})

test_that("angle works", {

    .angleMatrix(1) |>
        expect_equal(
            structure(c(135, 180, 225, 90, NaN, 270, 45, 0, 315),
                      dim = c(3L, 3L))
        )
    .angleMatrix(2) |>
        expect_equal(
            structure(c(135, 153.434948822922, 180,
                        206.565051177078, 225, 116.565051177078,
                        135, 180, 225, 243.434948822922, 90, 90,
                        NaN, 270, 270, 63.434948822922, 45, 0,
                        315, 296.565051177078, 45, 26.565051177078,
                        0, 333.434948822922, 315),
                      dim = c(5L, 5L))
        )

    .angleMatrix(3) |>
        expect_equal(
            structure(c(135, 146.30993247402, 161.565051177078,
                        180, 198.434948822922, 213.69006752598,
                        225, 123.69006752598, 135, 153.434948822922,
                        180, 206.565051177078, 225, 236.30993247402,
                        108.434948822922, 116.565051177078, 135, 180,
                        225, 243.434948822922, 251.565051177078,
                        90, 90, 90, NaN, 270, 270, 270, 71.565051177078,
                        63.434948822922, 45, 0, 315, 296.565051177078,
                        288.434948822922, 56.3099324740202, 45,
                        26.565051177078, 0, 333.434948822922, 315,
                        303.69006752598, 45, 33.6900675259798,
                        18.434948822922, 0, 341.565051177078,
                        326.30993247402, 315),
                      dim = c(7L, 7L))
        )

})


test_that(".selection helper functions work", {

    .selectSubmatrix(m = matrix(rep(c(1,0,1), 3), nrow=3, ncol=3),
                     invert=FALSE) |>
        .showSelection(buffer = 1)
    .selectSubmatrix(m = matrix(rep(c(1,0,1), 3), nrow=3, ncol=3),
                     invert=TRUE) |>
        .showSelection(buffer = 1)
    img <- read.table(text="
                   0 0 0 0 0 0 0 0 0
                   0 0 0 0 0 0 0 0 0
                   0 1 0 1 0 1 1 1 0
                   0 1 0 1 0 0 1 0 0
                   0 1 1 1 0 0 1 0 0
                   0 1 0 1 0 0 1 0 0
                   0 1 0 1 0 1 1 1 0
                   0 0 0 0 0 0 0 0 0
                   0 0 0 0 0 0 0 0 0
                   ")
    .selectSubmatrix(m=img, invert = FALSE) |>
        .showSelection(buffer= 4)
    .selectSubmatrix(m=img, invert = TRUE) |>
        .showSelection(buffer= 4)

    .selectCoordinates(buffer=5, rowInd=1:3, colInd=1:3, invert=FALSE) |>
        .showSelection(buffer=5)
    .selectCoordinates(buffer=5, rowInd=1:3, colInd=1:3, invert=TRUE) |>
        .showSelection(buffer=5)

    .selectCoordinates(buffer=5,
                       rowInd=c(6, 4), colInd=c(8, 5), invert=FALSE) |>
        .showSelection(buffer=5)
    .selectCoordinates(buffer=5,
                       rowInd=c(6, 4), colInd=c(8, 5), invert=TRUE) |>
        .showSelection(buffer=5)

    .selectBlock(buffer=5, rowInd=3:6, colInd=4:7, invert=FALSE) |>
        .showSelection(buffer=5)
    .selectBlock(buffer=5, rowInd=3:6, colInd=4:7, invert=TRUE) |>
        .showSelection(buffer=5)

    .selectBlock(buffer=5, rowInd=100, colInd=1:3, invert=FALSE) |>
        .showSelection(buffer=5) |>
        expect_error(".*subscript out of bounds")


    .selectTopLeft(buffer=5, n=3, invert=FALSE) |>
        .showSelection(buffer=5)
    .selectTopLeft(buffer=5, n=3, invert=TRUE) |>
        .showSelection(buffer=5)

    .selectTopRight(buffer=5, n=3, invert=FALSE) |>
        .showSelection(buffer=5)
    .selectTopRight(buffer=5, n=3, invert=TRUE) |>
        .showSelection(buffer=5)

    .selectBottomRight(buffer=5, n=3, invert=FALSE) |>
        .showSelection(buffer=5)
    .selectBottomRight(buffer=5, n=3, invert=TRUE) |>
        .showSelection(buffer=5)

    .selectBottomLeft(buffer=5, n=3, invert=FALSE) |>
        .showSelection(buffer=5)
    .selectBottomLeft(buffer=5, n=3, invert=TRUE) |>
        .showSelection(buffer=5)


    .selectCorners(buffer=5, n=3, invert=FALSE) |>
        .showSelection(buffer=5)
    .selectCorners(buffer=5, n=3, invert=TRUE) |>
        .showSelection(buffer=5)



    .selectInner(buffer=5, n=3, invert=FALSE) |>
        .showSelection(buffer=5)

    .selectInner(buffer=5, n=3, invert=TRUE) |>
        .showSelection(buffer=5)

    .selectRows(buffer=5, rows=1:3, invert=FALSE) |>
        .showSelection(buffer=5)
    .selectRows(buffer=5, rows=1:3, invert=TRUE) |>
        .showSelection(buffer=5)

    .selectCols(buffer=5, cols=1:3, invert=FALSE) |>
        .showSelection(buffer=5)
    .selectCols(buffer=5, cols=1:3, invert=TRUE) |>
        .showSelection(buffer=5)

    .selectOuter(buffer=5, n=1, invert=FALSE) |>
        .showSelection(buffer=5)
    .selectOuter(buffer=5, n=1, invert=TRUE) |>
        .showSelection(buffer=5)

    .selectOuter(buffer=5, n=3, invert=FALSE) |>
        .showSelection(buffer=5)
    .selectOuter(buffer=5, n=3, invert=TRUE) |>
        .showSelection(buffer=5)


    intersect(.selectOuter(buffer=5, n=1, invert=TRUE),
              .selectOuter(buffer=5, n=2, invert=FALSE)) |>
        .showSelection(buffer=5)

    .selectRadius(c(1,3,5), buffer=5, invert=FALSE) |>
        .showSelection(buffer=5)

    .angleMatrix(3)
    .manhattanMatrix(3)
    .selectCenterPixel(0, 5, FALSE)
    .selectCenterPixel(1:2, 5, FALSE) |>
        .showSelection(5)


    ## Argument for excluding center pixel (for all functions)
    ## TRUE for selectCenterPixel

    .selectRadius(1:3, 4, FALSE) |>
        .showSelection(4)

})


test_that("calcLoopEnrichment", {


    ## Doesn't accept loops of variable resolution
    loops |>
        binPairs(c(5e03, 10e03)) |>
        calcLoopEnrichment(files=hicFiles) |>
        expect_error("All ranges.*must be equal widths.")

    ## nBlocks must be > 0
    calcLoopEnrichment(loops, hicFiles, nBlocks=0) |>
        expect_error("`nBlocks` must be > 0.")

    res <-
        calcLoopEnrichment(x=loops[1:10] |> binPairs(100e03),
                           files=hicFiles)

    res


    .diagonalMatrix(loops[1:3], 1)

})


