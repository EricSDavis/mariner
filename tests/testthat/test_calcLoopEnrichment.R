library(mariner)
library(marinerData)

## Shared objects --------------------------------------------------------------
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


test_that("selection helper functions work", {

    ## Select radius (manhattan distance) --------------------------------------

    exp <- c(6L, 16L, 18L, 26L, 28L, 30L, 36L, 38L, 40L, 42L, 46L, 48L,
             50L, 52L, 54L, 56L, 58L, 60L, 62L, 64L, 66L, 68L, 70L, 72L, 74L,
             76L, 80L, 82L, 84L, 86L, 92L, 94L, 96L, 104L, 106L, 116L)
    .selectRadius(c(1,3,5), buffer=5, invert=FALSE) |>
        expect_identical(exp)
    selectRadius(c(1,3,5), buffer=5)$x |>
        expect_identical(exp)
    exp <- c(1L, 2L, 3L, 4L, 5L, 7L, 8L, 9L, 10L, 11L, 12L, 13L, 14L, 15L,
             17L, 19L, 20L, 21L, 22L, 23L, 24L, 25L, 27L, 29L, 31L, 32L, 33L,
             34L, 35L, 37L, 39L, 41L, 43L, 44L, 45L, 47L, 49L, 51L, 53L, 55L,
             57L, 59L, 61L, 63L, 65L, 67L, 69L, 71L, 73L, 75L, 77L, 78L, 79L,
             81L, 83L, 85L, 87L, 88L, 89L, 90L, 91L, 93L, 95L, 97L, 98L, 99L,
             100L, 101L, 102L, 103L, 105L, 107L, 108L, 109L, 110L, 111L, 112L,
             113L, 114L, 115L, 117L, 118L, 119L, 120L, 121L)
    selectRadius(c(1,3,5), buffer=5, invert=TRUE)$x |>
        expect_identical(exp)

    ## Select center pixel -----------------------------------------------------

    exp <- c(39L, 49L, 50L, 51L, 59L, 60L, 61L, 62L, 63L, 71L, 72L, 73L, 83L)
    selectCenterPixel(1:2, 5)$x |>
        expect_identical(exp)
    .selectCenterPixel(0, 5, FALSE) |>
        expect_identical(61L)
    selectCenterPixel(0, 5)$x |>
        expect_identical(61L)

    ## Select submatrix --------------------------------------------------------

    exp <- c(1L, 3L, 4L, 6L, 7L, 9L)
    .selectSubmatrix(m = matrix(rep(c(1,0,1), 3), nrow=3, ncol=3),
                     invert=FALSE) |>
        expect_identical(exp)
    selectSubmatrix(m = matrix(rep(c(1,0,1), 3), nrow=3, ncol=3)) |>
        expect_identical(exp)
    exp <- c(2L, 5L, 8L)
    selectSubmatrix(m = matrix(rep(c(1,0,1), 3), nrow=3, ncol=3),
                    invert=TRUE) |>
        expect_identical(exp)
    img <- as.matrix(read.table(text="
                   0 0 0 0 0 0 0 0 0
                   0 0 0 0 0 0 0 0 0
                   0 1 0 1 0 1 1 1 0
                   0 1 0 1 0 0 1 0 0
                   0 1 1 1 0 0 1 0 0
                   0 1 0 1 0 0 1 0 0
                   0 1 0 1 0 1 1 1 0
                   0 0 0 0 0 0 0 0 0
                   0 0 0 0 0 0 0 0 0
                   "))
    exp <- c(12L, 13L, 14L, 15L, 16L, 23L, 30L, 31L, 32L, 33L, 34L, 48L,
             52L, 57L, 58L, 59L, 60L, 61L, 66L, 70L)
    selectSubmatrix(m=img) |> expect_identical(exp)
    exp <- c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 11L, 17L, 18L, 19L,
             20L, 21L, 22L, 24L, 25L, 26L, 27L, 28L, 29L, 35L, 36L, 37L, 38L,
             39L, 40L, 41L, 42L, 43L, 44L, 45L, 46L, 47L, 49L, 50L, 51L, 53L,
             54L, 55L, 56L, 62L, 63L, 64L, 65L, 67L, 68L, 69L, 71L, 72L, 73L,
             74L, 75L, 76L, 77L, 78L, 79L, 80L, 81L)
    selectSubmatrix(m=img, invert=TRUE) |>
        expect_identical(exp)

    ## Select coordinates ------------------------------------------------------

    exp <- c(1L, 13L, 25L)
    .selectCoordinates(buffer=5, rowInd=1:3, colInd=1:3, invert=FALSE) |>
        expect_identical(exp)
    selectCoordinates(buffer=5, rowInd=1:3, colInd=1:3)$x |>
        expect_identical(exp)

    exp <- c(2L, 3L, 4L, 6L, 7L, 8L)
    .selectCoordinates(buffer=1, rowInd=1:3, colInd=1:3, invert=TRUE) |>
        expect_identical(exp)
    selectCoordinates(buffer=1, rowInd=1:3, colInd=1:3, invert=TRUE)$x |>
        expect_identical(exp)

    exp <- c(8L)
    .selectCoordinates(buffer=1, rowInd=2, colInd=3, invert=FALSE) |>
        expect_identical(exp)
    selectCoordinates(buffer=1, rowInd=2, colInd=3)$x |>
        expect_identical(exp)

    ## Select block ------------------------------------------------------------

    exp <- c(1L, 2L, 3L, 6L, 7L, 8L, 11L, 12L, 13L)
    .selectBlock(buffer=2, rowInd=1:3, colInd=1:3, invert=FALSE) |>
        expect_identical(exp)
    selectBlock(buffer=2, rowInd=1:3, colInd=1:3)$x |>
        expect_identical(exp)

    exp <- c(5L, 10L, 15L, 20L, 21L, 22L, 23L, 24L, 25L)
    .selectBlock(buffer=2, rowInd=1:4, colInd=1:4, invert=TRUE) |>
        expect_identical(exp)
    selectBlock(buffer=2, rowInd=1:4, colInd=1:4, invert=TRUE)$x |>
        expect_identical(exp)

    exp <- c(12L, 13L, 17L, 18L, 22L, 23L)
    .selectBlock(buffer=2, rowInd=2:3, colInd=3:5, invert=FALSE) |>
        expect_identical(exp)
    selectBlock(buffer=2, rowInd=2:3, colInd=3:5)$x |>
        expect_identical(exp)

    exp <- c(1L, 2L, 3L, 4L, 5L, 6L, 10L, 11L, 15L, 16L, 20L, 21L, 25L)
    .selectBlock(buffer=2, rowInd=2:4, colInd=2:5, invert=TRUE) |>
        expect_identical(exp)
    selectBlock(buffer=2, rowInd=2:4, colInd=2:5, invert=TRUE)$x |>
        expect_identical(exp)

    selectBlock(buffer=5, rowInd=100, colInd=1:3) |>
        expect_error(".*subscript out of bounds")


    ## Select TopLeft ----------------------------------------------------------

    exp <- c(1L, 2L, 3L, 12L, 13L, 14L, 23L, 24L, 25L)
    .selectTopLeft(buffer=5, n=3, inset=0, invert=FALSE) |>
        expect_identical(exp)
    selectTopLeft(buffer=5, n=3)$x |>
        expect_identical(exp)

    exp <- c(13L, 14L, 15L, 24L, 25L, 26L, 35L, 36L, 37L)
    .selectTopLeft(buffer=5, n=3, inset=1, invert=FALSE) |>
        expect_identical(exp)
    selectTopLeft(buffer=5, n=3, inset=1)$x |>
        expect_identical(exp)

    exp <- c(1L, 2L, 12L, 13L)
    .selectTopLeft(buffer=5, n=3, inset=-1, invert=FALSE) |>
        expect_identical(exp)
    selectTopLeft(buffer=5, n=3, inset=-1)$x |>
        expect_identical(exp)

    exp <- c(5L, 10L, 15L, 20L, 21L, 22L, 23L, 24L, 25L)
    .selectTopLeft(buffer=2, n=4, inset=0, invert=TRUE) |>
        expect_identical(exp)
    selectTopLeft(buffer=2, n=4, invert=TRUE)$x |>
        expect_identical(exp)

    exp <- c(1L, 2L, 3L, 4L, 5L, 6L, 11L, 16L, 21L)
    .selectTopLeft(buffer=2, n=4, inset=1, invert=TRUE) |>
        expect_identical(exp)
    selectTopLeft(buffer=2, n=4, inset=1, invert=TRUE)$x |>
        expect_identical(exp)

    selectTopLeft(buffer=5, n=3, inset=100) |>
        expect_error(".*subscript out of bounds")

    ## Select TopRight ---------------------------------------------------------

    exp <- c(89L, 90L, 91L, 100L, 101L, 102L, 111L, 112L, 113L)
    .selectTopRight(buffer=5, n=3, inset=0, invert=FALSE) |>
        expect_identical(exp)
    selectTopRight(buffer=5, n=3)$x |>
        expect_identical(exp)

    exp <- c(79L, 80L, 81L, 90L, 91L, 92L, 101L, 102L, 103L)
    .selectTopRight(buffer=5, n=3, inset=1, invert=FALSE) |>
        expect_identical(exp)
    selectTopRight(buffer=5, n=3, inset=1)$x |>
        expect_identical(exp)

    exp <- c(1L, 2L, 3L, 4L, 5L, 10L, 15L, 20L, 25L)
    .selectTopRight(buffer=2, n=4, inset=0, invert=TRUE) |>
        expect_identical(exp)
    selectTopRight(buffer=2, n=4, invert=TRUE)$x |>
        expect_identical(exp)

    exp <- c(1L, 6L, 11L, 16L, 21L, 22L, 23L, 24L, 25L)
    .selectTopRight(buffer=2, n=4, inset=1, invert=TRUE) |>
        expect_identical(exp)
    selectTopRight(buffer=2, n=4, inset=1, invert=TRUE)$x |>
        expect_identical(exp)

    selectTopRight(buffer=5, n=3, inset=100) |>
        expect_error(".*subscript out of bounds")

    ## Select BottomRight ------------------------------------------------------

    exp <- c(97L, 98L, 99L, 108L, 109L, 110L, 119L, 120L, 121L)
    .selectBottomRight(buffer=5, n=3, inset=0, invert=FALSE) |>
        expect_identical(exp)
    selectBottomRight(buffer=5, n=3)$x |>
        expect_identical(exp)

    exp <- c(85L, 86L, 87L, 96L, 97L, 98L, 107L, 108L, 109L)
    .selectBottomRight(buffer=5, n=3, inset=1, invert=FALSE) |>
        expect_identical(exp)
    selectBottomRight(buffer=5, n=3, inset=1)$x |>
        expect_identical(exp)

    exp <- c(1L, 2L, 3L, 4L, 5L, 6L, 11L, 16L, 21L)
    .selectBottomRight(buffer=2, n=4, inset=0, invert=TRUE) |>
        expect_identical(exp)
    selectBottomRight(buffer=2, n=4, invert=TRUE)$x |>
        expect_identical(exp)

    exp <- c(5L, 10L, 15L, 20L, 21L, 22L, 23L, 24L, 25L)
    .selectBottomRight(buffer=2, n=4, inset=1, invert=TRUE) |>
        expect_identical(exp)
    selectBottomRight(buffer=2, n=4, inset=1, invert=TRUE)$x |>
        expect_identical(exp)

    selectBottomRight(buffer=5, n=3, inset=-1) |>
        expect_error(".*subscript out of bounds")

    ## Select BottomLeft -------------------------------------------------------

    exp <- c(9L, 10L, 11L, 20L, 21L, 22L, 31L, 32L, 33L)
    .selectBottomLeft(buffer=5, n=3, inset=0, invert=FALSE) |>
        expect_identical(exp)
    selectBottomLeft(buffer=5, n=3)$x |>
        expect_identical(exp)

    exp <- c(19L, 20L, 21L, 30L, 31L, 32L, 41L, 42L, 43L)
    .selectBottomLeft(buffer=5, n=3, inset=1, invert=FALSE) |>
        expect_identical(exp)
    selectBottomLeft(buffer=5, n=3, inset=1)$x |>
        expect_identical(exp)

    exp <- c(1L, 6L, 11L, 16L, 21L, 22L, 23L, 24L, 25L)
    .selectBottomLeft(buffer=2, n=4, inset=0, invert=TRUE) |>
        expect_identical(exp)
    selectBottomLeft(buffer=2, n=4, invert=TRUE)$x |>
        expect_identical(exp)

    exp <- c(1L, 2L, 3L, 4L, 5L, 10L, 15L, 20L, 25L)
    .selectBottomLeft(buffer=2, n=4, inset=1, invert=TRUE) |>
        expect_identical(exp)
    selectBottomLeft(buffer=2, n=4, inset=1, invert=TRUE)$x |>
        expect_identical(exp)

    selectBottomLeft(buffer=5, n=3, inset=-1) |>
        expect_error(".*subscript out of bounds")
    selectBottomLeft(buffer=5, n=3, inset=100) |>
        expect_error(".*subscript out of bounds")


    ## Select corners ----------------------------------------------------------

    exp <- c(1L, 2L, 3L, 12L, 13L, 14L, 23L, 24L, 25L, 89L, 90L, 91L, 100L,
             101L, 102L, 111L, 112L, 113L, 97L, 98L, 99L, 108L, 109L, 110L,
             119L, 120L, 121L, 9L, 10L, 11L, 20L, 21L, 22L, 31L, 32L, 33L)
    .selectCorners(buffer=5, n=3, inset=0, invert=FALSE) |>
        expect_identical(exp)
    selectCorners(buffer=5, n=3)$x |>
        expect_identical(exp)

    exp <- c(13L, 14L, 15L, 24L, 25L, 26L, 35L, 36L, 37L, 79L, 80L, 81L,
             90L, 91L, 92L, 101L, 102L, 103L, 85L, 86L, 87L, 96L, 97L, 98L,
             107L, 108L, 109L, 19L, 20L, 21L, 30L, 31L, 32L, 41L, 42L, 43L)
    .selectCorners(buffer=5, n=3, inset=1, invert=FALSE) |>
        expect_identical(exp)
    selectCorners(buffer=5, n=3, inset=1)$x |>
        expect_identical(exp)

    exp <- c(3L, 8L, 11L, 12L, 13L, 14L, 15L, 18L, 23L)
    .selectCorners(buffer=2, n=2, inset=0, invert=TRUE) |>
        expect_identical(exp)
    selectCorners(buffer=2, n=2, invert=TRUE)$x |>
        expect_identical(exp)

    exp <- c(1L, 2L, 3L, 4L, 5L, 6L, 10L, 11L, 15L, 16L, 20L, 21L, 22L,
             23L, 24L, 25L)
    .selectCorners(buffer=2, n=2, inset=1, invert=TRUE) |>
        expect_identical(exp)
    selectCorners(buffer=2, n=2, inset=1, invert=TRUE)$x |>
        expect_identical(exp)

    selectCorners(buffer=5, n=3, inset=-1) |>
        expect_error(".*subscript out of bounds")
    selectCorners(buffer=5, n=3, inset=100) |>
        expect_error(".*subscript out of bounds")

    ## Select rows -------------------------------------------------------------

    exp <- c(1L, 6L, 11L, 16L, 21L)
    .selectRows(buffer=2, rows=1, invert=FALSE) |>
        expect_identical(exp)
    selectRows(buffer=2, rows=1)$x |>
        expect_identical(exp)

    exp <- c(1L, 3L, 6L, 8L, 11L, 13L, 16L, 18L, 21L, 23L)
    .selectRows(buffer=2, rows=c(1,3), invert=FALSE) |>
        expect_identical(exp)
    selectRows(buffer=2, rows=c(1,3))$x |>
        expect_identical(exp)

    exp <- c(1L, 2L, 3L, 6L, 7L, 8L, 11L, 12L, 13L, 16L, 17L, 18L, 21L,
             22L, 23L)
    .selectRows(buffer=2, rows=c(1:3), invert=FALSE) |>
        expect_identical(exp)
    selectRows(buffer=2, rows=c(1:3))$x |>
        expect_identical(exp)

    exp <- c(2L, 3L, 4L, 5L, 7L, 8L, 9L, 10L, 12L, 13L, 14L, 15L, 17L, 18L,
             19L, 20L, 22L, 23L, 24L, 25L)
    .selectRows(buffer=2, rows=1, invert=TRUE) |>
        expect_identical(exp)
    selectRows(buffer=2, rows=1, invert=TRUE)$x |>
        expect_identical(exp)

    exp <- c(2L, 4L, 5L, 7L, 9L, 10L, 12L, 14L, 15L, 17L, 19L, 20L, 22L,
             24L, 25L)
    .selectRows(buffer=2, rows=c(1,3), invert=TRUE) |>
        expect_identical(exp)
    selectRows(buffer=2, rows=c(1,3), invert=TRUE)$x |>
        expect_identical(exp)

    exp <- c(4L, 5L, 9L, 10L, 14L, 15L, 19L, 20L, 24L, 25L)
    .selectRows(buffer=2, rows=c(1:3), invert=TRUE) |>
        expect_identical(exp)
    selectRows(buffer=2, rows=c(1:3), invert=TRUE)$x |>
        expect_identical(exp)


    ## Select cols -------------------------------------------------------------

    exp <- 1:5L
    .selectCols(buffer=2, cols=1, invert=FALSE) |>
        expect_identical(exp)
    selectCols(buffer=2, cols=1)$x |>
        expect_identical(exp)

    exp <- c(1L, 2L, 3L, 4L, 5L, 11L, 12L, 13L, 14L, 15L)
    .selectCols(buffer=2, cols=c(1,3), invert=FALSE) |>
        expect_identical(exp)
    selectCols(buffer=2, cols=c(1,3))$x |>
        expect_identical(exp)

    exp <- 1:15L
    .selectCols(buffer=2, cols=c(1:3), invert=FALSE) |>
        expect_identical(exp)
    selectCols(buffer=2, cols=c(1:3))$x |>
        expect_identical(exp)

    exp <- 6:25L
    .selectCols(buffer=2, cols=1, invert=TRUE) |>
        expect_identical(exp)
    selectCols(buffer=2, cols=1, invert=TRUE)$x |>
        expect_identical(exp)

    exp <- c(6L, 7L, 8L, 9L, 10L, 16L, 17L, 18L, 19L, 20L, 21L, 22L, 23L,
             24L, 25L)
    .selectCols(buffer=2, cols=c(1,3), invert=TRUE) |>
        expect_identical(exp)
    selectCols(buffer=2, cols=c(1,3), invert=TRUE)$x |>
        expect_identical(exp)

    exp <- 16:25L
    .selectCols(buffer=2, cols=c(1:3), invert=TRUE) |>
        expect_identical(exp)
    selectCols(buffer=2, cols=c(1:3), invert=TRUE)$x |>
        expect_identical(exp)


    ## Select inner ------------------------------------------------------------


    exp <- c(17L, 18L, 19L, 24L, 25L, 26L, 31L, 32L, 33L)
    .selectInner(buffer=3, n=1, invert=FALSE) |>
        expect_identical(exp)
    selectInner(buffer=3, n=1)$x |>
        expect_identical(exp)

    exp <- c(9L, 10L, 11L, 12L, 13L, 16L, 17L, 18L, 19L, 20L, 23L, 24L,
             25L, 26L, 27L, 30L, 31L, 32L, 33L, 34L, 37L, 38L, 39L, 40L, 41L)
    .selectInner(buffer=3, n=2, invert=FALSE) |>
        expect_identical(exp)
    selectInner(buffer=3, n=2)$x |>
        expect_identical(exp)

    exp <- c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 11L, 12L, 13L, 14L,
             15L, 16L, 20L, 21L, 22L, 23L, 27L, 28L, 29L, 30L, 34L, 35L, 36L,
             37L, 38L, 39L, 40L, 41L, 42L, 43L, 44L, 45L, 46L, 47L, 48L, 49L)
    .selectInner(buffer=3, n=1, invert=TRUE) |>
        expect_identical(exp)
    selectInner(buffer=3, n=1, TRUE)$x |>
        expect_identical(exp)

    exp <- c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 14L, 15L, 21L, 22L, 28L, 29L,
             35L, 36L, 42L, 43L, 44L, 45L, 46L, 47L, 48L, 49L)
    .selectInner(buffer=3, n=2, invert=TRUE) |>
        expect_identical(exp)
    selectInner(buffer=3, n=2, TRUE)$x |>
        expect_identical(exp)


    ## Select outer ------------------------------------------------------------

    exp <- c(1L, 8L, 15L, 22L, 29L, 36L, 43L, 7L, 14L, 21L, 28L, 35L, 42L,
             49L, 2L, 3L, 4L, 5L, 6L, 44L, 45L, 46L, 47L, 48L)
    .selectOuter(buffer=3, n=1, invert=FALSE) |>
        expect_identical(exp)
    selectOuter(buffer=3, n=1)$x |>
        expect_identical(exp)

    exp <- c(1L, 2L, 8L, 9L, 15L, 16L, 22L, 23L, 29L, 30L, 36L, 37L, 43L,
             44L, 7L, 6L, 14L, 13L, 21L, 20L, 28L, 27L, 35L, 34L, 42L, 41L,
             49L, 48L, 3L, 4L, 5L, 10L, 11L, 12L, 45L, 46L, 47L, 38L, 39L,
             40L)
    .selectOuter(buffer=3, n=2, invert=FALSE) |>
        expect_identical(exp)
    selectOuter(buffer=3, n=2)$x |>
        expect_identical(exp)

    exp <- c(17L, 18L, 19L, 24L, 25L, 26L, 31L, 32L, 33L)
    .selectOuter(buffer=3, n=2, invert=TRUE) |>
        expect_identical(exp)
    selectOuter(buffer=3, n=2, invert=TRUE)$x |>
        expect_identical(exp)

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

    ## Show multiselection
    .showMultiSelection(fg=selectRadius(0:1, 3)$x,
                        bg=selectRadius(0:1, 3, TRUE)$x,
                        buffer=3) |>
        expect_error(NA)


    ## Test whole function
    exp <- new("DelayedMatrix",
               seed = structure(
                   c(1, 1, 1.57894736842105, 1.33333333333333,
                     0.736842105263158, 1.33333333333333, 0.75,
                     1.76470588235294, 0.818181818181818,
                     0.333333333333333, 1.66666666666667, 1.4,
                     2.16666666666667, 1, 0.777777777777778, 1,
                     1, 1.44827586206897, 0.909090909090909, 0.8),
                   dim = c(10L, 2L),
                   dimnames = list(NULL, c("FS","WT"))
               ))
    res <- calcLoopEnrichment(x=loops[1:10] |> binPairs(100e03),
                              files=hicFiles)
    expect_equal(res, exp)
    expect_s4_class(res, "DelayedMatrix")
    expect_equal(dim(res), dim(exp))


    ## Selection concatenation
    (selectInner(1, 4) - selectCenterPixel(0, 4))$x |>
        expect_identical(c(31L, 32L, 33L, 40L, 42L, 49L, 50L, 51L))


    exp <- c(18L, 24L, 26L, 32L, 11L, 17L, 19L, 23L, 27L, 31L, 33L, 39L)
    (selectRadius(x=1, buffer=3) + selectRadius(x=2, buffer=3))$x |>
        expect_identical(exp)

    exp <- c(11L, 17L, 19L, 23L, 27L, 31L, 33L, 39L)
    (selectRadius(x=0:2, buffer=3) - selectRadius(x=0:1, buffer=3))$x |>
        expect_identical(exp)

    (selectRadius(x=1, buffer=3) + selectRadius(x=2, buffer=4)) |>
        expect_error("`buffer` must be the same")

})

test_that("InteractionArray method for calcLoopEnrichment runs correctly",{
 
   ## Test non-square matrices 
  mats <- 
    binPairs(loops[1:100],binSize = c(100e3, 50e3)) |>
    pullHicMatrices(
      files=hicFiles,
      binSize=10e3,
      half="upper",
      norm = "KR"
    )
  
  expect_error(calcLoopEnrichment(mats), "must have square")
  
  ## Test even matrices
  mats <- 
    binPairs(loops[1:100],binSize = 100e3) |>
    pullHicMatrices(
      files=hicFiles,
      binSize=50e3,
      half="upper",
      norm = "KR"
    )
  
  expect_error(calcLoopEnrichment(mats), "must be odd")
  
  ## Test non-matching buffers
  mats <- 
    binPairs(loops[1:100],100e3) |>
    pixelsToMatrices(buffer=5) |>
    pullHicMatrices(
      files=hicFiles,
      binSize=100e3,
      half="upper",
      norm = "KR"
    )
  
  calcLoopEnrichment(
    x = mats,
    fg = selectCenterPixel(mhDist=1, buffer = 10)
  ) |> expect_error("`buffer` must be the same")
  
  mats <- 
    binPairs(loops[1:100],100e3) |>
    pixelsToMatrices(buffer=10) |>
    pullHicMatrices(
      files=hicFiles,
      binSize=100e3,
      half="upper",
      norm = "KR"
    )
  
  calcLoopEnrichment(
    x = mats,
    fg = selectCenterPixel(mhDist=1, buffer = 5)
  ) |> expect_error("`buffer` must be the same")
  
  ## Test auto-matching buffers
  ## no foreground
  expect_s4_class(calcLoopEnrichment(x = mats,
                                    bg = selectTopLeft(n=4, buffer=10) +
                                      selectBottomRight(n=4, buffer=10)),
                  "DelayedMatrix")
  
  ## no background
  expect_s4_class(
    calcLoopEnrichment(x = mats, 
                       fg = selectCenterPixel(mhDist=1, buffer = 10)),
    "DelayedMatrix")
  
  ## fg is not a selection
  expect_error(calcLoopEnrichment(x = mats, fg = 20),
               "`fg` must be a MatrixSelection")
  
  ## bg is not a selection
  expect_error(calcLoopEnrichment(x = mats, bg = "hello"),
               "`bg` must be a MatrixSelection")
  
  ## no foreground & no background
  expect_s4_class(calcLoopEnrichment(mats), "DelayedMatrix")
  
  ## Test whole function
  exp <- new("DelayedMatrix",
             seed = structure(
               c(1, 1, 1.57894736842105, 1.33333333333333,
                 0.736842105263158, 1.33333333333333, 0.75,
                 1.76470588235294, 0.818181818181818,
                 0.333333333333333, 1.66666666666667, 1.4,
                 2.16666666666667, 1, 0.777777777777778, 1,
                 1, 1.44827586206897, 0.909090909090909, 0.8),
               dim = c(10L, 2L),
               dimnames = list(NULL, c("FS","WT"))
             ))
  
  mats <- 
    binPairs(loops[1:10],100e3) |>
    pixelsToMatrices(buffer=5) |>
    pullHicMatrices(
      files=hicFiles,
      binSize=100e3
    )

  res <- calcLoopEnrichment(mats)
  expect_equal(res, exp)
  expect_s4_class(res, "DelayedMatrix")
  expect_equal(dim(res), dim(exp))
})

test_that("calcLoopEnrichment function can use local environment objects",{
    
  simpleFun <- function(x, y, z=5) { x + y + z }
  ## test defaults
  testFun1 <- deparse1(.modifyEnrichFun(simpleFun, eframe=parent.frame()))
  expFun1 <- deparse1(function(x1, x, y, z) {x1[x] + x1[y] + 5})
  expect_equal(testFun1, expFun1)
  
  ## test custom environment
  en <- new_environment(data=list(z=15))
  testFun2 <- deparse1(.modifyEnrichFun(simpleFun, eframe=en))
  expFun2 <- deparse1(function(x1, x, y, z) {x1[x] + x1[y] + 15})
  expect_equal(testFun2, expFun2)
  
  ## test current environment
  z <- 10
  testFun3 <- deparse1(.modifyEnrichFun(simpleFun, eframe=parent.frame()))
  expFun3 <- deparse1(function(x1, x, y, z) {x1[x] + x1[y] + 10})
  expect_equal(testFun3, expFun3)
  
  ## test within a function
  testModifyFun <- function(f){
    z <- 35
    .modifyEnrichFun(f, environment()) # use environment instead of parent frame
  }
  testFun4 <- deparse1(testModifyFun(simpleFun))
  expFun4 <- deparse1(function(x1, x, y, z) {x1[x] + x1[y] + 35})
  expect_equal(testFun4, expFun4)
  
  ## test with a non-function
  expect_error(.modifyEnrichFun(12), "`FUN` must be a function")
  
  ## using a function name as a function argument
  f <- function(x, y, sum) { sum(x,y) + sum }
  expect_error(.modifyEnrichFun(f),"share name.*with a function")
  
  ## no default and not in environment
  simpleFun <- function(x, y, z, w) { x + y + z }
  rm(z)
  expect_error(.modifyEnrichFun(simpleFun), "object.*not found")
  
  mats <- 
    binPairs(loops[1:10],100e3) |>
    pixelsToMatrices(buffer=5) |>
    pullHicMatrices(
      files=hicFiles,
      binSize=100e3
    )
  
  ## User defined function that returns more than one value
  customFun <- function(fg, bg, pseudo=10) {
    (fg + pseudo) / (bg + pseudo)
  }
  expect_warning(calcLoopEnrichment(mats, FUN = customFun, verbose = F),
                 "return more than one value")
  
  ## User defined function with default
  scoringFunction <- function(fg,bg, pseudo=10){
    median((fg + 1 + pseudo)) / median((bg + 1 + pseudo))
  }
  
  scores10 <- calcLoopEnrichment(mats, FUN = scoringFunction)
  expect_s4_class(scores10, "DelayedMatrix")
  expect_equal(dim(mats), dim(scores10))
  
  ## overwrite default value for pseudo with local variable
  pseudo <- 25
  scores25 <- calcLoopEnrichment(mats, FUN = scoringFunction)
  expect_s4_class(scores25, "DelayedMatrix")
  expect_equal(dim(mats), dim(scores25))
  
  ## User defined function that calls calcLoopEnrichment
  scoreLoops <- function(pseudo = 45){
    calcLoopEnrichment(mats, FUN = scoringFunction, verbose = F)
  }
  
  scoresDefault <- scoreLoops() # using default  pseudo <- 10
  expect_s4_class(scoresDefault, "DelayedMatrix")
  expect_equal(dim(mats), dim(scoresDefault))
  
  scoresLocal <- scoreLoops() # using local variable to overwrite default
  expect_s4_class(scoresLocal, "DelayedMatrix")
  expect_equal(dim(mats), dim(scoresLocal))
  
  # setting pseudo for function overwrites local variable and default
  scoresSet <- scoreLoops(pseudo = 50) 
  expect_s4_class(scoresSet, "DelayedMatrix")
  expect_equal(dim(mats), dim(scoresSet))
})


