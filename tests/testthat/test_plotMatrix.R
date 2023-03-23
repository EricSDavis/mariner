library(mariner)

test_that(".parseParams", {

    ## Embed param parsing function
    testFun <- function(params=NULL, x=1, y=NULL, z="a") {

        ## Get parsed arguments
        parsedArgs <- .parseParams(
            params = params,
            defaultArgs = formals(eval(match.call()[[1]])),
            declaredArgs = lapply(match.call()[-1], eval.parent, n=2)
        )

        ## Evaluate parsed arguments
        parsedArgs <- lapply(parsedArgs, eval)

        return(parsedArgs)
    }

    library(plotgardener)
    testFun() |>
        expect_identical(list(x = 1, y = NULL, z = "a"))
    testFun(params=pgParams(z="b")) |>
        expect_identical(list(x = 1, y = NULL, z = "b"))
    testFun(params=pgParams(y=2)) |>
        expect_identical(list(x = 1, y = 2, z = "a"))
    testFun(params=pgParams(y=2), y=3) |>
        expect_identical(list(x = 1, y = 3, z = "a"))
    testFun(params = "a") |>
        expect_error("`params` must be a `pgParams` object.")
})

test_that(".checkZrange", {
    .checkZrange(zrange=NULL) |> expect_null()
    .checkZrange(zrange=c(0,1)) |> expect_null()
    .checkZrange(zrange=1) |>
        expect_error(".*length 2 vector.")
    .checkZrange(zrange=c(1,1)) |>
        expect_error("`zrange[2]` must be > `zrange[1]`.", fixed=TRUE)
    .checkZrange(zrange=c("a", 1)) |>
        expect_error(".*vector of two numbers.")
})

test_that(".setZrange", {

    ## Positive sequential
    x <- structure(
        .Data=list(
            data = matrix(1:25, 5, 5),
            zrange = NULL
        ),
        class = "MatrixPlot"
    )
    expect_identical(.setZrange(x)$zrange, c(0, 25))

    ## Negative sequential
    x <- structure(
        .Data=list(
            data = matrix(-c(1:25), 5, 5),
            zrange = NULL
        ),
        class = "MatrixPlot"
    )
    expect_identical(.setZrange(x)$zrange, c(-25, 0))

    ## Divergent (center on 0 by default)
    x <- structure(
        .Data=list(
            data = matrix((-5):5, 11, 11),
            zrange = NULL
        ),
        class = "MatrixPlot"
    )
    expect_identical(.setZrange(x)$zrange, c(-5, 5))

})

test_that("plotMatrix", {

    ## Accepts matrix
    p <- plotMatrix(data=matrix(1:25, 5, 5), draw=FALSE)
    expect_s3_class(p, "MatrixPlot")

    ## Accepts DelayedMatrix
    library(DelayedArray)
    p <- plotMatrix(data=DelayedArray(matrix(1:25, 5, 5)), draw=FALSE)
    expect_s3_class(p, "MatrixPlot")

    ## Accepts converted matrix
    p <- DelayedArray(matrix(1:25, 5, 5)) |>
        as.matrix() |>
        plotMatrix(draw=FALSE)
    expect_s3_class(p, "MatrixPlot")

    ## Accepts non-square data
    p <- plotMatrix(data=matrix(1:24, 3, 8), draw=FALSE)
    expect_s3_class(p, "MatrixPlot")

})

test_that("plotMatrix accepts matrices with NA values", {

    m <- matrix(1:25, 5, 5)
    m[1,] <- NA
    p <- plotMatrix(data=m, na.color="grey", draw=FALSE)
    expect_identical(p$na.color, "grey")

})
