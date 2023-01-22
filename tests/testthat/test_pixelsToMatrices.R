library(mariner)

## Example GInteractions object
gi <- data.frame(
    chr1 = "chr1",
    start1 = c(100, 4000),
    end1 = c(125, 4025),
    chr2 = "chr1",
    start2 = c(1000, 6000),
    end2 = c(1025, 6025)
) |> as_ginteractions()


test_that("metadata columns are kept", {

    gi$col1 <- c(1,2)
    gi$col2 <- paste0("range", seq_along(gi))

    expect_identical(
        mcols(pixelsToMatrices(x=gi, buffer=3)),
        mcols(gi)
    )

})

test_that("throws error when ranges are not binned", {
    binPairs(gi, binSize=c(10, 20)) |>
        pixelsToMatrices(buffer=3) |>
        expect_error(".*must be the same width.*")
})
