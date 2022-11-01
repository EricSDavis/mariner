

test_that("", {

    ## Write tests for enhancePixels
    x |>
        binPairs(binSize = 5e03) |>
        pullHicMatrices(binSize = 10e03, files = hicFiles[1]) |>
        selectMaxPixel()

    x |>
        binPairs(binSize = 5e03) |>
        pullHicMatrices(binSize = 10e03, files = hicFiles[1])

    tmp <- enhancePixels(x, namedHicFiles, from=10e03, to=5e03)
    tmp <- enhancePixels(x, namedHicFiles, from=5e03, to=10e03)
    expect_identical(tmp |>
                         pullHicPixels(10e03, namedHicFiles) |>
                         counts() |>
                         rowSums(),
                     tmp$countSum)
})
