

test_that("", {
    ## Write tests for aggHicMatrices
    aggHicMatrices(x = iarr, FUN = median, nBlocks = 5,
                   BPPARAM = BiocParallel::SerialParam())

    aggHicMatrices(x = iarr, by="files", FUN = median, nBlocks = 5,
                   BPPARAM = BiocParallel::SerialParam())

    aggHicMatrices(x = iarr, by="interactions", FUN = sum, nBlocks = 5,
                   BPPARAM = BiocParallel::SerialParam())


    aggHicMatrices(iarr, BPPARAM = BiocParallel::SerialParam())
    aggHicMatrices(iarr, by = "files", BPPARAM = BiocParallel::SerialParam())
    aggHicMatrices(iarr, by = "interactions")

    iarr[1:20, 1] |>
        aggHicMatrices(by="interactions", nBlocks = 1)

    iarr[1:20, 2] |>
        aggHicMatrices(by="interactions", nBlocks = 1, verbose = FALSE)

    iarr[1,] |>
        aggHicMatrices(by="interactions")


    identical(iarr[1,] |>
                  aggHicMatrices(by="files", verbose=FALSE) |>
                  as.array(),
              drop(counts(iarr[1,])) |>
                  as.array())
})
