library(mariner)
library(data.table, include.only = "fread")

## Shared objects --------------------------------------------------------------


## Reference BEDPE files (loops called with SIP)
bedpeFiles <-
    system.file("extdata", package = "mariner") |>
    list.files(pattern = "Loops.txt", full.names = TRUE)

## Read in bedpeFiles as a list of GInteractions
giList <-
    lapply(bedpeFiles, fread) |>
    lapply(as_ginteractions)

## Merge pairs
mgi <- mergePairs(x = giList,
                  radius = 50e03)

## Test pullHicMatrices --------------------------------------------------------

## TODO: Change .pullHicMatrices to pullHicMatrices
test_that("pullHicMatrices checks for binning", {

    ## Inform if not binned
    .handleBinning(x = mgi, binSize = 50e03) |>
        expect_message("Pairs are not binned.*")

    ## Still returns binned object
    .handleBinning(x = mgi, binSize = 50e03) |>
        suppressMessages() |>
        expect_s4_class("MergedGInteractions")

    ## Accepts binned pairs
    .handleBinning(x = mgi |> binPairs(50e03),
                   binSize = 50e03) |>
        expect_s4_class("MergedGInteractions")

    ## Works for pullHicMatrices
    .pullHicMatrices(x = mgi, binSize = 50e03) |>
        expect_message("Pairs are not binned.*")

})

test_that("develop the function", {

    .pullHicMatrices(x = mgi, binSize = 50e03)

})
