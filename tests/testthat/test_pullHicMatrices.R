library(mariner)
library(data.table, include.only = "fread")

## Shared objects --------------------------------------------------------------

## Test .hic files
hicFiles <-
    system.file("extdata/test_hic", package = "mariner") |>
    list.files(full.names = TRUE)

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

## Bin MergedGInteractions
bgi <- binPairs(x = mgi, binSize = 50e03)

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
    .handleBinning(x = bgi,
                   binSize = 50e03) |>
        expect_s4_class("MergedGInteractions")

    ## Throws warning when using a different binSize
    .handleBinning(x = bgi,
                   binSize = 10e03) |>
        expect_message(".*Binning with binSize=10000.*")

    ## Works for pullHicMatrices
    .pullHicMatrices(x = mgi, binSize = 50e03) |>
        expect_message("Pairs are not binned.*")

})

test_that("straw returns data in expected order", {

    ## Load libraries
    library(glue, include.only = "glue")
    library(strawr, include.only = c("straw", "readHicChroms"))
    library(dplyr, include.only = "arrange")

    ## Define chromosome map & arrange by index
    chrMap <-
        readHicChroms(hicFiles[1]) |>
        arrange(index) |>
        head(n=25)

    ## Check every combination of chromosomes
    for (i in seq(2, nrow(chrMap))) {
        for (j in seq(2, nrow(chrMap))) {
            chr1loc <- glue('{chrMap[i, "name"]}:0:0')
            chr2loc <- chrMap[j, "name"]
            sparseMat <-
                straw(norm = "NONE",
                      fname = hicFiles[1],
                      chr1loc = chr1loc,
                      chr2loc = chr2loc,
                      unit = "BP",
                      binsize = 2500000,
                      matrix = "observed")

            ## All values in column "x" should be 0
            if (nrow(sparseMat) != 0) {
                if (i <= j) {
                    expect_true(all(sparseMat$x == 0),
                                label = paste(i,j))
                } else {
                    expect_false(all(sparseMat$x == 0),
                                 label = paste(i,j))
                }
            }
        }
    }

})

test_that("develop the function", {

    .pullHicMatrices(x = mgi, binSize = 50e03)

})
