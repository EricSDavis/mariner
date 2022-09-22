library(mariner)
library(data.table, include.only = "fread")
library(InteractionSet, include.only = c("regions", "GInteractions"))
library(glue, include.only = "glue")
library(strawr, include.only = c("straw", "readHicChroms"))
library(dplyr, include.only = "arrange")
library(GenomeInfoDb)
library(GenomicRanges)

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

## Create out-of-ordered, interchromosomal ranges
## (for testing)
giList <-
    lapply(giList, \(x) {
        ## Create a sample of regions to permute
        set.seed(123)
        s <- sample(x = seq_len(length(regions(x))),
                    size = length(regions(x)),
                    replace = FALSE)

        ## Scramble regions
        regions(x) <- regions(x)[s]
        x
    })

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

})

test_that("straw returns data in expected order", {

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
                strawr::straw(norm = "NONE",
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

test_that("Check chromosomes in .hic file", {

    ## Assign to x (to avoid modifying in place)
    x <- bgi

    ## Seqnames mismatch throws error
    .checkHicChroms(x, hicFiles[1]) |>
        expect_error("There's.*craziness.*")

    ## Corrected seqnames don't throw error
    seqlevelsStyle(x) <- "ENSEMBL"
    .checkHicChroms(x, hicFiles[1]) |>
        expect_null()

})


test_that("develop the function", {

    ## Assign to x (to avoid modifying in place)
    x <- bgi
    seqlevelsStyle(x) <- "ENSEMBL"
    binSize = 50e03
    file = hicFiles[1]

    ## Bin GInteractions if necessary
    x <- .handleBinning(x, binSize)

    ## Ensure seqnames are properly formatted
    .checkHicChroms(x, file)

    ## Convert to short format and sort interactions
    x <- .GInteractionsToShortFormat(x, file)



    ## Don't need to cluster - use grids
    # .clusterPairs(x=x,
    #               radius = binSize,
    #               method = 'manhattan',
    #               pos = 'center')




})
