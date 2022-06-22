library(mariner)
## Shared objects --------------------------------------------------------------

## Create some in-memory data.frames
dfs <-
    list(
        data.frame(seqnames1="chr1", start1=c(1,10), end1=c(9, 19),
                   seqnames2="chr1", start2=c(1,10), end2=c(9, 19)),
        data.frame(seqnames1="chr1", start1=c(5,10), end1=c(15, 19),
                   seqnames2="chr1", start2=c(5,10), end2=c(15, 19))
    )

## Reference BEDPE files (loops called with SIP)
bedpeFiles <-
    system.file("extdata", package = "mariner") |>
    list.files(pattern = "Loops.txt", full.names = TRUE)


## Test getters and setters ----------------------------------------------------

test_that("selectionMethod accessor works", {
    x <- mergePairs(bedpeFiles)
    selectionMethod(x) |>
        expect_identical("Mean of modes")

    x <- mergePairs(bedpeFiles, column = "APScoreAvg")
    selectionMethod(x) |>
        expect_identical("Select by column 'APScoreAvg'")
})

test_that("allPairs accessor works", {

    ## Merge pairs and add names
    x <- mergePairs(bedpeFiles)
    names(x) <- paste0("loop", 1:length(x))

    x[3:1] |>
        allPairs() |>
        names() |>
        expect_identical(paste0("loop", 3:1))

    x[1:3] |>
        allPairs() |>
        names() |>
        expect_identical(paste0("loop", 1:3))

    expect_equal(length(allPairs(x)), length(x))


    ## Use data.frame input
    x <- mergePairs(dfs, binSize = 5)
    expect_equal(length(allPairs(x)), length(x))
})
