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

    ## Merge pairs from files
    x <- mergePairs(bedpeFiles)

    ## No columns returns all
    allPairs(x) |>
        length() |>
        expect_equal(20661)

    ## Selected columns
    allPairs(x, source = "FS_5kbLoops.txt") |>
        length() |>
        expect_equal(8566)

    allPairs(x, source = "WT_5kbLoops.txt") |>
        length() |>
        expect_equal(12095)

    ## Wrong column throws error
    allPairs(x, source = "FS_5kbLoops") |>
        expect_error("Source.*not found.")

    ## Merge pairs from objects
    x <- mergePairs(dfs, binSize = 5) |>
        suppressMessages()

    allPairs(x) |>
        length() |>
        expect_equal(4)

    allPairs(x, source = 1) |>
        length() |>
        expect_equal(2)

    allPairs(x, source = 2) |>
        length() |>
        expect_equal(2)

    allPairs(x, source = 3) |>
        expect_error("Source.*not found.")

})
