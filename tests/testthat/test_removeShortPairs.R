library(mariner)
library(marinerData)

test_that("removeShortPairs works as expected", {

    ## Better example
    ## Overlapping or equal are always no
    ## Padding should be less than or equal to
    gi <- as_ginteractions(read.table(
        text="
        seqnames1 start1 end1 seqnames2 start2 end2 keep
        chr1 300 400 chr1 300 400 'no'
        chr1 100 200 chr1 300 400 'yes_padding<=100'
        chr1 300 400 chr1 100 200 'yes_padding<=100'
        chr1 300 400 chr2 300 400 'yes'
        chr1 250 350 chr1 300 400 'no'
        chr1 300 400 chr1 250 350 'no'
        chr1 300 400 chr1 700 800 'yes_padding<=300'
        ",
        header=TRUE
    ))


    ## Everything with no padding
    ans <- removeShortPairs(gi)
    expect_identical(
        object=ans$keep,
        expected=gi$keep[startsWith(gi$keep, "yes")]
    )

    ## Accepts padding equal to break point
    ans <- removeShortPairs(gi, padding=100)
    expect_identical(
        object=ans$keep,
        expected=gi$keep[startsWith(gi$keep, "yes")]
    )
    ans <- removeShortPairs(gi, padding=300)
    expect_identical(
        object=ans$keep,
        expected=grep("yes$|.*300", gi$keep, value=TRUE)
    )

    ## Accepts padding less than break point
    ans <- removeShortPairs(gi, padding=99)
    expect_identical(
        object=ans$keep,
        expected=gi$keep[startsWith(gi$keep, "yes")]
    )
    ans <- removeShortPairs(gi, padding=299)
    expect_identical(
        object=ans$keep,
        expected=grep("yes$|.*300", gi$keep, value=TRUE)
    )

    ## Removes pairs greater than all breakpoints
    ans <- removeShortPairs(gi, padding=301)
    expect_identical(
        object=ans$keep,
        expected="yes"
    )

})

test_that("complex test for removeShortPairs", {

    ## Hi-C files
    hicFiles <- c(
        LEUK_HEK_PJA27_inter_30.hic(),
        LEUK_HEK_PJA30_inter_30.hic()
    ) |> setNames(c("FS", "WT")) |>
        suppressMessages()

    ## Loops from each condition
    loopList <- c(FS_5kbLoops.txt(), WT_5kbLoops.txt()) |>
        setNames(c("FS", "WT")) |>
        lapply(read.table, header=TRUE, nrows=1000) |>
        lapply(as_ginteractions)

    ## Merge
    loops <-
        mergePairs(loopList, radius=100e3) |>
        GenomeInfoDb::`seqlevelsStyle<-`('ENSEMBL') |>
        assignToBins(binSize=100e3)

    ## Expand to apa regions
    regions <- pixelsToMatrices(loops, buffer=5)

    ## All regions
    everything <-
        regions |>
        pullHicMatrices(
            files=hicFiles,
            binSize=100e3,
            half='upper'
        )

    ## Filter out short pairs
    filtered <-
        regions |>
        removeShortPairs() |>
        pullHicMatrices(
            files=hicFiles,
            binSize=100e3,
            half='upper'
        )

    ## Expect NA values
    expect_true(anyNA(aggHicMatrices(everything, verbose=FALSE)))

    ## Expect no NAs
    expect_false(anyNA(aggHicMatrices(filtered, verbose=FALSE)))

    ## Loops are filtered out correctly
    (
        loops |> length() -
            loops |> removeShortPairs(padding=40e3) |> length() ==
            which(pairdist(loops, type="gap") <= 40e3) |> length()
    ) |>
        expect_true()

})
