library(mariner)
library(marinerData)

test_that("removeShortPairs works as expected", {

    ## Example GInteractions object
    gi <- as_ginteractions(read.table(
        text="
        seqnames1 start1 end1 seqnames2 start2 end2 keep
        chr1 300 400 chr1 300 400 'no'
        chr1 100 200 chr1 300 400 'yes'
        chr1 300 400 chr1 100 200 'yes'
        chr1 300 400 chr2 300 400 'yes'
        chr1 250 350 chr1 300 400 'only_with_padding_50'
        chr1 300 400 chr1 250 350 'only_with_padding_50'
        ",
        header=TRUE
    ))

    ## Test internal function
    ans <- .removeShortPairs(gi, padding=0)
    expect_identical(
        object = ans$keep,
        expected = rep("yes", 3)
    )

    ## Default of padding=0
    ans <- removeShortPairs(gi)
    expect_identical(
        object = ans$keep,
        expected = rep("yes", 3)
    )

    ## Not enough padding
    ans <- removeShortPairs(gi, padding=10)
    expect_identical(
        object = ans$keep,
        expected = rep("yes", 3)
    )

    ## Enough padding
    ans <-removeShortPairs(gi, padding=50)
    expect_identical(
        object = ans$keep,
        expected = c(
            rep("yes", 3),
            rep("only_with_padding_50", 2)
        )
    )

    ## Plenty of padding
    ans <- removeShortPairs(gi, padding=100)
    expect_identical(
        object = ans$keep,
        expected = c(
            rep("yes", 3),
            rep("only_with_padding_50", 2)
        )
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
        binPairs(binSize=100e3)

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

})
