library(mariner)
library(GenomicRanges)
library(InteractionSet)
library(data.table)
library(S4Vectors)

## Shared objects --------------------------------------------------------------

## Define example anchor regions
gr1 <-
    GRanges(seqnames = "chr1",
            ranges = IRanges(start = c(30,40,40,70,80),
                             end = c(40,50,50,80,90)))
gr2 <-
    GRanges(seqnames = "chr1",
            ranges = IRanges(start = c(30,30,50,10,30),
                             end = c(40,40,60,20,40)))

## Form GInteractions and convert to data.table
dt <- GInteractions(gr1, gr2) |> as.data.table()

## Split into two files
dts <- split(dt, c(rep(1,3), rep(2, 2)))

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


    ## Use data.table input
    x <- mergePairs(dts, binSize = 5)
    expect_equal(length(allPairs(x)), length(x))
})


test_that("deNovo method works (and sources accessor)", {

    ## Merge pairs
    x <- mergePairs(dts, binSize = 10, radius = 2)

    ## deNovo pairs
    lapply(deNovo(x), length) |>
        unname() |>
        unlist() |>
        expect_equal(c(2,1))


    ## LIMA loops test (subsample these for faster tests)
    loopFiles <- list.files("inst/extdata/lima_loops",
                            full.names = TRUE)
    x <- mergePairs(loopFiles, binSize = 10e03, radius = 5, column = "APScoreAvg")

    ## Test sources accessor
    sources(x)

    ## deNovo method dispatch
    deNovo(x)
    deNovo(x)[sources(x)[1]]
    deNovo(x, include = sources(x)[1])
    deNovo(x, exclude = sources(x)[1])

    ## deNovo(x) produces same result as using include & exclude
    expect_equal(length(deNovo(x)[[sources(x)[1]]]),
                 length(deNovo(x,
                               include = sources(x)[1],
                               exclude = sources(x)[2:8])))

    ## deNovo(x) == include/exclude for all sources
    expect_equal(
        lapply(seq_along(sources(x)), \(i){
            j <- seq_along(sources(x))[-i]
            deNovo(x = x, include = sources(x)[i], exclude = sources(x)[j])
        }) |>
            lapply(length) |>
            unlist(),
        deNovo(x) |> lapply(length) |> unlist() |> unname()
    )


    ## Include/exclude work as expected
    mgi <- deNovo(x, include = sources(x)[1:2], exclude = sources(x)[3:8])
    expect_identical(allPairs(mgi) |>
                         lapply(`[[`, "src") |>
                         unlist() |>
                         unique(),
                     sources(x)[1:2])


    ## Smaller test set
    # gi <-
    #     c(GInteractions(GRanges("chr1:20-30"), GRanges("chr1:40-50"), name=1),
    #       GInteractions(GRanges("chr1:20-30"), GRanges("chr1:30-40"), name=2),
    #       GInteractions(GRanges("chr1:50-60"), GRanges("chr1:60-70"), name=4),
    #       GInteractions(GRanges("chr1:60-70"), GRanges("chr1:60-70"), name=5),
    #       GInteractions(GRanges("chr1:80-90"), GRanges("chr1:30-40"), name=6))
    # gi2 <-
    #     GInteractions(GRanges("chr1:30-40"), GRanges("chr1:40-50"), name=3)
    # x <- mergePairs(list(gi, gi2), binSize = 10)
    # allPairs(x)


})


test_that("Aggregating metadata columns works", {

    ## Merge pairs
    x <- mergePairs(x = bedpeFiles,
                    binSize = 5e03,
                    radius = 0)

    aggPairMcols(x, columns = c("APScoreAvg", "avg"), fun = "mean")

    ## Testing character vs function input
    aggPairMcols(x, columns = "APScoreAvg", funs = "mean") |>
        mcols() |>
        colnames() |>
        expect_identical("mean.APScoreAvg")

    aggPairMcols(x, columns = "APScoreAvg", funs = mean) |>
        mcols() |>
        colnames() |>
        expect_identical("fun1.APScoreAvg")

    ## Testing values
    expect_identical(object = aggPairMcols(x[1:10],
                                           columns = "APScoreAvg",
                                           funs = "mean")$mean.APScoreAvg,
                     expected = allPairs(x[1:10]) |>
                         lapply(`[[`, "APScoreAvg") |>
                         lapply(mean) |>
                         unlist())

    ## Multiple input types
    aggPairMcols(x,
                 columns = "APScoreAvg",
                 fun = c("mean", min, \(x) mean(x))) |>
        mcols() |>
        colnames() |>
        expect_identical(c("mean.APScoreAvg",
                           "fun2.APScoreAvg",
                           "fun3.APScoreAvg"))

    aggPairMcols(x,
                 columns = c("APScoreAvg", "avg"),
                 fun = "mean") |>
        mcols() |>
        colnames() |>
        expect_identical(c("mean.APScoreAvg",
                           "mean.avg"))

    ## Throw errors
    aggPairMcols(x, columns = "APScoreAvg", funs = \(x) x) |>
        expect_error("Improper aggregation function.")

    aggPairMcols(x, columns = c("blah", "foo"), funs = "mean") |>
        expect_error("^Column.*not exist.$")

})
