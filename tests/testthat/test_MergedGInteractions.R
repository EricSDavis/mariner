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

## Read in bedpeFiles as a list of GInteractions
giList <-
    lapply(bedpeFiles, fread) |>
    lapply(as_ginteractions)


## Test getters and setters ----------------------------------------------------

test_that("selectionMethod accessor works", {

    x <- mergePairs(x = giList, radius = 10e03)
    selectionMethod(x) |>
        expect_identical("Mean of modes")

    x <- mergePairs(x = giList, radius = 10e03, column = "APScoreAvg")
    selectionMethod(x) |>
        expect_identical("Select by column 'APScoreAvg'")
})

test_that("getPairClusters accessor works", {

    ## Merge pairs and add names
    x <- mergePairs(x = giList, radius = 10e03)
    names(x) <- paste0("loop", 1:length(x))

    x[3:1] |>
        getPairClusters() |>
        names() |>
        expect_identical(paste0("loop", 3:1))

    x[1:3] |>
        getPairClusters() |>
        names() |>
        expect_identical(paste0("loop", 1:3))

    expect_equal(length(getPairClusters(x)), length(x))

})


test_that("subsetBySource method works (and sources accessor)", {

    ## LIMA loops test
    loopFiles <-
        system.file("extdata/lima_loops",
                    package = "mariner") |>
        list.files(full.names = TRUE)

    ## Subsample to 1/4 rows for faster tests
    set.seed(123)
    loops <-
        lapply(loopFiles, fread) |>
        lapply(\(x) x[sample(1:nrow(x), nrow(x)/4, replace = FALSE)]) |>
        lapply(as_ginteractions) |>
        setNames(basename(loopFiles))

    ## Merge loops
    x <- mergePairs(loops, radius = 25e03*5)

    ## Test sources accessor
    expect_identical(sources(x), basename(loopFiles))

    ## subsetBySource method dispatch
    expect_equal(subsetBySource(x) |> length(), 8)
    expect_equal(subsetBySource(x)[[sources(x)[1]]] |> length(), 846)
    expect_equal(subsetBySource(x, include = sources(x)[1]) |> length(), 3982)
    expect_equal(subsetBySource(x, exclude = sources(x)[1]) |> length(), 9890)

    ## Handles cases where none are found
    subsetBySource(x,
           include = sources(x)[1],
           exclude = sources(x)[1]) |>
        length() |>
        expect_equal(0)
    subsetBySource(x,
           include = sources(x)[1:2],
           exclude = sources(x)[1]) |>
        length() |>
        expect_equal(0)
    subsetBySource(x, exclude = sources(x)) |>
        length() |>
        expect_equal(0)

    ## Handles improper source name
    expect_error(subsetBySource(x, include = "foo"),
                 ".*foo.*not source option")
    expect_error(subsetBySource(x, exclude = "foo"),
                 ".*foo.*not source option")
    expect_error(subsetBySource(x, include = ""),
                 ".*not source option")
    expect_error(subsetBySource(x, exclude = ""),
                 ".*not source option")
    expect_error(subsetBySource(x, include = "foo", exclude = "foo"),
                 ".*foo.*not source option")
    expect_error(subsetBySource(x, include = "foo", exclude = "bar"),
                 ".*foo.*bar.*not source option")
    expect_error(subsetBySource(x,
                        include = c(sources(x)[1], "foo"),
                        exclude = "bar"),
                 ".*foo.*bar.*not source option")

    ## subsetBySource(x) produces same result as
    ## using include & exclude
    expect_equal(length(subsetBySource(x)[[sources(x)[1]]]),
                 length(subsetBySource(x,
                               include = sources(x)[1],
                               exclude = sources(x)[2:8])))

    ## subsetBySource(x) produces same result as
    ## include/exclude for all sources
    expect_equal(
        lapply(seq_along(sources(x)), \(i){
            j <- seq_along(sources(x))[-i]
            subsetBySource(x = x,
                           include = sources(x)[i],
                           exclude = sources(x)[j])
        }) |>
            lapply(length) |>
            unlist(),
        subsetBySource(x) |> lapply(length) |> unlist() |> unname()
    )

    ## Results track with getPairClusters(x)
    mgi <- subsetBySource(x,
                          include = sources(x)[1:2],
                          exclude = sources(x)[3:8])
    expect_identical(getPairClusters(mgi) |>
                         lapply(`[[`, "src") |>
                         unlist() |>
                         unique() |>
                         as.character(),
                     sources(x)[1:2])

})


test_that("Aggregating metadata columns works", {

    ## Merge pairs
    x <- mergePairs(x = giList,
                    radius = 5e03)

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
                     expected = getPairClusters(x[1:10]) |>
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
