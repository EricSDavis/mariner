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


## Test .checkListFormat() -----------------------------------------------------

test_that("Empty lists are not accepted", {
    .checkListFormat(list()) |>
        expect_error("^.*length > 0.")
})

test_that("Other types types are not accepted", {

    .checkListFormat(list(NA)) |>
        expect_error("^.*Your list is type logical.")

    .checkListFormat(list(NULL)) |>
        expect_error("^.*Your list is type NULL.")

    .checkListFormat(list(NA_character_)) |>
        expect_error("^.*Your list is type character.")

    .checkListFormat(list(NA_integer_)) |>
        expect_error("^.*Your list is type integer.")

    .checkListFormat(list(1,2,3)) |>
        expect_error("^.*Your list is type numeric.")
})

test_that("List type can return data.frame-like objects", {

    .checkListFormat(list(data.frame())) |>
        expect_null()

    .checkListFormat(list(data.table::data.table())) |>
        expect_null()

    .checkListFormat(list(S4Vectors::DataFrame())) |>
        expect_null()
})

test_that("List type can return GInteractions", {

    .checkListFormat(list(GInteractions())) |>
        expect_null()

    .checkListFormat(list(GInteractions(), GInteractions())) |>
        expect_null()
})

test_that("Type mixture is not accepted", {

    .checkListFormat(list(data.frame(), GInteractions())) |>
        expect_error("^.*must.*be.*same type.*")

    .checkListFormat(list(NULL, NA)) |>
        expect_error("^.*Your list is type NULL, logical.")
})


## Test .readBedpeFromList() ---------------------------------------------------

test_that("data.frame lists can be read", {

    .readBedpeFromList(dfs) |>
        expect_snapshot_output()

})

test_that("data.table lists can be read", {

    library(data.table)
    lapply(dfs, as.data.table) |>
        .readBedpeFromList() |>
        expect_snapshot_output()
})

test_that("DataFrame lists can be read", {

    library(S4Vectors)
    lapply(dfs, DataFrame) |>
        .readBedpeFromList() |>
        expect_snapshot_output()
})

test_that("GInteractions list can be read", {

    dfs |>
        lapply(as_ginteractions) |>
        .readBedpeFromList() |>
        expect_snapshot_output()

})

## Test mergePairs dispatch methods --------------------------------------------

## Need to write source accessor
# test_that("Source column naming works for .mergePairsCharacter", {
#
#     ## Unnamed list uses indices
#     .mergePairsCharacter(x = bedpeFiles, binSize = 3, radius = 2) |>
#         {\(x) unique(x$source)}() |>
#         suppressWarnings() |>
#         expect_equal(basename(bedpeFiles))
#
#     .mergePairsCharacter(bedpeFiles, 5000, 2) |>
#         {\(x) unique(x$source)}() |>
#         expect_equal(basename(bedpeFiles))
#
# })

test_that("mergePairs warns about binning", {

    .mergePairsCharacter(x = bedpeFiles, binSize = 5e03, radius = 2) |>
        expect_warning(NA)

    mergePairs(x = bedpeFiles, binSize = 10e03, radius = 2) |>
        expect_message("^.*not binned to `binSize`.*")

    mergePairs(x = dfs, binSize = 2, radius = 2) |>
        expect_message("^.*not binned to `binSize`.*")

})

test_that("Find overlaps by manhattan distance works", {
    library(InteractionSet)
    library(GenomicRanges)
    library(data.table)

    ## Define anchor regions to test distance clustering
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

    ## Manhattan distance (radius) of 0
    .findClusters(x = dt[,c('start1','start2')],
                  radius = 0,
                  binSize = 10) |>
        expect_equal(c(0,0,0,0,0))

    ## Manhattan distance (radius) of 1
    .findClusters(x = dt[,c('start1','start2')],
                  radius = 1,
                  binSize = 10) |>
        expect_equal(c(1,1,0,0,0))

    ## Manhattan distance (radius) of 2
    .findClusters(x = dt[,c('start1','start2')],
                  radius = 2,
                  binSize = 10) |>
        expect_equal(c(1,1,1,0,0))
})

test_that("Find overlaps by group with data.table", {
    library(InteractionSet)
    library(GenomicRanges)
    library(data.table)

    ## Define anchor regions to test distance clustering
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

    ## Manhattan distance (radius) of 1
    dt[,clst := .findClusters(x = .SD[,c('start1','start2')],
                              radius = 1,
                              binSize = 10),
       by = .(seqnames1, seqnames2)]
    dt$clst |> expect_equal(c(1,1,0,0,0))

    ## Manhattan distance (radius) of 2
    dt <- GInteractions(gr1, gr2) |> as.data.table()
    dt[,clst := .findClusters(x = .SD[,c('start1','start2')],
                              radius = 2,
                              binSize = 10),
       by = .(seqnames1, seqnames2)]
    dt$clst |> expect_equal(c(1,1,1,0,0))
})

test_that("Handle interchromosomal by group", {
    library(InteractionSet)
    library(GenomicRanges)
    library(data.table)

    ## Define anchor regions to test distance clustering
    gr1 <-
        GRanges(seqnames = "chr1",
                ranges = IRanges(start = c(30,40,40,70,80,30,30),
                                 end = c(40,50,50,80,90,40,40)))
    gr2 <-
        GRanges(seqnames = c(rep("chr1",5),"chr2", "chr2"),
                ranges = IRanges(start = c(30,30,50,10,30,30,30),
                                 end = c(40,40,60,20,40,40,40)))

    ## Form GInteractions and convert to data.table
    dt <- GInteractions(gr1, gr2) |> as.data.table()

    ## Manhattan distance (radius) of 1
    dt[,clst := .findClusters(x = .SD[,c('start1','start2')],
                              radius = 1,
                              binSize = 10),
       by = .(seqnames1, seqnames2)]
    dt[,grp := .GRP, by = .(seqnames1, seqnames2)]
    dt$clst |> expect_equal(c(1,1,0,0,0,1,1))
    dt$grp |> expect_equal(c(1,1,1,1,1,2,2))
})

test_that("Removing changed column names works", {
    .renameCols(c("id_1", "src_1", "grp_1", "clst_1")) |>
        expect_identical(c("id", "src", "grp", "clst"))
})

test_that("id, src, grp, and clst column names can be used", {

    ## Add an id column to merged pairs
    mp <-
        bedpeFiles |>
        lapply(fread) |>
        lapply(\(x) {x$id <- rev(seq_len(nrow(x))); x}) |>
        mergePairs(binSize = 5000, radius = 0)

    expect_length(mp$id, length(mp))

})
