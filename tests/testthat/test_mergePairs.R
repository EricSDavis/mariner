library(mariner)
library(data.table)

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

giList <-
    lapply(bedpeFiles, fread) |>
    lapply(as_ginteractions)

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

test_that("List type can not return data.frame-like objects", {

    .checkListFormat(list(data.frame())) |>
        expect_error("^.*Your list is type data.frame.")

    .checkListFormat(list(data.table::data.table())) |>
        expect_error("^.*Your list is type data.table, data.frame.")

    .checkListFormat(list(S4Vectors::DataFrame())) |>
        expect_error("^.*Your list is type DFrame.")
})

test_that("List type can return GInteractions", {

    .checkListFormat(list(GInteractions())) |>
        expect_null()

    .checkListFormat(list(GInteractions(), GInteractions())) |>
        expect_null()
})

test_that("Type mixture is not accepted", {

    .checkListFormat(list(data.frame(), GInteractions())) |>
        expect_error("^.*Your list is type data.frame, GInteractions.")

    .checkListFormat(list(NULL, NA)) |>
        expect_error("^.*Your list is type NULL, logical.")
})


## Test mergePairs dispatch methods --------------------------------------------

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

    ## Manhattan distance of 0 and binSize of 10 is radius = 0
    .findClusters(x = dt[,c('start1','start2')],
                  radius = 0,
                  method = "manhattan") |>
        expect_equal(c(0,0,0,0,0))

    ## Manhattan distance of 1 and binSize of 10 is radius = 10
    .findClusters(x = dt[,c('start1','start2')],
                  radius = 10,
                  method = "manhattan") |>
        expect_equal(c(1,1,0,0,0))

    ## Manhattan distance of 2 and binSize of 10 is radius = 20
    .findClusters(x = dt[,c('start1','start2')],
                  radius = 20,
                  method = "manhattan") |>
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

    ## Manhattan distance of 0 and binSize of 10 is radius = 0
    dt[,clst := .findClusters(x = .SD[,c('start1','start2')],
                              radius = 10,
                              method = "manhattan"),
       by = .(seqnames1, seqnames2)]
    dt$clst |> expect_equal(c(1,1,0,0,0))

    ## Manhattan distance of 2 and binSize of 10 is radius = 20
    dt <- GInteractions(gr1, gr2) |> as.data.table()
    dt[,clst := .findClusters(x = .SD[,c('start1','start2')],
                              radius = 20,
                              method = "manhattan"),
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

    ## Manhattan distance of 1 and binSize of 10 is radius = 10
    dt[,clst := .findClusters(x = .SD[,c('start1','start2')],
                              radius = 10,
                              method = 'manhattan'),
       by = .(seqnames1, seqnames2)]
    dt[,grp := .GRP, by = .(seqnames1, seqnames2)]
    dt$clst |> expect_equal(c(1,1,0,0,0,1,1))
    dt$grp |> expect_equal(c(1,1,1,1,1,2,2))
})

test_that("Removing changed column names works", {
    .renameCols(c("id_1", "src_1", "grp_1", "clst_1", "mid1_1", "mid2_1")) |>
        expect_identical(c("id", "src", "grp", "clst", "mid1", "mid2"))
})

test_that("id, src, grp, clst, mid1, and mid2 column names can be used", {

    ## Add an id column to merged pairs
    mp <-
        bedpeFiles |>
        lapply(fread) |>
        lapply(\(x) {
            x$id <- rev(seq_len(nrow(x)))
            x$src <- "src"
            x$grp <- "grp"
            x$clst <- "clst"
            x$mid1 <- rowMeans(x[,c('x1', 'x2')])
            x$mid2 <- rowMeans(x[,c('y1', 'y2')])
            x
        }) |>
        lapply(as_ginteractions) |>
        mergePairs(radius = 0, column = "APScoreAvg")

    grep("id|src|grp|clst|mid1|mid2", colnames(mcols(mp)), value = TRUE) |>
        expect_identical(c("id", "src", "grp", "clst", "mid1", "mid2"))

    expect_length(mp$id, length(mp))

})

test_that("Bad column name throws error.", {
    mergePairs(x = giList, radius = 0, column = "foo") |>
        expect_error("^Column.*does not exist.")
})

test_that("Remove metadata when using mean of modes but not column.", {
    x <- mergePairs(x = giList, radius = 0)
    expect_equal(ncol(mcols(x)), 0)

    x <- mergePairs(x = giList, radius = 0, column = "APScoreAvg")
    expect_equal(ncol(mcols(x)), 9)
})

test_that("selectMax parameter works", {

    gi <-
        GInteractions(anchor1 = GRanges(seqnames = "chr1",
                                        ranges = IRanges(start = c(1, 1),
                                                         end = c(10, 10))),
                      anchor2 = GRanges(seqnames = "chr1",
                                        ranges = IRanges(start = c(1, 1),
                                                         end = c(10, 10))),
                      value = c(10, 5))

    mergePairs(x=gi, radius = 0,
               column = "value", selectMax = TRUE)$value |>
        expect_equal(10)

    mergePairs(x=gi, radius = 0,
               column = "value", selectMax = FALSE)$value |>
        expect_equal(5)

})

test_that("pos parameter works", {

    gi <-
        GInteractions(anchor1 = GRanges(seqnames = "chr1",
                                        ranges = IRanges(start = c(1, 1, 1),
                                                         end = c(10,
                                                                 100,
                                                                 1e03))),
                      anchor2 = GRanges(seqnames = "chr1",
                                        ranges = IRanges(start = c(1, 1, 1),
                                                         end = c(10, 10, 10))))
    mergePairs(x = gi,
               radius = 50,
               method = "manhattan",
               pos = "start") |>
        length() |>
        expect_equal(1)

    mergePairs(x = gi,
               radius = 50,
               method = "manhattan",
               pos = "center") |>
        length() |>
        expect_equal(2)

    mergePairs(x = gi,
               radius = 50,
               method = "manhattan",
               pos = "end") |>
        length() |>
        expect_equal(3)
})
