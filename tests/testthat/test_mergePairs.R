library(mariner)
library(marinerData)
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
bedpeFiles <- c(
    FS_5kbLoops.txt(),
    WT_5kbLoops.txt()
)
names(bedpeFiles) <- c("FS", "WT")

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

test_that("S4Vectors List() and SimpleList() types are accepted", {

    ## Build GInteractions object
    gr1 <-
        GRanges(seqnames = "chr1",
                ranges = IRanges(start = c(30,40,40,70,80),
                                 end = c(40,50,50,80,90)))
    gr2 <-
        GRanges(seqnames = "chr1",
                ranges = IRanges(start = c(30,30,50,10,30),
                                 end = c(40,40,60,20,40)))
    gi <- GInteractions(gr1, gr2)

    ## Test List() and SimpleList()
    gil <- .readGInteractionsList(list(gi, gi), column=NULL)
    .readGInteractionsList(list(gi, gi), column=NULL) |>
        expect_identical(gil)
    .readGInteractionsList(S4Vectors::SimpleList(gi, gi), column=NULL) |>
        expect_identical(gil)

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

test_that("Mean of modes option doesn't alter original ranges", {

    ## Define start and end of grid,
    ## binSize, and n interactions
    binSize <- 5
    n <- 10
    breaks <- seq(0, 60, binSize)

    ## Generate grid of options to avoid
    ## duplicating pixels
    opts <- expand.grid(
        bin1 = breaks[-length(breaks)],
        bin2 = breaks[-length(breaks)]
    )

    ## Randomly select pixels
    set.seed(123)
    s <- sample(seq_len(nrow(opts)), n, replace = FALSE)
    cnts <- sample(breaks, n, replace = TRUE)
    r1 <- opts[s, 1]
    r2 <- opts[s, 2]

    ## Construct GInteractions object
    library(InteractionSet)
    gi <- GInteractions(
        anchor1 = GRanges(
            seqnames = "chr1",
            ranges = IRanges(
                start = r1,
                width = binSize + 1
            )
        ),
        anchor2 = GRanges(
            seqnames = "chr1",
            ranges = IRanges(
                start = r2,
                width = binSize + 1
            )
        ),
        count = cnts
    )

    mgi1 <- mergePairs(x=gi, radius=binSize*2, column='count')
    mgi2 <- mergePairs(x=gi, radius=binSize*2)
    expect_identical(
        clusters(mgi1),
        clusters(mgi2)
    )

})

test_that("radius and method columns don't cause issues", {

    ## Create example giList
    giList <- list(
        read.table(text="chr1 10 20 chr1 50 60") |>
            as_ginteractions(),
        read.table(text="chr1 5 15 chr1 45 55") |>
            as_ginteractions()
    )

    ## Add "radius" column to a giList
    giList[[1]]$radius <- 1
    giList[[2]]$radius <- 2
    giList[[1]]$method <- 2
    giList[[2]]$method <- 1
    mgi <- mergePairs(x=giList, radius=10, column="radius")
    expect_identical(mgi$radius, 2)
    expect_identical(mgi$method, 1)

    mgi <- mergePairs(x=giList, radius=10, column="method")
    expect_identical(mgi$radius, 1)
    expect_identical(mgi$method, 2)

})
