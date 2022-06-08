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

## Helper function to extract the source column
extractSource <- \(x) {
    lapply(x, `[[`, 'source') |>
        lapply(unique) |>
        unlist() |>
        unname()
}


## Test .checkListFormat() -------------------------------------------------------

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


## Test .readBedpeFromList() ------------------------------------------

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

## Test .mergePairs dispatch methods -------------------------------------------

test_that("Source column naming works for .mergePairsList", {

    ## Unnamed list uses indices
    .mergePairsList(x = dfs, binSize = 3) |>
        extractSource() |>
        expect_equal(seq_along(dfs))

    ## Named list uses names
    dfs |>
        `names<-`(value = c("set1", "set2")) |>
        .mergePairsList(binSize = 3) |>
        extractSource() |>
        expect_equal(c('set1', 'set2'))
})

test_that("Source column naming works for .mergePairsCharacter", {

    ## Unnamed list uses indices
    .mergePairsCharacter(x = bedpeFiles, binSize = 3) |>
        extractSource() |>
        expect_equal(basename(bedpeFiles))

})
