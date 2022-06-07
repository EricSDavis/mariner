library(mariner)
## Shared objects --------------------------------------------------------------


## Test .checkListType() -------------------------------------------------------

test_that("Empty lists are not accepted", {
    .checkListType(list()) |>
        expect_error("^.*length > 0.")
})

test_that("Other types types are not accepted", {

    .checkListType(list(NA)) |>
        expect_error("^.*Your list is type logical.")

    .checkListType(list(NULL)) |>
        expect_error("^.*Your list is type NULL.")

    .checkListType(list(NA_character_)) |>
        expect_error("^.*Your list is type character.")

    .checkListType(list(NA_integer_)) |>
        expect_error("^.*Your list is type integer.")

    .checkListType(list(1,2,3)) |>
        expect_error("^.*Your list is type numeric.")
})

test_that("List type can return data.frame-like objects", {

    .checkListType(list(data.frame())) |>
        expect_identical('data.frame')

    .checkListType(list(data.table::data.table())) |>
        expect_identical('data.frame')

    .checkListType(list(S4Vectors::DataFrame())) |>
        expect_identical('data.frame')
})

test_that("List type can return GInteractions", {

    .checkListType(list(GInteractions())) |>
        expect_identical("GInteractions")

    .checkListType(list(GInteractions(), GInteractions())) |>
        expect_identical("GInteractions")
})

test_that("Type mixture is not accepted", {

    .checkListType(list(data.frame(), GInteractions())) |>
        expect_error("^.*must.*be.*same type.*")

    .checkListType(list(NULL, NA)) |>
        expect_error("^.*Your list is type NULL, logical.")
})
