library(mariner)

test_that("plotBullseye works on a 21x21 matrix", {

    plotBullseye(x = matrix(1:(21*21), 21, 21)) |>
        expect_error(NA)

})
