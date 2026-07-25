library(mariner)
library(rhdf5)

## Fixtures are built by helper-coolFixtures.R. See that file for why these
## are synthesized rather than pulled from marinerData.

coolFile <- makeCoolFile(binSize = 100L)
mcoolFile <- makeMcoolFile(binSizes = c(100L, 200L))

test_that("Cool and mcool files are detected", {
    expect_identical(.checkIfCool(coolFile), ".cool")
    expect_identical(.checkIfCool(mcoolFile), ".mcool")
})

test_that("Non-cooler files are rejected", {
    ## Not HDF5 at all
    txt <- tempfile(fileext = ".txt")
    writeLines("not an hdf5 file", txt)
    expect_error(.checkIfCool(txt), "HDF5")

    ## HDF5, but without the datasets a cooler file must have
    empty <- tempfile(fileext = ".h5")
    h5createFile(empty)
    h5write(1:10, empty, "somethingElse")
    h5closeAll()
    expect_error(.checkIfCool(empty), "expected datasets")
})

test_that("Resolutions are read from cool and mcool files", {
    ## .cool files carry one resolution, inferred from the first bin
    expect_identical(readCoolBpResolutions(coolFile), 100L)

    ## .mcool files list theirs under /resolutions
    expect_identical(readCoolBpResolutions(mcoolFile), c(100L, 200L))
})

test_that("Chromosomes are read from cool files", {
    chroms <- readCoolChroms(coolFile, 100L)
    expect_identical(chroms$name, c("1", "2"))
    expect_identical(chroms$length, c(1000L, 600L))
    expect_identical(chroms$index, 1:2)

    ## mcool resolutions describe the same chromosomes
    expect_identical(
        readCoolChroms(mcoolFile, 100L)[, c("name", "length")],
        readCoolChroms(mcoolFile, 200L)[, c("name", "length")]
    )

    ## invalid resolution is rejected
    expect_error(readCoolChroms(mcoolFile, 12345L), "resolution")
})

test_that("Normalizations are read from cool files", {
    ## `weight` is reported as BALANCE, and NONE is always available
    expect_setequal(readCoolNormTypes(coolFile, 100L), c("BALANCE", "NONE"))
    expect_setequal(readCoolNormTypes(mcoolFile, 200L), c("BALANCE", "NONE"))

    expect_error(readCoolNormTypes(mcoolFile, 12345L), "resolution")
})

test_that("coolStraw returns the expected counts", {
    ## Each row is a query shape. The off-diagonal and interchromosomal cases
    ## matter most: a truncation bug that dropped the last row of the matrix
    ## was invisible for square on-diagonal queries.
    queries <- list(
        list("1", 0, 1000, "1", 0, 1000),    # whole chromosome, end == length
        list("1", 0, 900, "1", 0, 900),      # end lands on a bin start
        list("1", 200, 500, "1", 200, 500),  # interior, on-diagonal
        list("1", 0, 300, "1", 400, 900),    # off-diagonal
        list("1", 0, 1000, "2", 0, 600),     # interchromosomal, full
        list("1", 300, 600, "2", 100, 400),  # interchromosomal, interior
        list("2", 0, 600, "2", 0, 600),      # second chromosome
        list("1", 0, 5000, "1", 0, 5000)     # end past chromosome length
    )

    for (q in queries) {
        label <- sprintf("%s:%s:%s x %s:%s:%s", q[[1]], q[[2]], q[[3]],
                         q[[4]], q[[5]], q[[6]])
        got <- coolStraw(
            norm = "NONE",
            fname = coolFile,
            chr1loc = paste(q[[1]], q[[2]], q[[3]], sep = ":"),
            chr2loc = paste(q[[4]], q[[5]], q[[6]], sep = ":"),
            binsize = 100L
        )
        got <- got[order(got$x, got$y), ]
        rownames(got) <- NULL

        want <- expectedCoolStraw(100L, q[[1]], q[[2]], q[[3]],
                                  q[[4]], q[[5]], q[[6]])
        rownames(want) <- NULL

        expect_equal(got, want, info = label)
    }
})

test_that("coolStraw agrees between cool and mcool files", {
    args <- list(norm = "NONE", chr1loc = "1:0:1000", chr2loc = "2:0:600",
                 binsize = 100L)
    fromCool <- do.call(coolStraw, c(list(fname = coolFile), args))
    fromMcool <- do.call(coolStraw, c(list(fname = mcoolFile), args))
    expect_equal(fromCool, fromMcool)
})

test_that("coolStraw honors a second mcool resolution", {
    got <- coolStraw("NONE", mcoolFile, "1:0:1000", "1:0:1000", 200L)
    got <- got[order(got$x, got$y), ]
    rownames(got) <- NULL

    want <- expectedCoolStraw(200L, "1", 0, 1000, "1", 0, 1000)
    rownames(want) <- NULL

    expect_equal(got, want)
})

test_that("coolStraw applies balancing weights", {
    binSize <- 100L
    got <- coolStraw("BALANCE", coolFile, "1:200:500", "1:200:500", binSize)
    got <- got[order(got$x, got$y), ]

    raw <- expectedCoolStraw(binSize, "1", 200, 500, "1", 200, 500)
    weights <- .coolWeights(nrow(.coolBins(binSize)))
    ## cooler defines the balanced value as A_ij * w_i * w_j, so weights are
    ## multiplied rather than divided out.
    ## bins 2..5 cover 200-500 on chromosome 1
    ids <- cbind(raw$x %/% binSize, raw$y %/% binSize)
    want <- raw$counts * (weights[ids[, 1] + 1L] * weights[ids[, 2] + 1L])

    expect_equal(got$counts, want)

    ## and that BALANCE actually differs from NONE
    expect_false(isTRUE(all.equal(got$counts, raw$counts)))
})

test_that("coolStraw validates its arguments", {
    expect_error(
        coolStraw("NONE", coolFile, "1:0:1000", "1:0:1000", 12345L),
        "binsize"
    )
    expect_error(
        coolStraw("NOPE", coolFile, "1:0:1000", "1:0:1000", 100L),
        "norm"
    )
    expect_error(
        coolStraw("NONE", coolFile, "not-a-location", "1:0:1000", 100L),
        "chr1loc"
    )
    expect_error(
        coolStraw("NONE", coolFile, "1:0:1000", "99:0:1000", 100L),
        "not found"
    )
})

test_that("Cool args are checked correctly", {
    ## check that it works
    .checkCoolArgs(
        files = coolFile,
        half = "both",
        norm = "BALANCE",
        binSize = 100L,
        matrix = "observed"
    ) |>
        expect_null()

    ## error if matrix is not "observed"
    .checkCoolArgs(
        files = coolFile,
        half = "both",
        norm = "BALANCE",
        binSize = 100L,
        matrix = "oe"
    ) |>
        expect_error("matrix")

    ## error if `half` is invalid
    .checkCoolArgs(
        files = coolFile,
        half = "neither",
        norm = "BALANCE",
        binSize = 100L,
        matrix = "observed"
    ) |>
        expect_error("half")

    ## error if `norm` is not in files
    .checkCoolArgs(
        files = coolFile,
        half = "both",
        norm = "KR",
        binSize = 100L,
        matrix = "observed"
    ) |>
        expect_error("norm")

    ## error if `binSize` is invalid
    .checkCoolArgs(
        files = coolFile,
        half = "both",
        norm = "BALANCE",
        binSize = 101L,
        matrix = "observed"
    ) |>
        expect_error("binSize")
})
