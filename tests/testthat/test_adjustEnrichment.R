library(mariner)

## Shared objects --------------------------------------------------------------
## Read .hic file paths
hicFiles <-
    system.file("extdata/test_hic", package="mariner") |>
    list.files(pattern=".hic", full.names=TRUE)

## Read in loops as GInteractions object
loops <-
    system.file("extdata", package="mariner") |>
    list.files(pattern="WT.*Loops.txt", full.names=TRUE) |>
    read.table(header=TRUE) |>
    as_ginteractions(keep.extra.columns=FALSE)

## Removes the "chr" prefix for compatibility
## with the preprocessed hic files
GenomeInfoDb::seqlevelsStyle(loops) <- 'ENSEMBL'

## Calculate loop enrichment
enrich <- calcLoopEnrichment(
    x=binPairs(loops, 100e03),
    files=hicFiles
)

test_that("adjusting enrichment for distance", {

    ans <- .rollEnrich(x=loops, scores=enrich[,1], k=50)
    expect_equal(nrow(ans), 754)

    .adjustEnrich(x=loops, scores=enrich[,1]) |>
    expect_error(NA)

    ans <- .adjustEnrichment(enrich, interactions=loops, k=25, nknots=10)
    expect_equal(dim(ans), c(12095, 2))

    ans <- adjustEnrichment(enrich, loops)
    expect_equal(dim(ans), c(12095, 2))

    .plotEnrich(enrich[,1], loops, k=25, nknots=10, plot=FALSE) |>
        expect_error(NA)

    plotEnrichment(enrich[,2], loops, plot=FALSE) |>
        expect_error(NA)

})
