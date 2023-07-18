library(mariner)
library(marinerData)

## Shared objects --------------------------------------------------------------
## Read .hic file paths
hicFiles <- c(
    LEUK_HEK_PJA27_inter_30.hic(),
    LEUK_HEK_PJA30_inter_30.hic()
)
names(hicFiles) <- c("FS", "WT")

## Read in loops as GInteractions object
loops <-
    WT_5kbLoops.txt() |>
    setNames("WT") |>
    read.table(header=TRUE) |>
    as_ginteractions(keep.extra.columns=FALSE)

## Removes the "chr" prefix for compatibility
## with the preprocessed hic files
GenomeInfoDb::seqlevelsStyle(loops) <- 'ENSEMBL'

## Calculate loop enrichment
enrich <- calcLoopEnrichment(
    x=assignToBins(loops, 100e03),
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
