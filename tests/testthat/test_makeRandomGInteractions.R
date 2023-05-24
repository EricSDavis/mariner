library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomeInfoDb)

test_that("Making random GRanges", {
    testthat::skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
    testthat::skip_if_not_installed("GenomeInfoDb")

    ## Define your seqinfo object w/chromosome info
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    si <- seqinfo(txdb)
    si <- keepStandardChromosomes(si)

    set.seed(123)
    a1 <- .makeRandomGRanges(si, 100)
    expect_identical(length(a1), 2L)
    a2 <- .makeRandomGRanges(si, 100, .rows=a1$rows)
    expect_identical(a1$rows, a2$rows)

    set.seed(123)
    gr <- makeRandomGRanges(seqinfo=si)
    expect_s4_class(gr, "GRanges")
    expect_identical(length(gr), 100L)


})

test_that("Making random GInteractions", {
    testthat::skip_if_not_installed("TxDb.Hsapiens.UCSC.hg38.knownGene")
    testthat::skip_if_not_installed("GenomeInfoDb")

    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    si <- seqinfo(txdb)
    si <- keepStandardChromosomes(si)

    set.seed(123)
    gi <- makeRandomGInteractions(si, n=100, interchromosomal=FALSE)
    expect_true(all(InteractionSet::intrachr(gi)))
    expect_identical(length(gi), 100L)
    expect_identical(seqinfo(gi), si)
    expect_s4_class(gi, "GInteractions")

    set.seed(123)
    gi <- makeRandomGInteractions(si, n=100)
    expect_false(all(InteractionSet::intrachr(gi)))

})
