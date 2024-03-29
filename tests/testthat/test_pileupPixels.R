test_that("pileupPixels works", {
    
    ## Read .hic file paths
    hicFiles <- c(
        marinerData::LEUK_HEK_PJA27_inter_30.hic(),
        marinerData::LEUK_HEK_PJA30_inter_30.hic()
    )
    names(hicFiles) <- c("FS", "WT")
    
    ## Loops
    loops <-
        WT_5kbLoops.txt() |>
        setNames("WT") |>
        read.table(header=TRUE, nrows=1000) |>
        as_ginteractions(keep.extra.columns=FALSE) |>
        assignToBins(binSize=5e3)

    ## Removes the "chr" prefix for compatibility
    ## with the preprocessed hic files
    GenomeInfoDb::seqlevelsStyle(loops) <- 'ENSEMBL'

    ## APA
    mat <- pileupPixels(
        x=loops,
        files=hicFiles[2],
        binSize=5e3,
        minPairDist=50e3,
        normalize=FALSE
    )
    expect_equal(sum(mat), 1133)
    
    mat <- pileupPixels(
        x=loops,
        files=hicFiles[2],
        binSize=5e3,
        minPairDist=50e3,
        normalize=TRUE
    )
    expect_equal(round(sum(mat), 2), 1.24)
    
    mat <- pileupPixels(
        x=loops,
        files=hicFiles[2],
        binSize=5e3,
        minPairDist=50e3,
        normalize=TRUE,
        half="lower"
    )
    expect_true(all(is.na(mat)))

})
