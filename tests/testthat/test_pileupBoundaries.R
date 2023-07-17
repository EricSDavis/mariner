test_that("pileupBoundaries works", {
    
    ## Read .hic file paths
    hicFile <- marinerData::LEUK_HEK_PJA27_inter_30.hic()
    names(hicFile) <- "FS"
    
    ## Loops
    loops <- 
        marinerData::FS_5kbLoops.txt() |>
        read.table(header=TRUE, nrows=100) |>
        as_ginteractions() |>
        GenomeInfoDb::`seqlevelsStyle<-`(value='ENSEMBL')
    
    ## Warn about small binSize
    pileupBoundaries(
        x=loops,
        files=hicFile,
        binSize=10e3,
    ) |>
        expect_message("Caution") |>
        suppressMessages() |>
        suppressWarnings()
    
    ## Expect DelayedMatrix out
    boundaries <- pileupBoundaries(
        x=loops,
        files=hicFile,
        binSize=10e3,
    ) |>
        suppressMessages() |>
        suppressWarnings()
    
    expect_s4_class(boundaries, "DelayedMatrix")
    
    ## Accepts GRanges
    pileupBoundaries(
        x=unique(first(loops)),
        files=hicFile,
        binSize=10e3
    ) |>
        suppressMessages() |>
        suppressWarnings() |>
        expect_no_error()
})
