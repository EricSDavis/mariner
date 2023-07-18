test_that("pileupDomains works", {
    
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
    pileupDomains(x=loops, files=hicFile, binSize=10e3, buffer=0,
              verbose=FALSE) |>
        expect_message("Caution") |>
        suppressMessages() |>
        suppressWarnings()
})
