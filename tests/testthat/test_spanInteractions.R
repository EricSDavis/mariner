test_that("spanInteractions works", {
    
    ## Example GInteractions
    gi <- as_ginteractions(read.table(text="
        1 10 40 1 50 80
        2 100 110 2 150 160
        1 10 40 2 50 80
        2 10 40 1 50 80
        1 10 40 1 50 80
    "))
    
    GenomeInfoDb::seqinfo(gi) <- 
        GenomeInfoDb::Seqinfo(
            seqnames=c("1", "2"),
            seqlengths=c(110, 200)
        )
    
    ## Warn about interchromosomal
    spanInteractions(gi) |>
        expect_warning("Inter.*not allowed.*3,4")
    
    ## Warn about inter & oob ranges
    spanInteractions(gi, buffer=2) |>
        expect_warning(".*3,4") |>
        expect_warning(".*1,2,5")
    
    ## Expect first 3
    exp <- GRanges(
        seqnames=c(1,2,1),
        ranges=IRanges(
            start=c(2,93,2),
            end=c(87,166,87)
        ),
        seqinfo=seqinfo(gi)
    )
    suppressWarnings({
        spanInteractions(gi, buffer=0.1) |>
            expect_identical(exp)
    })
    
    ## No warnings
    expect_no_warning(spanInteractions(gi[c(1,2,5)]))
    
})
    