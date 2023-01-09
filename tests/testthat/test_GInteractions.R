library(mariner)

## Shared objects --------------------------------------------------------------

## Create example reference interactions objects
sn1 <- c("chr1", "chr2", "chr1")
s1 <- c(10L, 30L, 50L)
e1 <- c(20L, 40L, 60L)
sn2 <- c("chr1", "chr2", "chr3")
s2 <- c(50L, 60L, 10L)
e2 <- c(60L, 70L, 20L)

gi <- as_ginteractions(data.frame(sn1, s1, e1, sn2, s2, e2))
iset <- InteractionSet(assays=matrix(nrow=3), interactions=gi)
imat <-
    InteractionMatrix(
        interactions=gi,
        assays=list(counts=array(0, dim=c(3,0)))
    )
iarr <-
    InteractionArray(
        interactions=gi,
        assays=list(
            counts=array(0, dim=c(3,0,0,0)),
            rownames=array(0, dim=c(3,0,0,0)),
            colnames=array(0, dim=c(3,0,0,0))
        )
    )

test_that("Accessors for GInteraction-like objects", {

    ## seqnames1
    expect_identical(seqnames1(gi), sn1)
    expect_identical(seqnames1(iset), sn1)
    expect_identical(seqnames1(imat), sn1)
    expect_identical(seqnames1(iarr), sn1)

    ## start1
    expect_identical(start1(gi), s1)
    expect_identical(start1(iset), s1)
    expect_identical(start1(imat), s1)
    expect_identical(start1(iarr), s1)

    ## end1
    expect_identical(end1(gi), e1)
    expect_identical(end1(iset), e1)
    expect_identical(end1(imat), e1)
    expect_identical(end1(iarr), e1)

    ## seqnames2
    expect_identical(seqnames2(gi), sn2)
    expect_identical(seqnames2(iset), sn2)
    expect_identical(seqnames2(imat), sn2)
    expect_identical(seqnames2(iarr), sn2)

    ## start2
    expect_identical(start2(gi), s2)
    expect_identical(start2(iset), s2)
    expect_identical(start2(imat), s2)
    expect_identical(start2(iarr), s2)

    ## end2
    expect_identical(end2(gi), e2)
    expect_identical(end2(iset), e2)
    expect_identical(end2(imat), e2)
    expect_identical(end2(iarr), e2)

})
