hicFile <- "inst/extdata/test.hic"

strawr::readHicNormTypes(hicFile)
strawr::readHicChroms(hicFile)
strawr::readHicBpResolutions(hicFile)

strawr::straw(norm = "NONE",
              fname = hicFile,
              chr1loc = "1",
              chr2loc = "1",
              unit = "BP",
              binsize = 2.5e06,
              matrix = "observed") |> head()


## Classic straw
chr1locs <- "1:2500000:2500000"
chr2locs <- "1:5000000:5000000"
strawr::straw(norm = "NONE",
              fname = hicFile,
              chr1loc = chr1locs,
              chr2loc = chr2locs,
              unit = "BP",
              binsize = 2.5e06,
              matrix = "observed")

## Vectorized straw
strawv <- Vectorize(strawr::straw, vectorize.args = c("chr1loc", "chr2loc"))

chr1locs <- rep(c("1:2500000:2500000"), 22000)
chr2locs <- rep(c("1:5000000:5000000"), 22000)
system.time({
    result1 <-
        strawv(norm = "NONE",
               fname = hicFile,
               chr1loc = chr1locs,
               chr2loc = chr2locs,
               unit = "BP",
               binsize = 2.5e06,
               matrix = "observed") |> t() |> unname()
})

## Compare with for-loop in R
system.time({
    result2 <-
        lapply(seq_along(chr1locs), \(i){
            strawr::straw(norm = "NONE",
                          fname = hicFile,
                          chr1loc = chr1locs[i],
                          chr2loc = chr2locs[i],
                          unit = "BP",
                          binsize = 2.5e06,
                          matrix = "observed")
        }) |> do.call(rbind, args = _)
})

## Compare with hictoolsr method
library(InteractionSet)
chr1locs <- gsub("(.*):", "\\1-", chr1locs)
chr2locs <- gsub("(.*):", "\\1-", chr2locs)
gi <- GInteractions(anchor1 = GRanges(paste0('chr', chr1locs)),
                    anchor2 = GRanges(paste0('chr', chr2locs))) |>
    hictoolsr::binBedpe(res = 2.5e06,
                        a1Pos = 'start',
                        a2Pos = 'start')

chroms <- grep("ALL",
               value = TRUE,
               invert = TRUE,
               strawr::readHicChroms(hicFile)$name)
system.time({
    result3 <-
        hictoolsr::extractCounts(bedpe = gi[1],
                                 chroms = chroms,
                                 hic = hicFile,
                                 res = 2.5e06,
                                 norm = "NONE",
                                 matrix = "observed")
})
