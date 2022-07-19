library(hictoolsr)
library(data.table)
library(purrr)
library(glue)

## Define paths to loop files
loopFiles <-
    system.file("extdata", package = "hictoolsr") |>
    list.files(pattern = "5kbLoops.txt", full.names = TRUE)

## Read in loop files & combine
loops <-
    lapply(loopFiles, fread) |>
    do.call(rbind, args = _) |>
    as_ginteractions()


# "../test/custom.hic" |>
#     strawr::straw(norm = "KR",
#                   fname = _,
#                   chr1loc = "WT_9:14460000:14765000",
#                   chr2loc = "WT_9:14460000:14765000",
#                   unit = "BP",
#                   binsize = 100e03,
#                   matrix = "observed")


# "../test/custom.hic" |>
#     strawr::straw(norm = "NONE",
#                   fname = _,
#                   chr1loc = "WT_9:9740000:19740000",
#                   chr2loc = "WT_9:9740000:9740000",
#                   unit = "BP",
#                   binsize = 5e03,
#                   matrix = "observed")


library(InteractionSet)

## Define resolution and buffer
res <- 5000
buffer <- 20

## Replace with expanded regions
regions(loops) <- (regions(loops) + res*buffer)
seqlevelsStyle(loops) <- "Ensembl"

## Convert back to data.table
## and put into "short format"
## https://github.com/aidenlab/juicer/wiki/Pre#short-format
loops |>
    as.data.table() |>
    {\(x) x[,c(1:3, 6:8)]}()

loops |>
    lapply(as.data.table()) |>
    do.call(rbind, args = _)
    {\(x) x[,c(1:3, 6:8)]}()

loops |>
    lapply(as.data.table) |>
    do.call(rbind, args = _)


width(loops)
