# Get selectionMethod from MergedGInteractions object

Get selectionMethod from MergedGInteractions object

## Usage

``` r
selectionMethod(x, ...)

# S4 method for class 'MergedGInteractions'
selectionMethod(x, ...)
```

## Arguments

- x:

  MergedGInteractions object.

- ...:

  Additional arguments.

## Value

A character vector describing which selection method was used for
merging.

## Examples

``` r
## Load required packages
library(data.table, include.only="fread")

## Load marinerData
if (!require("marinerData", quietly = TRUE))
    BiocManager::install("marinerData")

## Reference BEDPE files (loops called with SIP)
bedpeFiles <- c(
    marinerData::FS_5kbLoops.txt(),
    marinerData::WT_5kbLoops.txt()
)
#> see ?marinerData and browseVignettes('marinerData') for documentation
#> loading from cache
#> see ?marinerData and browseVignettes('marinerData') for documentation
#> loading from cache
names(bedpeFiles) <- c("FS", "WT")

## Read in bedpeFiles as a list of GInteractions
## Use only first 1000 rows for fast example
giList <-
    lapply(bedpeFiles, fread, nrows=1000) |>
    lapply(as_ginteractions)

## Cluster & merge pairs
x <- mergePairs(x = giList,
                radius = 10e03,
                column = "APScoreAvg")

selectionMethod(x)
#> Select by column 'APScoreAvg'
```
