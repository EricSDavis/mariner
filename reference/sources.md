# Accessor for sources

Access the names or source files of a \`MergedGInteractions\` object.

## Usage

``` r
sources(x)

# S4 method for class 'MergedGInteractions'
sources(x)
```

## Arguments

- x:

  MergedGInteractions object.

## Value

A character vector of names or source files of a \`MergedGInteractions\`
object.

## Examples

``` r
## Load required packages
library(data.table, include.only="fread")

## Load marinerData
if (!require("marinerData", quietly = TRUE))
    BiocManager::install("marinerData")

## Reference BEDPE files (loops called with SIP)
loopFiles <- c(
    marinerData::FS_5kbLoops.txt(),
    marinerData::WT_5kbLoops.txt()
)
#> see ?marinerData and browseVignettes('marinerData') for documentation
#> loading from cache
#> see ?marinerData and browseVignettes('marinerData') for documentation
#> loading from cache
names(loopFiles) <- c("FS", "WT")

## Read in loopFiles as a list of GInteractions
## Use only first 1000 rows for fast example
giList <-
    lapply(loopFiles, fread, nrows=1000) |>
    lapply(as_ginteractions)

## Cluster & merge pairs
x <- mergePairs(x = giList,
                radius = 10e03)

sources(x)
#> [1] "FS" "WT"
```
