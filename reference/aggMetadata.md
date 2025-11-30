# Aggregate the metadata columns of merged pairs

Aggregate the metadata columns of merged pairs

## Usage

``` r
aggMetadata(x, columns, funs)

# S4 method for class 'MergedGInteractions,character,character_OR_function_OR_list'
aggMetadata(x, columns, funs)
```

## Arguments

- x:

  MergedGInteractions object.

- columns:

  Character vector of columns to aggregate.

- funs:

  Character vector of functions to apply to \`columns\`.

## Value

\`x\` with aggregated metadata columns

## Examples

``` r
## Load marinerData
if (!require("marinerData", quietly = TRUE))
    BiocManager::install("marinerData")

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
    lapply(bedpeFiles, read.table, header=TRUE, nrows=1000) |>
    lapply(as_ginteractions) |>
    setNames(gsub("^.*extdata/(.{2}).*$", "\\1", bedpeFiles))

## Add names describing the source and loop
giList <- lapply(seq_along(giList), \(i) {
    x <- giList[[i]]
    x$name <- paste0(names(giList)[i], "_loop_", length(x))
    return(x)
})

## Cluster & merge pairs
x <- mergePairs(x = giList,
                radius = 5e03)

## List loop names
aggMetadata(x, columns = "name", fun = "list")
#> MergedGInteractions object with 1786 interactions and 1 metadata column:
#>          seqnames1             ranges1     seqnames2             ranges2 |
#>              <Rle>           <IRanges>         <Rle>           <IRanges> |
#>      [1]      chr9 118645000-118650000 ---      chr9 119330000-119335000 |
#>      [2]      chr9   15280000-15285000 ---      chr9   15405000-15410000 |
#>      [3]      chr9   17815000-17820000 ---      chr9   18205000-18210000 |
#>      [4]      chr9 110180000-110185000 ---      chr9 111520000-111525000 |
#>      [5]      chr9   80375000-80380000 ---      chr9   80650000-80655000 |
#>      ...       ...                 ... ...       ...                 ... .
#>   [1782]      chr7   33795000-33800000 ---      chr7   33932500-33937500 |
#>   [1783]      chr7   25295000-25300000 ---      chr7   25887500-25892500 |
#>   [1784]      chr7   86567500-86572500 ---      chr7   86775000-86780000 |
#>   [1785]      chr7 116622500-116627500 ---      chr7 116760000-116765000 |
#>   [1786]      chr7 140215000-140220000 ---      chr7 140347500-140352500 |
#>                                              list.name
#>                                                 <list>
#>      [1]                        /home/runner/.cache/..
#>      [2]                        /home/runner/.cache/..
#>      [3]                        /home/runner/.cache/..
#>      [4]                        /home/runner/.cache/..
#>      [5]                        /home/runner/.cache/..
#>      ...                                           ...
#>   [1782] /home/runner/.cache/..,/home/runner/.cache/..
#>   [1783] /home/runner/.cache/..,/home/runner/.cache/..
#>   [1784] /home/runner/.cache/..,/home/runner/.cache/..
#>   [1785] /home/runner/.cache/..,/home/runner/.cache/..
#>   [1786] /home/runner/.cache/..,/home/runner/.cache/..
#>   -------
#>   regions: 2962 ranges and 0 metadata columns
#>   seqinfo: 3 sequences from an unspecified genome; no seqlengths

## Aggregate values
aggMetadata(x, columns = c("APScoreAvg"), fun = "mean")
#> MergedGInteractions object with 1786 interactions and 1 metadata column:
#>          seqnames1             ranges1     seqnames2             ranges2 |
#>              <Rle>           <IRanges>         <Rle>           <IRanges> |
#>      [1]      chr9 118645000-118650000 ---      chr9 119330000-119335000 |
#>      [2]      chr9   15280000-15285000 ---      chr9   15405000-15410000 |
#>      [3]      chr9   17815000-17820000 ---      chr9   18205000-18210000 |
#>      [4]      chr9 110180000-110185000 ---      chr9 111520000-111525000 |
#>      [5]      chr9   80375000-80380000 ---      chr9   80650000-80655000 |
#>      ...       ...                 ... ...       ...                 ... .
#>   [1782]      chr7   33795000-33800000 ---      chr7   33932500-33937500 |
#>   [1783]      chr7   25295000-25300000 ---      chr7   25887500-25892500 |
#>   [1784]      chr7   86567500-86572500 ---      chr7   86775000-86780000 |
#>   [1785]      chr7 116622500-116627500 ---      chr7 116760000-116765000 |
#>   [1786]      chr7 140215000-140220000 ---      chr7 140347500-140352500 |
#>          mean.APScoreAvg
#>                <numeric>
#>      [1]         2.61103
#>      [2]         2.45301
#>      [3]         4.73751
#>      [4]         3.40635
#>      [5]         2.09352
#>      ...             ...
#>   [1782]         2.89949
#>   [1783]         5.40040
#>   [1784]         4.60903
#>   [1785]         2.86398
#>   [1786]         4.48766
#>   -------
#>   regions: 2962 ranges and 0 metadata columns
#>   seqinfo: 3 sequences from an unspecified genome; no seqlengths
aggMetadata(x, columns = c("APScoreAvg", "avg"), fun = "mean")
#> MergedGInteractions object with 1786 interactions and 2 metadata columns:
#>          seqnames1             ranges1     seqnames2             ranges2 |
#>              <Rle>           <IRanges>         <Rle>           <IRanges> |
#>      [1]      chr9 118645000-118650000 ---      chr9 119330000-119335000 |
#>      [2]      chr9   15280000-15285000 ---      chr9   15405000-15410000 |
#>      [3]      chr9   17815000-17820000 ---      chr9   18205000-18210000 |
#>      [4]      chr9 110180000-110185000 ---      chr9 111520000-111525000 |
#>      [5]      chr9   80375000-80380000 ---      chr9   80650000-80655000 |
#>      ...       ...                 ... ...       ...                 ... .
#>   [1782]      chr7   33795000-33800000 ---      chr7   33932500-33937500 |
#>   [1783]      chr7   25295000-25300000 ---      chr7   25887500-25892500 |
#>   [1784]      chr7   86567500-86572500 ---      chr7   86775000-86780000 |
#>   [1785]      chr7 116622500-116627500 ---      chr7 116760000-116765000 |
#>   [1786]      chr7 140215000-140220000 ---      chr7 140347500-140352500 |
#>          mean.APScoreAvg  mean.avg
#>                <numeric> <numeric>
#>      [1]         2.61103   2.60512
#>      [2]         2.45301   2.73756
#>      [3]         4.73751   4.68239
#>      [4]         3.40635   4.49898
#>      [5]         2.09352   2.20826
#>      ...             ...       ...
#>   [1782]         2.89949   2.08193
#>   [1783]         5.40040   3.58636
#>   [1784]         4.60903   3.34010
#>   [1785]         2.86398   2.43657
#>   [1786]         4.48766   3.98591
#>   -------
#>   regions: 2962 ranges and 0 metadata columns
#>   seqinfo: 3 sequences from an unspecified genome; no seqlengths
aggMetadata(x, columns = c("APScoreAvg"), fun = c("mean", "median"))
#> MergedGInteractions object with 1786 interactions and 2 metadata columns:
#>          seqnames1             ranges1     seqnames2             ranges2 |
#>              <Rle>           <IRanges>         <Rle>           <IRanges> |
#>      [1]      chr9 118645000-118650000 ---      chr9 119330000-119335000 |
#>      [2]      chr9   15280000-15285000 ---      chr9   15405000-15410000 |
#>      [3]      chr9   17815000-17820000 ---      chr9   18205000-18210000 |
#>      [4]      chr9 110180000-110185000 ---      chr9 111520000-111525000 |
#>      [5]      chr9   80375000-80380000 ---      chr9   80650000-80655000 |
#>      ...       ...                 ... ...       ...                 ... .
#>   [1782]      chr7   33795000-33800000 ---      chr7   33932500-33937500 |
#>   [1783]      chr7   25295000-25300000 ---      chr7   25887500-25892500 |
#>   [1784]      chr7   86567500-86572500 ---      chr7   86775000-86780000 |
#>   [1785]      chr7 116622500-116627500 ---      chr7 116760000-116765000 |
#>   [1786]      chr7 140215000-140220000 ---      chr7 140347500-140352500 |
#>          mean.APScoreAvg median.APScoreAvg
#>                <numeric>         <numeric>
#>      [1]         2.61103           2.61103
#>      [2]         2.45301           2.45301
#>      [3]         4.73751           4.73751
#>      [4]         3.40635           3.40635
#>      [5]         2.09352           2.09352
#>      ...             ...               ...
#>   [1782]         2.89949           2.89949
#>   [1783]         5.40040           5.40040
#>   [1784]         4.60903           4.60903
#>   [1785]         2.86398           2.86398
#>   [1786]         4.48766           4.48766
#>   -------
#>   regions: 2962 ranges and 0 metadata columns
#>   seqinfo: 3 sequences from an unspecified genome; no seqlengths

## Custom functions
aggMetadata(x, columns = c("APScoreAvg"), fun = \(x) {
    ifelse(is.na(sd(x)), 0, sd(x))
})
#> MergedGInteractions object with 1786 interactions and 1 metadata column:
#>          seqnames1             ranges1     seqnames2             ranges2 |
#>              <Rle>           <IRanges>         <Rle>           <IRanges> |
#>      [1]      chr9 118645000-118650000 ---      chr9 119330000-119335000 |
#>      [2]      chr9   15280000-15285000 ---      chr9   15405000-15410000 |
#>      [3]      chr9   17815000-17820000 ---      chr9   18205000-18210000 |
#>      [4]      chr9 110180000-110185000 ---      chr9 111520000-111525000 |
#>      [5]      chr9   80375000-80380000 ---      chr9   80650000-80655000 |
#>      ...       ...                 ... ...       ...                 ... .
#>   [1782]      chr7   33795000-33800000 ---      chr7   33932500-33937500 |
#>   [1783]      chr7   25295000-25300000 ---      chr7   25887500-25892500 |
#>   [1784]      chr7   86567500-86572500 ---      chr7   86775000-86780000 |
#>   [1785]      chr7 116622500-116627500 ---      chr7 116760000-116765000 |
#>   [1786]      chr7 140215000-140220000 ---      chr7 140347500-140352500 |
#>          fun1.APScoreAvg
#>                <numeric>
#>      [1]               0
#>      [2]               0
#>      [3]               0
#>      [4]               0
#>      [5]               0
#>      ...             ...
#>   [1782]      0.03410843
#>   [1783]      0.96721962
#>   [1784]      0.00133997
#>   [1785]      0.15114917
#>   [1786]      0.04735593
#>   -------
#>   regions: 2962 ranges and 0 metadata columns
#>   seqinfo: 3 sequences from an unspecified genome; no seqlengths
```
