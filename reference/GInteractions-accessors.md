# Access each portion of a GInteractions-like object

Access each portion of a GInteractions-like object

## Usage

``` r
seqnames1(x, ...)

seqnames2(x, ...)

start1(x, ...)

end1(x, ...)

start2(x, ...)

end2(x, ...)

# S4 method for class 'GInteractions_OR_InteractionSet'
seqnames1(x)

# S4 method for class 'GInteractions_OR_InteractionSet'
seqnames2(x)

# S4 method for class 'GInteractions_OR_InteractionSet'
start1(x)

# S4 method for class 'GInteractions_OR_InteractionSet'
end1(x)

# S4 method for class 'GInteractions_OR_InteractionSet'
start2(x)

# S4 method for class 'GInteractions_OR_InteractionSet'
end2(x)
```

## Arguments

- x:

  GInteractions object.

- ...:

  Additional arguments.

## Value

A vector of values corresponding to the requested component of a
GInteractions-like object. For seqnames1 and seqnames2 the RLE is
coerced to a character vector.

## Examples

``` r
library(InteractionSet)
#> Loading required package: GenomicRanges
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: ‘generics’
#> The following objects are masked from ‘package:base’:
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: ‘BiocGenerics’
#> The following objects are masked from ‘package:stats’:
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from ‘package:base’:
#> 
#>     Filter, Find, Map, Position, Reduce, anyDuplicated, aperm, append,
#>     as.data.frame, basename, cbind, colnames, dirname, do.call,
#>     duplicated, eval, evalq, get, grep, grepl, is.unsorted, lapply,
#>     mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     rank, rbind, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: ‘S4Vectors’
#> The following object is masked from ‘package:utils’:
#> 
#>     findMatches
#> The following objects are masked from ‘package:base’:
#> 
#>     I, expand.grid, unname
#> Loading required package: IRanges
#> Loading required package: Seqinfo
#> Loading required package: SummarizedExperiment
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: ‘MatrixGenerics’
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#>     colWeightedMeans, colWeightedMedians, colWeightedSds,
#>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#>     rowWeightedSds, rowWeightedVars
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: ‘Biobase’
#> The following object is masked from ‘package:MatrixGenerics’:
#> 
#>     rowMedians
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     anyMissing, rowMedians
## Create example reference interactions objects
gi <- read.table(text="
    chr1 10 20 chr1 50 60
    chr2 30 40 chr2 60 70
    chr1 50 60 chr3 10 20") |>
    as_ginteractions()

iset <- InteractionSet(assays=matrix(nrow=3),
                       interactions=gi)

## Access vectors of values
seqnames1(gi)
#> [1] "chr1" "chr2" "chr1"
start1(gi)
#> [1] 10 30 50
end1(gi)
#> [1] 20 40 60
seqnames2(gi)
#> [1] "chr1" "chr2" "chr3"
start2(gi)
#> [1] 50 60 10
end2(gi)
#> [1] 60 70 20

## Also works for InteractionSet-like objects
seqnames1(iset)
#> [1] "chr1" "chr2" "chr1"
start1(iset)
#> [1] 10 30 50
end1(iset)
#> [1] 20 40 60
seqnames2(iset)
#> [1] "chr1" "chr2" "chr3"
start2(iset)
#> [1] 50 60 10
end2(iset)
#> [1] 60 70 20
```
