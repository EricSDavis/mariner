# Get clustered pairs from MergedGInteractions object

Returns the clustered pairs associated with each range in the
\`MergedGInteractions\` object. Order always follows the indices of the
\`MergedGInteractions\` object.

## Usage

``` r
clusters(x, ...)

# S4 method for class 'MergedGInteractions'
clusters(x)
```

## Arguments

- x:

  MergedGInteractions object.

- ...:

  Additional arguments.

## Value

A list of data.tables cooresponding to each pair in \`x\`.

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
    lapply(bedpeFiles, fread, nrows = 1000) |>
    lapply(as_ginteractions)

## Cluster & merge pairs
x <- mergePairs(x = giList,
                radius = 10e03,
                column = "APScoreAvg")

## Access pair clusters
clusters(x[1:3])
#> [[1]]
#>    seqnames1    start1      end1 width1 strand1 seqnames2    start2      end2
#>       <fctr>     <int>     <int>  <int>  <fctr>    <fctr>     <int>     <int>
#> 1:      chr9 118645000 118650000   5001       *      chr9 119330000 119335000
#>    width2 strand2  color APScoreAvg ProbabilityofEnrichment RegAPScoreAvg
#>     <int>  <fctr> <char>      <num>                   <num>         <num>
#> 1:   5001       *  0,0,0   2.611034               0.9860443      1.414376
#>    Avg_diffMaxNeihgboor_1 Avg_diffMaxNeihgboor_2      avg       std    value
#>                     <num>                  <num>    <num>     <num>    <num>
#> 1:               1.875585               2.121637 2.605124 0.7942202 4.272311
#>       src
#>    <char>
#> 1:     FS
#> 
#> [[2]]
#>    seqnames1   start1     end1 width1 strand1 seqnames2   start2     end2
#>       <fctr>    <int>    <int>  <int>  <fctr>    <fctr>    <int>    <int>
#> 1:      chr9 15280000 15285000   5001       *      chr9 15405000 15410000
#>    width2 strand2  color APScoreAvg ProbabilityofEnrichment RegAPScoreAvg
#>     <int>  <fctr> <char>      <num>                   <num>         <num>
#> 1:   5001       *  0,0,0   2.453013               0.9828022      1.543698
#>    Avg_diffMaxNeihgboor_1 Avg_diffMaxNeihgboor_2      avg      std    value
#>                     <num>                  <num>    <num>    <num>    <num>
#> 1:               1.491568               1.607766 2.737556 0.856876 4.063394
#>       src
#>    <char>
#> 1:     FS
#> 
#> [[3]]
#>    seqnames1    start1      end1 width1 strand1 seqnames2    start2      end2
#>       <fctr>     <int>     <int>  <int>  <fctr>    <fctr>     <int>     <int>
#> 1:      chr9 110180000 110185000   5001       *      chr9 111520000 111525000
#>    width2 strand2  color APScoreAvg ProbabilityofEnrichment RegAPScoreAvg
#>     <int>  <fctr> <char>      <num>                   <num>         <num>
#> 1:   5001       *  0,0,0   3.406346               0.9965453      1.802207
#>    Avg_diffMaxNeihgboor_1 Avg_diffMaxNeihgboor_2      avg      std    value
#>                     <num>                  <num>    <num>    <num>    <num>
#> 1:               3.768934               3.891115 4.498982 1.492786 7.849146
#>       src
#>    <char>
#> 1:     FS
#> 
clusters(x[3:1])
#> [[1]]
#>    seqnames1    start1      end1 width1 strand1 seqnames2    start2      end2
#>       <fctr>     <int>     <int>  <int>  <fctr>    <fctr>     <int>     <int>
#> 1:      chr9 110180000 110185000   5001       *      chr9 111520000 111525000
#>    width2 strand2  color APScoreAvg ProbabilityofEnrichment RegAPScoreAvg
#>     <int>  <fctr> <char>      <num>                   <num>         <num>
#> 1:   5001       *  0,0,0   3.406346               0.9965453      1.802207
#>    Avg_diffMaxNeihgboor_1 Avg_diffMaxNeihgboor_2      avg      std    value
#>                     <num>                  <num>    <num>    <num>    <num>
#> 1:               3.768934               3.891115 4.498982 1.492786 7.849146
#>       src
#>    <char>
#> 1:     FS
#> 
#> [[2]]
#>    seqnames1   start1     end1 width1 strand1 seqnames2   start2     end2
#>       <fctr>    <int>    <int>  <int>  <fctr>    <fctr>    <int>    <int>
#> 1:      chr9 15280000 15285000   5001       *      chr9 15405000 15410000
#>    width2 strand2  color APScoreAvg ProbabilityofEnrichment RegAPScoreAvg
#>     <int>  <fctr> <char>      <num>                   <num>         <num>
#> 1:   5001       *  0,0,0   2.453013               0.9828022      1.543698
#>    Avg_diffMaxNeihgboor_1 Avg_diffMaxNeihgboor_2      avg      std    value
#>                     <num>                  <num>    <num>    <num>    <num>
#> 1:               1.491568               1.607766 2.737556 0.856876 4.063394
#>       src
#>    <char>
#> 1:     FS
#> 
#> [[3]]
#>    seqnames1    start1      end1 width1 strand1 seqnames2    start2      end2
#>       <fctr>     <int>     <int>  <int>  <fctr>    <fctr>     <int>     <int>
#> 1:      chr9 118645000 118650000   5001       *      chr9 119330000 119335000
#>    width2 strand2  color APScoreAvg ProbabilityofEnrichment RegAPScoreAvg
#>     <int>  <fctr> <char>      <num>                   <num>         <num>
#> 1:   5001       *  0,0,0   2.611034               0.9860443      1.414376
#>    Avg_diffMaxNeihgboor_1 Avg_diffMaxNeihgboor_2      avg       std    value
#>                     <num>                  <num>    <num>     <num>    <num>
#> 1:               1.875585               2.121637 2.605124 0.7942202 4.272311
#>       src
#>    <char>
#> 1:     FS
#> 
clusters(x[c(3, 1, 2)])
#> [[1]]
#>    seqnames1    start1      end1 width1 strand1 seqnames2    start2      end2
#>       <fctr>     <int>     <int>  <int>  <fctr>    <fctr>     <int>     <int>
#> 1:      chr9 110180000 110185000   5001       *      chr9 111520000 111525000
#>    width2 strand2  color APScoreAvg ProbabilityofEnrichment RegAPScoreAvg
#>     <int>  <fctr> <char>      <num>                   <num>         <num>
#> 1:   5001       *  0,0,0   3.406346               0.9965453      1.802207
#>    Avg_diffMaxNeihgboor_1 Avg_diffMaxNeihgboor_2      avg      std    value
#>                     <num>                  <num>    <num>    <num>    <num>
#> 1:               3.768934               3.891115 4.498982 1.492786 7.849146
#>       src
#>    <char>
#> 1:     FS
#> 
#> [[2]]
#>    seqnames1    start1      end1 width1 strand1 seqnames2    start2      end2
#>       <fctr>     <int>     <int>  <int>  <fctr>    <fctr>     <int>     <int>
#> 1:      chr9 118645000 118650000   5001       *      chr9 119330000 119335000
#>    width2 strand2  color APScoreAvg ProbabilityofEnrichment RegAPScoreAvg
#>     <int>  <fctr> <char>      <num>                   <num>         <num>
#> 1:   5001       *  0,0,0   2.611034               0.9860443      1.414376
#>    Avg_diffMaxNeihgboor_1 Avg_diffMaxNeihgboor_2      avg       std    value
#>                     <num>                  <num>    <num>     <num>    <num>
#> 1:               1.875585               2.121637 2.605124 0.7942202 4.272311
#>       src
#>    <char>
#> 1:     FS
#> 
#> [[3]]
#>    seqnames1   start1     end1 width1 strand1 seqnames2   start2     end2
#>       <fctr>    <int>    <int>  <int>  <fctr>    <fctr>    <int>    <int>
#> 1:      chr9 15280000 15285000   5001       *      chr9 15405000 15410000
#>    width2 strand2  color APScoreAvg ProbabilityofEnrichment RegAPScoreAvg
#>     <int>  <fctr> <char>      <num>                   <num>         <num>
#> 1:   5001       *  0,0,0   2.453013               0.9828022      1.543698
#>    Avg_diffMaxNeihgboor_1 Avg_diffMaxNeihgboor_2      avg      std    value
#>                     <num>                  <num>    <num>    <num>    <num>
#> 1:               1.491568               1.607766 2.737556 0.856876 4.063394
#>       src
#>    <char>
#> 1:     FS
#> 
clusters(x) |> length()
#> [1] 1681
```
