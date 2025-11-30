# Merge sets of paired interactions

Sets of paired range objects (i.e., \`GInteractions\`) are first
clustered by genomic distance with \`dbscan\`, then a representative
interaction is selected for each cluster.

## Usage

``` r
mergePairs(
  x,
  radius,
  method = "manhattan",
  column = NULL,
  selectMax = TRUE,
  pos = "center"
)

# S4 method for class 'list_OR_SimpleList_OR_GInteractions,numeric'
mergePairs(
  x,
  radius,
  method = "manhattan",
  column = NULL,
  selectMax = TRUE,
  pos = "center"
)
```

## Arguments

- x:

  List of \`GInteractions\` or \`data.frame\`-like objects.

- radius:

  Numeric describing the distance in base pairs used to define a cluster
  or pairs.

- method:

  Character describing the distance measure to be used. This must be one
  of "euclidean", "maximum", "manhattan", "canberra", "binary" or
  "minkowski". Any unambiguous substring can be given. Default is
  "manhattan".

- column:

  Character denoting the column to be used to select among clustered
  interactions.

- selectMax:

  Logical. TRUE (default) uses \`which.max()\` to select the interaction
  pair. FALSE uses \`which.min()\`. Only applicable when \`column\` is
  specified.

- pos:

  Positions used for clustering pairs. Must be one of "start", "end" or
  "center". Default is "center".

## Value

Returns a \`MergedGInteractions\` object.

## Details

Interactions are clustered into groups using the provided base pair
\`radius\`, and distance \`method\` with \`dbscan()\`. Representative
interactions are selected for each group by one of two methods. If
\`column\` and \`selectMax\` arguments are provided, the representative
interaction with the maximum (or minimum) value in \`column\` is
returned for each cluster. If these parameters are missing, new ranges
for each pair are returned by calculating the median of modes for each
cluster.

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
    lapply(as_ginteractions)

## Cluster & merge pairs
x <- mergePairs(x = giList,
                radius = 10e03,
                column = "APScoreAvg")
x
#> MergedGInteractions object with 1681 interactions and 9 metadata columns:
#>          seqnames1             ranges1     seqnames2             ranges2 |
#>              <Rle>           <IRanges>         <Rle>           <IRanges> |
#>      [1]      chr9 118645000-118650000 ---      chr9 119330000-119335000 |
#>      [2]      chr9   15280000-15285000 ---      chr9   15405000-15410000 |
#>      [3]      chr9 110180000-110185000 ---      chr9 111520000-111525000 |
#>      [4]      chr9   80375000-80380000 ---      chr9   80650000-80655000 |
#>      [5]      chr9 108380000-108385000 ---      chr9 108475000-108480000 |
#>      ...       ...                 ... ...       ...                 ... .
#>   [1677]      chr7   27235000-27240000 ---      chr7   27815000-27820000 |
#>   [1678]      chr7   86565000-86570000 ---      chr7   86775000-86780000 |
#>   [1679]      chr7   95695000-95700000 ---      chr7   96660000-96665000 |
#>   [1680]      chr7 116620000-116625000 ---      chr7 116760000-116765000 |
#>   [1681]      chr7 140215000-140220000 ---      chr7 140350000-140355000 |
#>                color APScoreAvg ProbabilityofEnrichment RegAPScoreAvg
#>          <character>  <numeric>               <numeric>     <numeric>
#>      [1]       0,0,0    2.61103                0.986044       1.41438
#>      [2]       0,0,0    2.45301                0.982802       1.54370
#>      [3]       0,0,0    3.40635                0.996545       1.80221
#>      [4]       0,0,0    2.09352                0.946897       1.40194
#>      [5]       0,0,0    2.14182                0.936547       1.41005
#>      ...         ...        ...                     ...           ...
#>   [1677]       0,0,0    5.93708                0.999598       3.72133
#>   [1678]       0,0,0    4.60998                0.993894       2.73905
#>   [1679]       0,0,0    5.34541                0.999958       2.47288
#>   [1680]       0,0,0    2.97086                0.969639       1.86991
#>   [1681]       0,0,0    4.52114                0.993302       2.96563
#>          Avg_diffMaxNeihgboor_1 Avg_diffMaxNeihgboor_2       avg       std
#>                       <numeric>              <numeric> <numeric> <numeric>
#>      [1]               1.875585                2.12164   2.60512  0.794220
#>      [2]               1.491568                1.60777   2.73756  0.856876
#>      [3]               3.768934                3.89112   4.49898  1.492786
#>      [4]               0.818514                1.10611   2.20826  0.589976
#>      [5]               0.907892                1.01824   1.95072  0.543293
#>      ...                    ...                    ...       ...       ...
#>   [1677]                2.28945                3.41471   5.78381  1.175767
#>   [1678]                2.07661                2.19518   3.25321  0.810704
#>   [1679]                5.50065                5.70786   5.18389  2.052481
#>   [1680]                1.32308                1.36217   2.31890  0.643926
#>   [1681]                1.04456                2.16908   4.07791  0.663425
#>              value
#>          <numeric>
#>      [1]   4.27231
#>      [2]   4.06339
#>      [3]   7.84915
#>      [4]   2.93583
#>      [5]   2.75773
#>      ...       ...
#>   [1677]   7.81888
#>   [1678]   5.09909
#>   [1679]  10.07336
#>   [1680]   3.49497
#>   [1681]   5.00641
#>   -------
#>   regions: 2812 ranges and 0 metadata columns
#>   seqinfo: 3 sequences from an unspecified genome; no seqlengths
```
