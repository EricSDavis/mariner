# Get the pixel representing the strongest or weakest interaction in an InteractionArray

Get the pixel representing the strongest or weakest interaction in an
InteractionArray

## Usage

``` r
selectPixel(
  x,
  aggFUN = sum,
  selectFUN = "which.max",
  nBlocks = 5,
  verbose = TRUE
)

# S4 method for class 'InteractionArray'
selectPixel(
  x,
  aggFUN = sum,
  selectFUN = "which.max",
  nBlocks = 5,
  verbose = TRUE
)
```

## Arguments

- x:

  InteractionArray object

- aggFUN:

  Function to use for aggregating across Hi-C files. Must be passable to
  \`which.max\` or \`which.min\`. Default is "sum".

- selectFUN:

  Function to use for selecting among aggregated interactions. Must be
  one of "which.max" or "which.min".

- nBlocks:

  Number of blocks for block-processing arrays. Default is 5. Increase
  this for large datasets. To read and process all data at once, set
  this value to 1.

- verbose:

  Boolean (TRUE or FALSE) describing whether to report block-processing
  progress. Default is TRUE.

## Value

A GInteractions object with the updated pixel interactions, along with a
column with the aggregated max/min value for that pixel.

## Examples

``` r
## Load marinerData
if (!require("marinerData", quietly = TRUE))
    BiocManager::install("marinerData")

## Read .hic file paths
hicFiles <- c(
    marinerData::LEUK_HEK_PJA27_inter_30.hic(),
    marinerData::LEUK_HEK_PJA30_inter_30.hic()
)
#> see ?marinerData and browseVignettes('marinerData') for documentation
#> loading from cache
#> see ?marinerData and browseVignettes('marinerData') for documentation
#> loading from cache
names(hicFiles) <- c("FS", "WT")

## Read in loops as GInteractions object
loops <-
    WT_5kbLoops.txt() |>
    setNames("WT") |>
    read.table(header=TRUE) |>
    as_ginteractions(keep.extra.columns=FALSE)
#> see ?marinerData and browseVignettes('marinerData') for documentation
#> loading from cache

## Removes the "chr" prefix for compatibility
## with the preprocessed hic files
GenomeInfoDb::seqlevelsStyle(loops) <- 'ENSEMBL'

## Rebin loops to 2.5e6 resolution
loops <- assignToBins(x=loops, binSize=2.5e06)

## Pull 5x5 matrices
iarr <- pullHicMatrices(x=loops[1:5],
                        files=hicFiles,
                        binSize=500e3,
                        norm="KR",
                        half='upper')

## Select pixel
selectPixel(iarr)
#> / Reading and realizing block 1/5 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 2/5 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 3/5 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 4/5 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 5/5 ... 
#> OK
#> \ Processing it ... 
#> OK
#> GInteractions object with 5 interactions and 1 metadata column:
#>       seqnames1             ranges1     seqnames2             ranges2 |
#>           <Rle>           <IRanges>         <Rle>           <IRanges> |
#>   [1]         9   14500000-15000000 ---         9   14500000-15000000 |
#>   [2]         9   89500000-90000000 ---         9   89500000-90000000 |
#>   [3]         9   23500000-24000000 ---         9   23500000-24000000 |
#>   [4]         9 128500000-129000000 ---         9 128500000-129000000 |
#>   [5]         9 113000000-113500000 ---         9 113000000-113500000 |
#>           value
#>       <numeric>
#>   [1]   543.407
#>   [2]   503.798
#>   [3]   464.689
#>   [4]   719.754
#>   [5]   540.652
#>   -------
#>   regions: 5 ranges and 0 metadata columns
#>   seqinfo: 23 sequences from an unspecified genome; no seqlengths
```
