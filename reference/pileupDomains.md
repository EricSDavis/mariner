# Pileup Hi-C domains

pileupDomains expands then extracts regions/domains from Hi-C files,
regularizes them so they are the same size, then aggregates them into a
single matrix. This can be used to perform aggregate TAD analysis.

## Usage

``` r
pileupDomains(
  x,
  files,
  binSize,
  buffer = 0.5,
  ndim = c(100, 100),
  scale = TRUE,
  normalize = TRUE,
  FUN = sum,
  nBlocks = 50,
  verbose = TRUE,
  BPPARAM = bpparam(),
  blockSize = 1e+06,
  ...
)

# S4 method for class 'GRanges_OR_GInteractions,character,numeric'
pileupDomains(
  x,
  files,
  binSize,
  buffer = 0.5,
  ndim = c(100, 100),
  scale = TRUE,
  normalize = TRUE,
  FUN = sum,
  nBlocks = 50,
  verbose = TRUE,
  BPPARAM = bpparam(),
  blockSize = 1e+06,
  ...
)
```

## Arguments

- x:

  GRanges or GInteractions object containing the TADs or Loops to be
  aggregated.

- files:

  Character file paths to \`.hic\` files.

- binSize:

  Integer (numeric) describing the resolution (range widths) of the
  paired data. Note that small values for this argument may lead to R
  session crashes.

- buffer:

  Fraction (length one numeric vector) pair-distance to expand around
  the resulting range.

- ndim:

  Numeric vector of length two describing the new dimensions of the
  output matrices.

- scale:

  Boolean (TRUE/FALSE) indicating whether the values in the new matrices
  should be scaled to the total signal in each matrix.

- normalize:

  Boolean, whether to normalize the aggregated values to the number of
  interactions.

- FUN:

  Function to use for aggregating.

- nBlocks:

  Number of blocks for block-processing arrays. Default is 50. Increase
  this for large datasets. To read and process all data at once, set
  this value to 1.

- verbose:

  Boolean (TRUE or FALSE) describing whether to report block-processing
  progress.

- BPPARAM:

  Parallelization params (passed to \`BiocParallel::bplapply()\`).
  Default is the result of \`BiocParallel::bpparams()\`. Parallel
  processing is not available when \`by=interactions\`.

- blockSize:

  Number (length one numeric vector) describing the size in base-pairs
  to pull from each \`.hic\` file. Default is 1e6. For large \`.hic\`
  files \`blockSize\` can be reduced to conserve the amount of data read
  in at a time. Larger \`blockSize\` values speed up performance, but
  use more memory.

- ...:

  Additional arguments passed to \`pullHicMatrices()\`.

## Value

A DelayedArray of aggregated counts.

## Details

It may be necessary to adjust the \`zrange\` in \`plotMatrix()\` since
the Hi-C diagonal will dominate the scale.

If interactions are passed to the function, only intrachromosomal ranges
are maintained.

Using small \`binSize\` values with large ranges may lead to pulling
very large sections of a Hi-C map that can crash your R session. If this
happens try increasing the \`binSize\` and \`nBlocks\` parameters, while
lower the \`blockSize\` parameter.

## Examples

``` r
## Load marinerData
if (!require("marinerData", quietly = TRUE))
    BiocManager::install("marinerData")

## Read .hic file paths
hicFile <- marinerData::LEUK_HEK_PJA27_inter_30.hic()
#> see ?marinerData and browseVignettes('marinerData') for documentation
#> loading from cache
names(hicFile) <- "FS"

## Loops
loops <- 
    marinerData::FS_5kbLoops.txt() |>
    read.table(header=TRUE, nrows=100) |>
    as_ginteractions() |>
    GenomeInfoDb::`seqlevelsStyle<-`(value='ENSEMBL')
#> see ?marinerData and browseVignettes('marinerData') for documentation
#> loading from cache

## Warn about small binSize
pileupDomains(x=loops, files=hicFile, binSize=50e3, buffer=0.25)
#> Pairs are not binned to `binSize`.
#> ℹ Snapping to binSize=50000, 
#> ℹ Use `assignToBins()` or `snapToBins()` for more binning options.
#> Warning: ✖ `blockSize` must be twice the longest range.
#> ℹ Setting `blockSize` = 5200002.
#> / Reading and realizing block 1/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> Loading required namespace: fields
#> OK
#> / Reading and realizing block 2/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 3/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 4/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 5/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 6/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 7/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 8/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 9/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 10/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 11/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 12/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 13/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 14/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 15/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 16/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 17/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 18/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 19/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 20/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 21/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 22/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 23/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 24/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 25/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 26/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 27/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 28/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 29/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 30/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 31/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 32/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 33/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 34/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 35/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 36/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 37/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 38/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 39/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 40/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 41/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 42/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 43/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 44/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 45/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 46/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 47/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 48/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 49/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 50/50 of file 1/1 ... 
#> OK
#> \ Processing it ... 
#> OK
#> <100 x 100> DelayedMatrix object of type "double":
#>                [,1]         [,2]         [,3] ...        [,99]       [,100]
#>   [1,] 0.0004685484 0.0004085836 0.0003488673   . 2.171695e-05 2.067000e-05
#>   [2,] 0.0004085836 0.0003805260 0.0003523911   . 2.214473e-05 2.150746e-05
#>   [3,] 0.0003488673 0.0003523911 0.0003555542   . 2.256156e-05 2.232234e-05
#>   [4,] 0.0002986538 0.0003217905 0.0003449622   . 2.288477e-05 2.277575e-05
#>   [5,] 0.0002618003 0.0002887709 0.0003161078   . 2.342872e-05 2.333199e-05
#>    ...            .            .            .   .            .            .
#>  [96,] 2.277187e-05 2.334181e-05 2.391175e-05   . 0.0002714690 0.0002505143
#>  [97,] 2.326133e-05 2.334866e-05 2.343599e-05   . 0.0003023433 0.0002843798
#>  [98,] 2.276390e-05 2.279363e-05 2.282336e-05   . 0.0003302839 0.0003293499
#>  [99,] 2.171695e-05 2.214473e-05 2.256156e-05   . 0.0003538560 0.0003773022
#> [100,] 2.067000e-05 2.150746e-05 2.232234e-05   . 0.0003773022 0.0004251867
```
