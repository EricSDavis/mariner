# Pileup Hi-C contacts around boundary regions

pileupBoundaries expands input loci to the specified \`width\`, extracts
then aggregates them into a single matrix. This can be used to aggregate
windows of interactions centered on a set of loci.

## Usage

``` r
pileupBoundaries(
  x,
  files,
  binSize,
  width = 5e+05,
  normalize = TRUE,
  FUN = sum,
  nBlocks = 50,
  verbose = TRUE,
  BPPARAM = bpparam(),
  blockSize = 1e+06,
  ...
)

# S4 method for class 'GRanges_OR_GInteractions,character,numeric'
pileupBoundaries(
  x,
  files,
  binSize,
  width = 5e+05,
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

  GRanges or GInteractions object containing the loci to be aggregated.
  GInteractions will be split into unique anchors.

- files:

  Character file paths to \`.hic\` files.

- binSize:

  Integer (numeric) describing the resolution (range widths) of the
  paired data. Note that small values for this argument may lead to R
  session crashes.

- width:

  Number of base pairs to expand the loci of interest in \`x\`.

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
pileupBoundaries(x=loops, files=hicFile, binSize=50e3)
#> Pairs are not binned to `binSize`.
#> ℹ Snapping to binSize=50000, 
#> ℹ Use `assignToBins()` or `snapToBins()` for more binning options.
#> Warning: ✖ `blockSize` must be twice the longest range.
#> ℹ Setting `blockSize` = 1000002.
#> <10 x 10> DelayedMatrix object of type "double":
#>             [,1]       [,2]       [,3] ...       [,9]      [,10]
#>  [1,] 19.1326531  6.9234694  3.1938776   .  0.8673469  0.8724490
#>  [2,]  6.9234694 19.0714286  7.0153061   .  1.0306122  0.8061224
#>  [3,]  3.1938776  7.0153061 18.0561224   .  1.2755102  0.9591837
#>  [4,]  2.0714286  2.9540816  6.2806122   .  1.1887755  1.2755102
#>  [5,]  1.9030612  2.2959184  2.7602041   .  1.6938776  1.4540816
#>  [6,]  1.3520408  1.6581633  2.2295918   .  2.1530612  1.7551020
#>  [7,]  1.2806122  1.2806122  1.6173469   .  3.1224490  2.1887755
#>  [8,]  1.0765306  1.2755102  1.2500000   .  6.6836735  3.0408163
#>  [9,]  0.8673469  1.0306122  1.2755102   . 18.4234694  6.6530612
#> [10,]  0.8724490  0.8061224  0.9591837   .  6.6530612 18.2704082
```
