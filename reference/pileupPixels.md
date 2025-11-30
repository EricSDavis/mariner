# Pileup Hi-C pixels

pileupPixels optionally removes short interactions that intersect the
diagonal before extracting then aggregating square regions around each
pixel from Hi-C files. This is also known as aggregate peak analysis
(APA)

## Usage

``` r
pileupPixels(
  x,
  files,
  binSize,
  buffer = 5,
  removeShort = TRUE,
  minPairDist = 0,
  normalize = TRUE,
  FUN = sum,
  nBlocks = 5,
  verbose = TRUE,
  BPPARAM = bpparam(),
  ...
)

# S4 method for class 'GInteractions,character,numeric'
pileupPixels(
  x,
  files,
  binSize,
  buffer = 5,
  removeShort = TRUE,
  minPairDist = 0,
  normalize = TRUE,
  FUN = sum,
  nBlocks = 5,
  verbose = TRUE,
  BPPARAM = bpparam(),
  ...
)
```

## Arguments

- x:

  GInteractions object containing interactions to extract from Hi-C
  files. These should be pixels of a single \`binSize\` in width.

- files:

  Character file paths to \`.hic\` files.

- binSize:

  Integer (numeric) describing the resolution (range widths) of the
  paired data.

- buffer:

  Integer indicating the buffer size, or number of pixels

- removeShort:

  Boolean, whether to remove short pairs (Default) or not.

- minPairDist:

  Pairs with a distance less than or equal to this value will be
  filtered out.

- normalize:

  Boolean, whether to normalize the aggregated values to the number of
  interactions (after filtering out short pairs - if applicable).

- FUN:

  Function to use for aggregating.

- nBlocks:

  Number of blocks for block-processing arrays. Default is 5. Increase
  this for large datasets. To read and process all data at once, set
  this value to 1.

- verbose:

  Boolean (TRUE or FALSE) describing whether to report block-processing
  progress.

- BPPARAM:

  Parallelization params (passed to \`BiocParallel::bplapply()\`).
  Default is the result of \`BiocParallel::bpparams()\`. Parallel
  processing is not available when \`by=interactions\`.

- ...:

  Additional arguments passed to \`pullHicMatrices()\`.

## Value

A DelayedMatrix of aggregated counts.

## Details

Note that pair distance filtering is done after expanding interactions
to matrices.

## Examples

``` r
## Load marinerData
if (!require("marinerData", quietly = TRUE))
    BiocManager::install("marinerData")

## Read .hic file paths
hicFile <- marinerData::LEUK_HEK_PJA30_inter_30.hic()
#> see ?marinerData and browseVignettes('marinerData') for documentation
#> loading from cache
names(hicFile) <- "WT"

## Loops
loops <-
    WT_5kbLoops.txt() |>
    setNames("WT") |>
    read.table(header=TRUE, nrows=1000) |>
    as_ginteractions(keep.extra.columns=FALSE) |>
    assignToBins(binSize=5e3)
#> see ?marinerData and browseVignettes('marinerData') for documentation
#> loading from cache

## Removes the "chr" prefix for compatibility
## with the preprocessed hic files
GenomeInfoDb::seqlevelsStyle(loops) <- 'ENSEMBL'

## APA
mat <- pileupPixels(
    x=loops,
    files=hicFile,
    binSize=5e3,
    minPairDist=50e3,
    normalize=FALSE
)

```
