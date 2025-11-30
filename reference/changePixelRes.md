# Change pixels from one resolution to another selecting the new pixel using Hi-C data.

A GInteractions object containing pixels of interest is resized to the
\`from\` resolution (if its not already), then count matrices are
extracted for each interaction and Hi-C file using the new \`to\`
resolution. Count matrices are aggregated by interactions with the
supplied \`aggFUN\`, and a new pixel is selected with the supplied
\`selectFUN\`. For large datasets, increase \`nBlocks\` to allow for
smaller blocks of data to be processed in memory.

## Usage

``` r
changePixelRes(
  x,
  files,
  from,
  to,
  aggFUN = sum,
  selectFUN = "which.max",
  nBlocks = 5,
  verbose = TRUE,
  norm = "KR",
  half = "upper",
  ...
)

# S4 method for class 'GInteractions,character'
changePixelRes(
  x,
  files,
  from,
  to,
  aggFUN = sum,
  selectFUN = "which.max",
  nBlocks = 5,
  verbose = TRUE,
  norm = "KR",
  half = "upper",
  ...
)
```

## Arguments

- x:

  GInteractions object.

- files:

  Character file paths to \`.hic\` files.

- from:

  Number (length one numeric vector) describing the resolution of \`x\`.
  Data will be binned to this value if it is not already binned.

- to:

  Number (length one numeric vector) describing the new resolution for
  the pixels.

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

- norm:

  String (length one character vector) describing the Hi-C normalization
  to apply. Use \`strawr::readHicNormTypes()\` to see accepted values
  for each file in \`files\`.

- half:

  String (character vector of length one) indicating whether to keep
  values for the upper triangular (\`half="upper"\`) where \`start1 \<
  start2\`, lower triangular (\`half="lower"\`) where \`start1 \>
  start2\`, or both (\`half="both"\`, default). When \`half="upper"\`
  all lower triangular values are \`NA\`. When \`half="lower"\` all
  upper triangular values are \`NA\`. When \`half="both"\` there are no
  \`NA\` values. For interchromosomal interactions there is no inherent
  directionality between chromosomes, so data is returned regardless of
  specified order.

- ...:

  Additional arguments passed to \`pullHicMatrices()\`. See
  ?\[\`pullHicMatrices\`\].

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

## Change pixel resolution from 2.5e6 to 500e3
changePixelRes(x=loops[1:5],
               files=hicFiles,
               from=2.5e6,
               to=500e3)
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
