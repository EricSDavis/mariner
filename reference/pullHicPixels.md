# Pull contact frequency from \`.hic\` files

Pull contact frequency from \`.hic\` files

## Usage

``` r
pullHicPixels(
  x,
  files,
  binSize,
  ...,
  h5File = tempfile(fileext = ".h5"),
  half = "both",
  norm = "NONE",
  matrix = "observed",
  blockSize = 248956422,
  onDisk = TRUE,
  compressionLevel = 0,
  chunkSize = 1
)

# S4 method for class 'GInteractions,character,numeric'
pullHicPixels(
  x,
  files,
  binSize,
  h5File,
  half,
  norm,
  matrix,
  blockSize,
  onDisk,
  compressionLevel,
  chunkSize
)
```

## Arguments

- x:

  GInteractions object containing interactions to extract from Hi-C
  files.

- files:

  Character file paths to \`.hic\` files.

- binSize:

  Integer (numeric) describing the resolution (range widths) of the
  paired data.

- ...:

  Additional arguments.

- h5File:

  Character file path to save \`.h5\` file.

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

- norm:

  String (length one character vector) describing the Hi-C normalization
  to apply. Use \`strawr::readHicNormTypes()\` to see accepted values
  for each file in \`files\`.

- matrix:

  String (length one character vector) Type of matrix to extract. Must
  be one of "observed", "oe", or "expected". "observed" is observed
  counts, "oe" is observed/expected counts, "expected" is expected
  counts.

- blockSize:

  Number (length one numeric vector) describing the size in base-pairs
  to pull from each \`.hic\` file. Default is 248956422 (the length of
  the longest chromosome in the human hg38 genome). For large \`.hic\`
  files \`blockSize\` can be reduced to conserve the amount of data read
  in at a time. Larger \`blockSize\` values speed up performance, but
  use more memory.

- onDisk:

  Boolean (length one logical vector that is not NA) indicating whether
  extracted data should be stored on disk in an HDF5 file. Default is
  TRUE.

- compressionLevel:

  Number (length one numeric vector) between 0 (Default) and 9
  indicating the compression level used on HDF5 file.

- chunkSize:

  Number (length one numeric vector) indicating how many values of \`x\`
  to chunk for each write to HDF5 stored data. This has downstream
  implications for accessing subsets later. For small
  \`compressionLevel\` values use smaller \`chunkSize\` values and for
  large \`compressionLevel\` values use large (i.e. \`length(x)\`)
  values to improve performance.

## Value

InteractionSet object with a 2-dimensional array of Hi-C interactions
(rows) and Hi-C sample (columns).

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

## Read in loop pixels as GInteractions object
pixels <-
  WT_5kbLoops.txt() |>
  setNames("WT") |>
  read.table(header=TRUE) |>
  as_ginteractions(keep.extra.columns=FALSE) |>
  assignToBins(binSize=100e3)
#> see ?marinerData and browseVignettes('marinerData') for documentation
#> loading from cache

## Removes the "chr" prefix for compatibility
## with the preprocessed hic files
GenomeInfoDb::seqlevelsStyle(pixels) <- 'ENSEMBL'

## Extract the first 100 pixels
imat <- pullHicPixels(x=pixels[1:100],
                      files=hicFiles,
                      binSize=100e3)
imat
#> class: InteractionMatrix 
#> dim: count matrix with 100 interactions and 2 file(s)
#> metadata(3): binSize norm matrix
#> assays(1): counts
#> rownames: NULL
#> rowData names(0):
#> colnames(2): FS WT
#> colData names(2): files fileNames
#> type: GInteractions
#> regions: 12176

## Access count matrix
counts(imat)
#> <100 x 2> DelayedMatrix object of type "double":
#>        FS WT
#>   [1,]  6  4
#>   [2,]  6  6
#>   [3,] 33 30
#>   [4,]  3  2
#>   [5,]  6  6
#>    ...  .  .
#>  [96,]  3  1
#>  [97,]  2  1
#>  [98,]  2  1
#>  [99,] 10  7
#> [100,]  3  5
```
