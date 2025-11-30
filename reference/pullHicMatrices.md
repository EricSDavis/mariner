# Pull submatrices from \`.hic\` files

The dimensions of the pulled submatrix is defined by dividing the widths
of anchors in \`x\` by the \`binSize\`. When the anchor widths are the
same for each interaction, an InteractionArray is returned. However, if
the anchor widths differ in \`x\`, an InteractionJaggedArray is returned
instead.

## Usage

``` r
pullHicMatrices(
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
pullHicMatrices(
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

InteractionSet object with a 4-dimensional array of Hi-C submatrices,
rownames, and colnames. Array is stored with the following dimensions:
Interactions in \`x\`, Hi-C \`files\`, rows of submatrix, columns of
submatrix. The submatrices returned have rows cooresponding to anchor1
of \`x\` and columns correspond to anchor2 of \`x\`.

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

## Expand pixels to regions for pulling
## Hi-C submatrices
regions <- pixelsToMatrices(x=pixels, buffer=5)

## Extract 11x11 count matrices from the
## first 100 regions and 2 Hi-C files
iarr <- pullHicMatrices(x=regions[1:100],
                        files=hicFiles,
                        binSize=100e3)
iarr
#> class: InteractionArray 
#> dim: 100 interaction(s), 2 file(s), 11x11 count matrix(es)
#> metadata(3): binSize norm matrix
#> assays(3): counts rownames colnames
#> rownames: NULL
#> rowData names(0):
#> colnames(2): FS WT
#> colData names(2): files fileNames
#> type: GInteractions
#> regions: 12176

## Access count matrices
counts(iarr)
#> <11 x 11 x 100 x 2> DelayedArray object of type "double":
#> ,,1,FS
#>        [,1]  [,2]  [,3] ... [,10] [,11]
#>  [1,]     4     3     4   .     0     0
#>  [2,]     6     9     3   .     2     1
#>   ...     .     .     .   .     .     .
#> [10,]     1     1     4   .     2     2
#> [11,]     1     1     0   .     5     2
#> 
#> ...
#> 
#> ,,100,WT
#>        [,1]  [,2]  [,3] ... [,10] [,11]
#>  [1,]     0     0     1   .     0     0
#>  [2,]     2     0     0   .     1     1
#>   ...     .     .     .   .     .     .
#> [10,]     5     8    32   .     0     5
#> [11,]     2     4     8   .     0     1
#> 

## Display the start bin of each
## interaction in the count
## matrices
counts(iarr, showDimnames=TRUE)
#> <11 x 11 x 100 x 2> DelayedArray object of type "double":
#> ,,1,FS
#>          14200000 14300000 14400000 ... 15100000 15200000
#> 13900000        4        3        4   .        0        0
#> 14000000        6        9        3   .        2        1
#> 14100000       11        8        4   .        1        1
#> 14200000       38        8        4   .        0        0
#> 14300000        8       31       11   .        0        1
#> 14400000        4       11       35   .        3        2
#> 14500000        4        1       12   .        1        0
#> 14600000        2        0        4   .        3        2
#> 14700000        1        2        6   .        3        3
#> 14800000        1        1        4   .        2        2
#> 14900000        1        1        0   .        5        2
#> 
#> ...
#> 
#> ,,100,WT
#>          16300000 16400000 16500000 ... 17200000 17300000
#> 15600000        0        0        1   .        0        0
#> 15700000        2        0        0   .        1        1
#> 15800000        1        2        1   .        0        1
#> 15900000        2        3        2   .        1        1
#> 16000000        1        2        0   .        1        1
#> 16100000        6        3        2   .        0        0
#> 16200000        4        4        3   .        1        0
#> 16300000       25        8        5   .        0        1
#> 16400000        8       19        8   .        0        0
#> 16500000        5        8       32   .        0        5
#> 16600000        2        4        8   .        0        1

## InteractionJaggedArray example
gi <- read.table(text="
            1 51000000 51300000 1 51000000 51500000
            2 52000000 52300000 3 52000000 52500000
            1 150000000 150500000 1 150000000 150300000
            2 52000000 52300000 2 52000000 52800000") |>
    as_ginteractions()

iarr <- pullHicMatrices(gi, hicFiles, 100e03, half="both")
iarr
#> class: InteractionJaggedArray
#> dim: 4 interaction(s), 2 file(s), variable count matrix(es)
#> metadata(3): binSize, norm, matrix
#> colData: FS, WT
#> colData names(2): files, fileNames
#> HDF5: /tmp/Rtmp9FSwoA/file7b57418ce8d6.h5

counts(iarr)
#> <n x m x 4 x 2> JaggedArray:
#> ,,1,1
#> <3 x 5> matrix
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]   53   15    5    1    4
#> [2,]   15   68   19    8    5
#> [3,]    5   19   69   12    2
#> 
#> ...
#> 
#> ,,4,2
#> <3 x 8> matrix
#>      [,1] [,2] [,3] ... [,7] [,8]
#> [1,]   31    7    2   .    0    0
#> [2,]    7   22    5   .    1    1
#> [3,]    2    5   26   .    0    2
#> 
```
