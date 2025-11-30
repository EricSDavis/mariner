# Regularize JaggedArray or InteractionJaggedArray objects

InteractionJaggedArray objects and their count matrices (JaggedArray
objects) contain variable dimension matrices. The \`regularize\`
function resizes these matrices to the new dimensions supplied in
\`ndim\`. The result is a DelayedArray object (for JaggedArray) or an
InteractionArray object (for InteractionJaggedArray).

## Usage

``` r
regularize(
  x,
  ndim = c(10, 10),
  h5File = tempfile(fileext = ".h5"),
  scale = TRUE,
  nBlocks = 5,
  verbose = TRUE,
  chunkSize = 1,
  compressionLevel = 0,
  ...
)

# S4 method for class 'JaggedArray'
regularize(
  x,
  ndim,
  h5File,
  scale,
  nBlocks,
  verbose,
  chunkSize,
  compressionLevel
)

# S4 method for class 'InteractionJaggedArray'
regularize(
  x,
  ndim,
  h5File,
  scale,
  nBlocks,
  verbose,
  chunkSize,
  compressionLevel
)
```

## Arguments

- x:

  A JaggedArray or InteractionJaggedArray object.

- ndim:

  Numeric vector of length two describing the new dimensions of the
  output matrices.

- h5File:

  Character file path to save \`.h5\` file.

- scale:

  Boolean (TRUE/FALSE) indicating whether the values in the new matrices
  should be scaled to the total signal in each matrix.

- nBlocks:

  Number of blocks for block-processing JaggedArrays. Default is 5.
  Increase this for large datasets. To read and process all data at
  once, set this value to 1.

- verbose:

  Boolean (TRUE or FALSE) describing whether to report block-processing
  progress.

- chunkSize:

  Number (length one numeric vector) indicating how many values of \`x\`
  to chunk for each write to HDF5 stored data. This has downstream
  implications for accessing subsets later. For small
  \`compressionLevel\` values use smaller \`chunkSize\` values and for
  large \`compressionLevel\` values use large (i.e. \`length(x)\`)
  values to improve performance.

- compressionLevel:

  Number (length one numeric vector) between 0 (Default) and 9
  indicating the compression level used on HDF5 file.

- ...:

  Additional arguments.

## Value

If \`x\` is a JaggedArray then \`regularize\` returns an HDF5-backed
4-dimensional DelayedArray object where the first and second dimensions
are the rows and columns of the count matrices (\`ndim\`), the third
dimension is the number of interactions and the fourth dimension is the
number of files. If \`x\` is an InteractionJaggedArray then an
InteractionArray object is returned where counts returns the object
described above.

## Details

Note that the interaction/binSize/count matrices relationship will be
inconsistent in the resulting InteractionArray object and the row/col
names will not be available.

## Examples

``` r
## Load marinerData
if (!require("marinerData", quietly = TRUE))
    BiocManager::install("marinerData")

## Read .hic file paths
hicFiles <- c(
    LEUK_HEK_PJA27_inter_30.hic(),
    LEUK_HEK_PJA30_inter_30.hic()
)
#> see ?marinerData and browseVignettes('marinerData') for documentation
#> loading from cache
#> see ?marinerData and browseVignettes('marinerData') for documentation
#> loading from cache
names(hicFiles) <- c("FS", "WT")

## Create test interactions
gi <- read.table(text="
            1 51000000 51300000 1 51000000 51500000
            2 52000000 52300000 3 52000000 52500000
            1 150000000 150500000 1 150000000 150300000
            2 52000000 52300000 2 52000000 52800000") |>
    as_ginteractions()
gi <- c(gi,gi) # make more interactions

## InteractionJaggedArray object
ija <- pullHicMatrices(gi, hicFiles, 100e03, half="both")

## Regularize InteractionJaggedArray
ia <- regularize(ija, ndim=c(5,5), nBlocks=1)
#> / Reading and realizing block 1/1 of file 1/2 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 1/1 of file 2/2 ... 
#> OK
#> \ Processing it ... 
#> OK
aggHicMatrices(ia, nBlocks=1)
#> / reading and realizing block 1/1 ... 
#> ok
#> \ processing it ... 
#> ok
#> <5 x 5> DelayedMatrix object of type "double":
#>           [,1]      [,2]      [,3]      [,4]      [,5]
#> [1,] 1.5333093 0.5283683 0.2446879 0.2068139 0.1215322
#> [2,] 0.9166665 0.8038605 0.5434594 0.3513887 0.1819409
#> [3,] 0.4196802 0.8874755 0.3388203 0.4250340 0.6039009
#> [4,] 0.2658102 0.7882093 1.1711284 0.2252012 0.2712338
#> [5,] 0.1107120 0.7057650 2.0383085 0.1545785 0.1621146

## Regularize JaggedArray
ja <- counts(ija)
regularize(ja, ndim=c(5,5), nBlocks=1)
#> / Reading and realizing block 1/1 of file 1/2 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 1/1 of file 2/2 ... 
#> OK
#> \ Processing it ... 
#> OK
#> <5 x 5 x 8 x 2> HDF5Array object of type "double":
#> ,,1,1
#>             [,1]        [,2]        [,3]        [,4]        [,5]
#> [1,] 0.104433498 0.029556650 0.009852217 0.001970443 0.007881773
#> [2,] 0.066995074 0.081773399 0.023645320 0.008866995 0.008866995
#> [3,] 0.029556650 0.133990148 0.037438424 0.015763547 0.009852217
#> [4,] 0.019704433 0.085714286 0.086699507 0.019704433 0.006896552
#> [5,] 0.009852217 0.037438424 0.135960591 0.023645320 0.003940887
#> 
#> ...
#> 
#> ,,8,2
#>             [,1]        [,2]        [,3]        [,4]        [,5]
#> [1,] 0.203112203 0.021294021 0.013104013 0.014742015 0.000000000
#> [2,] 0.124488124 0.040950041 0.018018018 0.013104013 0.003276003
#> [3,] 0.045864046 0.060606061 0.022932023 0.011466011 0.006552007
#> [4,] 0.029484029 0.098280098 0.034398034 0.010647011 0.009828010
#> [5,] 0.013104013 0.135954136 0.045864046 0.009828010 0.013104013
#> 
```
