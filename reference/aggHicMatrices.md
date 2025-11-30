# Aggregate count matrices from InteractionArray objects

Aggregation of count matrices is done blocks to avoid large memory
usage. Use \`nBlocks\` to control the number of blocks read into memory
at once. Blocks are defined as \`length(interactions(x))/nBlocks\`.

## Usage

``` r
aggHicMatrices(
  x,
  by = NULL,
  FUN = sum,
  nBlocks = 5,
  verbose = TRUE,
  BPPARAM = bpparam(),
  compressionLevel = 0
)

# S4 method for class 'InteractionArray'
aggHicMatrices(
  x,
  by = NULL,
  FUN = sum,
  nBlocks = 5,
  verbose = TRUE,
  BPPARAM = bpparam(),
  compressionLevel = 0
)
```

## Arguments

- x:

  InteractionArray object.

- by:

  String (length one character vector) describing whether to aggregate
  by interactions, files, or neither (i.e. NULL as default).

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

- compressionLevel:

  Number (length one numeric vector) between 0 (Default) and 9
  indicating the compression level used on HDF5 file.

## Value

An aggregated \`DelayedArray\` object. If \`by=interactions\` or
\`by=files\` then a 3-dimensional \`DelayedArray\` is returned. If
\`by=NULL\` (default) then A 2-dimensional \`DelayedMatrix\` is
returned.

## Details

Since interactions are typically the largest dimension in an
InteractionArray, using \`by=interactions\` creates an HDF5-backed array
to store these large arrays. Currently parallel processing for
HDF5-backed arrays are not supported regardless of the value of
\`BPPARAM\`.

Both \`by=NULL\` and \`by=files\` support parallel processing.

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

## Expand pixel ranges with a 5 pixel buffer on either side
loops <-
    assignToBins(loops, binSize=100e3) |>
    pixelsToMatrices(buffer=5)

## Extract 10, 11x11 count matrices from 2 hic files
iarr <-
    loops[1:10] |>
    pullHicMatrices(binSize=100e3,
                    files=hicFiles)

## Aggregate all, by files, or by interactions
aggHicMatrices(x=iarr)
#> <11 x 11> DelayedMatrix object of type "double":
#>        [,1]  [,2]  [,3] ... [,10] [,11]
#>  [1,]   250    85    44   .     9    27
#>  [2,]   203   201    87   .    14    23
#>  [3,]   231   210   255   .    18    24
#>  [4,]   320   225   201   .    25    28
#>  [5,]   123   279   177   .    22    23
#>  [6,]   254   170   294   .    43    38
#>  [7,]   116   216   169   .    40    30
#>  [8,]    93   104   276   .    55    42
#>  [9,]    56   118   110   .    89    61
#> [10,]    42    53   114   .   175    85
#> [11,]    34    46    46   .   200   190
aggHicMatrices(x=iarr, by="files")
#> <11 x 11 x 2> DelayedArray object of type "double":
#> ,,FS
#>        [,1]  [,2]  [,3] ... [,10] [,11]
#>  [1,]   143    51    27   .     6    13
#>  [2,]   108   110    44   .     8    15
#>   ...     .     .     .   .     .     .
#> [10,]    25    26    63   .    71    54
#> [11,]    11    24    18   .    95    97
#> 
#> ,,WT
#>        [,1]  [,2]  [,3] ... [,10] [,11]
#>  [1,]   107    34    17   .     3    14
#>  [2,]    95    91    43   .     6     8
#>   ...     .     .     .   .     .     .
#> [10,]    17    27    51   .   104    31
#> [11,]    23    22    28   .   105    93
#> 
aggHicMatrices(x=iarr, by="interactions")
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
#> <11 x 11 x 10> DelayedArray object of type "double":
#> ,,1
#>        [,1]  [,2]  [,3] ... [,10] [,11]
#>  [1,]     7     4     6   .     0     0
#>  [2,]     8    11     5   .     3     1
#>   ...     .     .     .   .     .     .
#> [10,]     2     2     8   .     8     2
#> [11,]     2     1     3   .     9     6
#> 
#> ...
#> 
#> ,,10
#>        [,1]  [,2]  [,3] ... [,10] [,11]
#>  [1,]     5     3     0   .     0     3
#>  [2,]     3     5     4   .     0     2
#>   ...     .     .     .   .     .     .
#> [10,]     8    11    53   .     4     4
#> [11,]     8     7    10   .     6     8
#> 
```
