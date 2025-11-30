# JaggedArray Class

The \`JaggedArray\` class creates a container for storing irregular or
jagged array data. This allows the storage of matrices with different
dimensions on-disk using HDF5.

Subset a JaggedArray by its interactions (\[i,\]) or its Hi-C files
(\[,j\]).

\`as.list\` reads the on-disk data and returns it as an in-memory list
of matrices.

## Usage

``` r
# S4 method for class 'JaggedArray'
show(object)

# S4 method for class 'JaggedArray,ANY,ANY,ANY'
x[i, j, ..., drop = TRUE]

# S4 method for class 'JaggedArray'
as.list(x)

# S4 method for class 'JaggedArray'
path(object)

# S4 method for class 'JaggedArray'
dim(x)
```

## Arguments

- object:

  JaggedArray object.

- x:

  JaggedArray object.

- i:

  Numeric vector indicating the indices of interactions to extract.

- j:

  Numeric vector indicating the indices of files to extract.

- ...:

  Additional indices for subsetting multidimensional arrays.

- drop:

  Not accepted for JaggedArray objects.

## Value

\`JaggedArray()\` creates a JaggedArray object.

Subsetting returns a JaggedArray or DelayedArray object (see Details).

\`as.list()\` returns a list of matrices.

\`path()\` returns a character vector with the path to the HDF5 file
with the JaggedArray data.

\`dim()\` returns a list of dimensions of the JaggedArray of rows, cols,
interactions and files.

## Details

NOTE: This class is designed specifically for holding a 4-dimensional
JaggedArray \<n x m x i x j\> where n x m are rows and cols of count
matrices, i is interactions, and j is Hi-C files.

The object returned will be a JaggedArray if the submatrices contain
different dimensions. However, the returned object will automatically be
coerced into a DelayedArray if possible (i.e. the dimensions of the rows
and columns are the same.)

The JaggedArray data is still stored on-disk in an HDF5 file until it is
coerced into a DelayedArray or realized as a list of matrices.

## Slots

- `h5File`:

  path to file for creating and storing data as an HDF5 file.

- `dim`:

  dimensions describing the number of matrices contained. dim\[1\] is
  the number of interactions, dim\[2\] is the number of files.

- `subList`:

  is a list of length 2 where the first position refers to interactions
  and the second refers to files. This list is used to record subsetting
  operations which are then later applied when accessing data stored in
  the HDF5 file.

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

## Create test interactions
gi <- read.table(text="
            1 51000000 51300000 1 51000000 51500000
            2 52000000 52300000 3 52000000 52500000
            1 150000000 150500000 1 150000000 150300000
            2 52000000 52300000 2 52000000 52800000") |>
    as_ginteractions()

## InteractionJaggedArray object
iarr <- pullHicMatrices(gi, hicFiles, 100e03, half="both")
arr <- counts(iarr)
arr
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

## Subsetting
arr[,,1,] # DelayedArray
#> <3 x 5 x 1 x 2> DelayedArray object of type "double":
#> ,,1,1
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]   53   15    5    1    4
#> [2,]   15   68   19    8    5
#> [3,]    5   19   69   12    2
#> 
#> ,,1,2
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]   49   27    4    5    3
#> [2,]   27   49   13    2    6
#> [3,]    4   13   56    7    8
#> 
arr[,,,1] # JaggedArray
#> <n x m x 4 x 1> JaggedArray:
#> ,,1,1
#> <3 x 5> matrix
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]   53   15    5    1    4
#> [2,]   15   68   19    8    5
#> [3,]    5   19   69   12    2
#> 
#> ...
#> 
#> ,,4,1
#> <3 x 8> matrix
#>      [,1] [,2] [,3] ... [,7] [,8]
#> [1,]   29    7    3   .    1    2
#> [2,]    7   28   12   .    0    0
#> [3,]    3   12   29   .    0    1
#> 

## Realize as list
as.list(arr)
#> [[1]]
#> [[1]][[1]]
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]   53   15    5    1    4
#> [2,]   15   68   19    8    5
#> [3,]    5   19   69   12    2
#> 
#> [[1]][[2]]
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    0    0    0    0    0
#> [2,]    0    0    0    0    0
#> [3,]    0    0    0    0    0
#> 
#> [[1]][[3]]
#>      [,1] [,2] [,3]
#> [1,]   63   25   15
#> [2,]   25   68   28
#> [3,]   15   28   87
#> [4,]    4    7   45
#> [5,]    2    3   11
#> 
#> [[1]][[4]]
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
#> [1,]   29    7    3    3    1    4    1    2
#> [2,]    7   28   12    4    2    3    0    0
#> [3,]    3   12   29   11    3    1    0    1
#> 
#> 
#> [[2]]
#> [[2]][[1]]
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]   49   27    4    5    3
#> [2,]   27   49   13    2    6
#> [3,]    4   13   56    7    8
#> 
#> [[2]][[2]]
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    0    0    0    0    0
#> [2,]    0    0    0    0    0
#> [3,]    0    0    1    0    0
#> 
#> [[2]][[3]]
#>      [,1] [,2] [,3]
#> [1,]   56   26    8
#> [2,]   26   60   14
#> [3,]    8   14   89
#> [4,]    9   10   22
#> [5,]    2    2   15
#> 
#> [[2]][[4]]
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
#> [1,]   31    7    2    2    2    3    0    0
#> [2,]    7   22    5    6    1    2    1    1
#> [3,]    2    5   26   10    4    2    0    2
#> 
#> 

## Find the data path
path(arr)
#> [1] "/tmp/Rtmp9FSwoA/file7b57b90f1b7.h5"

## Find the data path
dim(arr)
#> $rows
#> [1] 3 3 5 3
#> 
#> $cols
#> [1] 5 5 3 8
#> 
#> $interactions
#> [1] 4
#> 
#> $files
#> [1] 2
#> 
```
