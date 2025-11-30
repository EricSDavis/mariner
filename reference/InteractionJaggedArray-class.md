# InteractionJaggedArray Class

The \`InteractionJaggedArray\` class creates a container for storing
interaction data alongside irregular arrays. This allows the storage of
matrices with different dimensions on-disk using HDF5.

Subset an InteractionJaggedArray by its interactions (\[i,\]) or its
Hi-C files (\[,j\]).

## Usage

``` r
# S4 method for class 'InteractionJaggedArray'
show(object)

# S4 method for class 'InteractionJaggedArray'
dim(x)

# S4 method for class 'InteractionJaggedArray'
interactions(x)

# S4 method for class 'InteractionJaggedArray'
metadata(x)

# S4 method for class 'InteractionJaggedArray'
colData(x)

# S4 method for class 'InteractionJaggedArray'
counts(object)

# S4 method for class 'InteractionJaggedArray'
path(object)

# S4 method for class 'InteractionJaggedArray'
length(x)

# S4 method for class 'InteractionJaggedArray,ANY,ANY,ANY'
x[i, j]
```

## Arguments

- object:

  InteractionJaggedArray object.

- x:

  An InteractionJaggedArray object.

- i:

  Numeric vector indicating the indices of interactions to extract.

- j:

  Numeric vector indicating the indices of files to extract.

## Value

\`InteractionJaggedArray()\` creates an InteractionJaggedArray object.

\`dim()\` returns a list of the dimensions of the interactions, files,
and count matrices.

\`interactions()\` returns the interactions.

\`metadata()\` returns the metadata.

\`colData()\` returns the column data.

\`counts()\` returns the JaggedArray object containing count matrix
information.

\`path()\` returns a character vector with the path to the HDF5 file
with the JaggedArray data.

\`length()\` returns an integer with the number of interactions in an
InteractionJaggedArray object.

Subsetting returns an InteractionJaggedArray or InteractionArray object
(see Details).

## Details

The object returned will be a InteractionJaggedArray if the submatrices
contain different dimensions. However, the returned object will
automatically be coerced into a InteractionArray if possible (i.e. the
dimensions of the rows and columns of submatrices are the same.)

## Slots

- `interactions`:

  A GInteractions object.

- `colData`:

  Column data describing Hi-C files.

- `counts`:

  A JaggedArray object with data.

- `metadata`:

  List of metadata describing the object.

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
iarr
#> class: InteractionJaggedArray
#> dim: 4 interaction(s), 2 file(s), variable count matrix(es)
#> metadata(3): binSize, norm, matrix
#> colData: FS, WT
#> colData names(2): files, fileNames
#> HDF5: /tmp/Rtmp9FSwoA/file7b572694fe34.h5

## Show dimensions
dim(iarr)
#> $interactions
#> [1] 4
#> 
#> $files
#> [1] 2
#> 
#> $rows
#> [1] 3 3 5 3
#> 
#> $cols
#> [1] 5 5 3 8
#> 

## Access interactions
interactions(iarr)
#> GInteractions object with 4 interactions and 0 metadata columns:
#>       seqnames1             ranges1     seqnames2             ranges2
#>           <Rle>           <IRanges>         <Rle>           <IRanges>
#>   [1]         1   51000000-51300000 ---         1   51000000-51500000
#>   [2]         2   52000000-52300000 ---         3   52000000-52500000
#>   [3]         1 150000000-150500000 ---         1 150000000-150300000
#>   [4]         2   52000000-52300000 ---         2   52000000-52800000
#>   -------
#>   regions: 7 ranges and 0 metadata columns
#>   seqinfo: 3 sequences from an unspecified genome; no seqlengths

## Access metadata
metadata(iarr)
#> $binSize
#> [1] 1e+05
#> 
#> $norm
#> [1] "NONE"
#> 
#> $matrix
#> [1] "observed"
#> 

## Access colData
colData(iarr)
#> DataFrame with 2 rows and 2 columns
#>                     files         fileNames
#>               <character>       <character>
#> FS /home/runner/.cache/.. 7b5756449908_8147
#> WT /home/runner/.cache/.. 7b575c1774a3_8148

## Access count matrices
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

## Access path to HDF5 data
path(iarr)
#> [1] "/tmp/Rtmp9FSwoA/file7b572694fe34.h5"

## length
length(iarr)
#> [1] 4

## Subsetting
iarr[1:3,1]
#> class: InteractionJaggedArray
#> dim: 3 interaction(s), 1 file(s), variable count matrix(es)
#> metadata(3): binSize, norm, matrix
#> colData: FS
#> colData names(2): files, fileNames
#> HDF5: /tmp/Rtmp9FSwoA/file7b572694fe34.h5
```
