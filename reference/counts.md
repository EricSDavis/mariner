# Access count matrices from InteractionArray or InteractionMatrix

Access count matrices from InteractionArray or InteractionMatrix

Access count matrices from InteractionArray or InteractionMatrix

Replace method for counts

## Usage

``` r
# S4 method for class 'InteractionArray'
counts(object, showDimnames = FALSE)

# S4 method for class 'InteractionMatrix'
counts(object)

# S4 method for class 'InteractionMatrix'
counts(object) <- value
```

## Arguments

- object:

  InteractionMatrix object

- showDimnames:

  Logical vector of length-one indicating whether to show dimensions of
  count matrices (default FALSE). Only applicable for InteractionArray
  objects.

- value:

  Value for replacement

## Value

For InteractionArray, a 4-dimensional DelayedArray of Hi-C submatrices
is returned with the following dimensions: rows of count matrix, columns
of count matrix, Interactions in \`object\`, Hi-C \`files\`.

For InteractionMatrix, a 2-dimensional DelayedArray is returned with
rows representing interactions in \`object\` and columns for each Hi-C
file in \`files\`.

For InteractionMatrix, the replace matrix replaces the counts assay with
matrix-like objects supplied in \`value\`.

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

######################################
## Accessing Hi-C count submatrices ##
######################################

## Create example interactions
x <- read.table(text="
        9 14435000 14490000 9 14740000 14795000
        9 89540000 89595000 9 89785000 89840000
        9 23700000 23755000 9 23760000 23815000")
x <- as_ginteractions(x)

## Extract 3, 11x11 count matrices from 2 hic files
iarr <- pullHicMatrices(x, hicFiles, 5e03)

## Access count matrices
counts(iarr)
#> <11 x 11 x 3 x 2> DelayedArray object of type "double":
#> ,,1,FS
#>        [,1]  [,2]  [,3] ... [,10] [,11]
#>  [1,]     0     0     0   .     0     0
#>  [2,]     0     0     0   .     0     0
#>   ...     .     .     .   .     .     .
#> [10,]     0     0     0   .     0     0
#> [11,]     0     0     0   .     0     0
#> 
#> ...
#> 
#> ,,3,WT
#>        [,1]  [,2]  [,3] ... [,10] [,11]
#>  [1,]     0     0     0   .     0     0
#>  [2,]     0     0     0   .     0     0
#>   ...     .     .     .   .     .     .
#> [10,]     0     0     0   .     0     0
#> [11,]     0     0     0   .     0     0
#> 
counts(iarr, FALSE)
#> <11 x 11 x 3 x 2> DelayedArray object of type "double":
#> ,,1,FS
#>        [,1]  [,2]  [,3] ... [,10] [,11]
#>  [1,]     0     0     0   .     0     0
#>  [2,]     0     0     0   .     0     0
#>   ...     .     .     .   .     .     .
#> [10,]     0     0     0   .     0     0
#> [11,]     0     0     0   .     0     0
#> 
#> ...
#> 
#> ,,3,WT
#>        [,1]  [,2]  [,3] ... [,10] [,11]
#>  [1,]     0     0     0   .     0     0
#>  [2,]     0     0     0   .     0     0
#>   ...     .     .     .   .     .     .
#> [10,]     0     0     0   .     0     0
#> [11,]     0     0     0   .     0     0
#> 

#################################
## Accessing Hi-C count matrix ##
#################################

## Create example interactions
x <- read.table(text="
        9 14000000 14500000 9 14500000 15000000
        9 89500000 90000000 9 89500000 90000000
        9 23500000 24000000 9 23500000 24000000")
x <- as_ginteractions(x)

## Extract 3 pixels from 2 hic files
imat <- pullHicPixels(x, hicFiles, 500e03)

## Access count matrix
counts(imat)
#> <3 x 2> DelayedMatrix object of type "double":
#>       FS  WT
#> [1,]  49  38
#> [2,] 314 300
#> [3,] 258 187

#################################
## Replacing Hi-C count matrix ##
#################################

## Realize as in-memory matrix
counts(imat) <- as.matrix(counts(imat))
counts(imat)
#>       FS  WT
#> [1,]  49  38
#> [2,] 314 300
#> [3,] 258 187
imat
#> class: InteractionMatrix 
#> dim: count matrix with 3 interactions and 2 file(s)
#> metadata(3): binSize norm matrix
#> assays(1): counts
#> rownames: NULL
#> rowData names(0):
#> colnames(2): FS WT
#> colData names(2): files fileNames
#> type: GInteractions
#> regions: 4
```
