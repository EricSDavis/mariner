# CountMatrix Class

A class for displaying dimnames associated with the count matrices
resulting from pullHicMatrices() \|\> counts(showDimnames=TRUE).

## Usage

``` r
# S4 method for class 'CountMatrix'
show(object)
```

## Arguments

- object:

  A CountMatrix object.

## Value

A CountMatrix object (clone of DelayedArray)

A DelayedArray object with dimnames for the first two dimensions.

## Details

This class is used only for attaching a "show" method.

## Slots

- `object`:

  InteractionArray object

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
#> 
#> see ?marinerData and browseVignettes('marinerData') for documentation
#> downloading 1 resources
#> retrieving 1 resource
#> 
#> loading from cache
#> see ?marinerData and browseVignettes('marinerData') for documentation
#> downloading 1 resources
#> retrieving 1 resource
#> 
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
#> downloading 1 resources
#> retrieving 1 resource
#> 
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
```
