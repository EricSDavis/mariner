# Accessor for h5File path from an InteractionMatrix

Returns the file path describing where the on-disk HDF5 data associated
with the InteractionMatrix object is stored.

This method circumvents the \`assays\<-\` and \`path\<-\` methods for
updating the HDF5 path because they are not accessible when the file
path is broken.

## Usage

``` r
# S4 method for class 'InteractionMatrix'
path(object)

# S4 method for class 'InteractionMatrix'
path(object) <- value
```

## Arguments

- object:

  InteractionMatrix object

- value:

  String (length-one character vector) to use for path replacement.

## Value

The path to the HDF5 file associated with the InteractionMatrix object.

Updates path to HDF5 file for the InteractionMatrix object.

## Details

If the file no longer exists, the path is returned along with a warning.

This allows the file path to be updated even if the original linked data
no longer exists.

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

#################################
## Accessing path to HDF5 data ##
#################################

## Create example interactions
x <- read.table(text="
        9 14000000 14500000 9 14500000 15000000
        9 89500000 90000000 9 89500000 90000000
        9 23500000 24000000 9 23500000 24000000")
x <- as_ginteractions(x)

## Extract 3 pixels from 2 hic files
imat <- pullHicPixels(x, hicFiles, 500e03)

## Access path
path(imat)
#> [1] "/tmp/Rtmp9FSwoA/file7b57719f0aa1.h5"

#################################
## Updating path to HDF5 data ##
################################

## Create example interactions
x <- read.table(text="
        9 14000000 14500000 9 14500000 15000000
        9 89500000 90000000 9 89500000 90000000
        9 23500000 24000000 9 23500000 24000000")
x <- as_ginteractions(x)

## Extract 3 pixels from 2 hic files
h5File <- tempfile(fileext=".h5")
imat <- pullHicPixels(x, hicFiles, 500e03, h5File=h5File)

## Move file to new location
newFile <- tempfile(fileext="_new.h5")
file.rename(from=h5File, to=newFile)
#> [1] TRUE

## Update path
path(imat) <- newFile
path(imat)
#> [1] "/tmp/Rtmp9FSwoA/file7b5766c480f6_new.h5"
```
