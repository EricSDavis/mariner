# HDF5-backed blockApply

Read in array data in blocks, apply function, and write back to an HDF5
file.

## Usage

``` r
hdf5BlockApply(x, FUN, sink, grid, sink_grid, verbose = TRUE)

# S4 method for class 'DelayedArray'
hdf5BlockApply(x, FUN, sink, grid, sink_grid, verbose = TRUE)
```

## Arguments

- x:

  Delayed Array object.

- FUN:

  Function that takes one argument 'block' and processes it.

- sink:

  HDF5RealizationSink object.

- grid:

  ArrayGrid over array \`x\`.

- sink_grid:

  ArrayGrid over \`sink\`.

- verbose:

  Logical - whether block processing progress should be displayed.

## Value

An HDF5Array object.

## Details

Implements an HDF5-backed option for block processing on DelayedArray
objects.

## Examples

``` r
## ################################################
## This function is intended for advanced users.
## To learn more about using DelayedArray
## or HDF5-backed objects, see ?DelayedArray or
## ?HDF5Array
###################################################

library(DelayedArray)
#> Loading required package: Matrix
#> 
#> Attaching package: ‘Matrix’
#> The following object is masked from ‘package:S4Vectors’:
#> 
#>     expand
#> Loading required package: S4Arrays
#> Loading required package: abind
#> 
#> Attaching package: ‘S4Arrays’
#> The following object is masked from ‘package:abind’:
#> 
#>     abind
#> The following object is masked from ‘package:base’:
#> 
#>     rowsum
#> Loading required package: SparseArray
#> 
#> Attaching package: ‘DelayedArray’
#> The following objects are masked from ‘package:base’:
#> 
#>     apply, scale, sweep
library(HDF5Array)
#> Loading required package: h5mread
#> Loading required package: rhdf5
#> 
#> Attaching package: ‘h5mread’
#> The following object is masked from ‘package:rhdf5’:
#> 
#>     h5ls
library(rhdf5)

## Create example array that is longer in the
## 3rd dimension (representing interactions)
dims <- c(11L, 11L, 100L, 2L)
a <- array(data=seq(1, prod(dims)), dim=dims)
a <- DelayedArray(a)

## Define spacings, breaking up the longest dim
## Here we are processing in blocks of 10
spacings <- dim(a)
spacings[3] <- ceiling(spacings[3]/10)

## Define storage dimensions (all except those
## over which the function is being applied)
storageDims <- dims[c(1,2,3)]

## Define chunk dimensions for writing to HDF5
chunkDims <- storageDims
chunkDims[3] <- spacings[3]

## Create grid for applying the data (grid)
## and grid for writing to the sink (sink_grid)
grid <- RegularArrayGrid(dims, spacings)
sink_grid <- RegularArrayGrid(storageDims, chunkDims)

## Create HDF5 file for writing
h5 <- tempfile(fileext = ".h5")
h5createFile(h5)

## Define compression for HDF5
compressionLevel <- 0

## Create HDF5-backed realization sink
sink <- HDF5RealizationSink(filepath=h5,
                            name="counts",
                            type="integer",
                            dim=storageDims,
                            chunkdim=chunkDims,
                            level=compressionLevel)

## Wrap function that operates on each block
## this can be anything, here it is sum
FUN <- \(block) apply(block, c(1,2,3), sum)

## Read, apply, and write to HDF5
ans <- hdf5BlockApply(x=a,
                      FUN=FUN,
                      sink=sink,
                      grid=grid,
                      sink_grid=sink_grid,
                      verbose=TRUE)
#> / Reading and realizing block 1/10 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 2/10 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 3/10 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 4/10 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 5/10 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 6/10 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 7/10 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 8/10 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 9/10 ... 
#> OK
#> \ Processing it ... 
#> OK
#> / Reading and realizing block 10/10 ... 
#> OK
#> \ Processing it ... 
#> OK
ans
#> <11 x 11 x 100> HDF5Array object of type "integer":
#> ,,1
#>        [,1]  [,2]  [,3]  [,4] ...  [,8]  [,9] [,10] [,11]
#>  [1,] 12102 12124 12146 12168   . 12256 12278 12300 12322
#>  [2,] 12104 12126 12148 12170   . 12258 12280 12302 12324
#>   ...     .     .     .     .   .     .     .     .     .
#> [10,] 12120 12142 12164 12186   . 12274 12296 12318 12340
#> [11,] 12122 12144 12166 12188   . 12276 12298 12320 12342
#> 
#> ...
#> 
#> ,,100
#>        [,1]  [,2]  [,3]  [,4] ...  [,8]  [,9] [,10] [,11]
#>  [1,] 36060 36082 36104 36126   . 36214 36236 36258 36280
#>  [2,] 36062 36084 36106 36128   . 36216 36238 36260 36282
#>   ...     .     .     .     .   .     .     .     .     .
#> [10,] 36078 36100 36122 36144   . 36232 36254 36276 36298
#> [11,] 36080 36102 36124 36146   . 36234 36256 36278 36300
#> 
```
