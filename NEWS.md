# mariner 1.2.1

New features:

Added `pileupBoundaries()` for visualizing pileup
Hi-C contacts in regions around genomic loci/boundary
elements.

# mariner 1.2.0

Breaking changes:

* Function renaming:
    * `subsetBySource()` -> `sets()`
    * `getPairClusters()` -> `clusters()`
    * `aggPairMcols()` -> `aggMetadata()`
    * `binPairs()` -> `assignToBins()`

* `sets()` function now returns all combinations
of sets from a `MergedGInteractions` object.

# mariner 1.1.3

New features:

Wrapper functions for performing pileup
analysis on pixels or domains.

* `pileupPixels()` for extracting and aggregating
Hi-C counts in square regions around pixels of
interest (aggregate peak analysis).

* `pileupDomains()` for extracting, resizing, and
aggregating Hi-C counts in domain regions of
interest (aggregate domain analysis).

Changes:

* Improve documentation for `removeShortPairs()`
to make it more clear that it should be run after
any function that resizes interactions.

* Add `verbose` option to `regularize()` to update
users on regularization of jagged arrays.

# mariner 1.1.2

Bug fixes:

* Fix id mapping in `aggPairMcols()`

# mariner 1.1.1

Bug fixes:

* Fix in `pullHic` functions:
After running `.prepareInputs()` which involves snapping ranges
to bins, the adjusted ranges were not being used to calculate the
expected matrix dimensions. This can sometimes cause a missmatch
between the data used for enumerating bins, and the data that is
extracted from the Hi-C file. Fixed by using the adjusted
interactions to set the matrix dimensions.

# mariner 1.1.0

New features:

* `JaggedArray` and `InteractionJaggedArray` classes
for irregular matrices.

* Functions for generating random `GRanges` and
`GInteractions` objects.

* `regularize` method for converting irregular (i.e. jagged)
to regular arrays.

* New `calcLoopEnrichment` method  for `InteractionArray`
objects.

* `defaultBuffer()` function for setting the buffer argument
from an `InteractionArray`.

Bug Fixes:

* Fix bug in `removeShortPairs` where padding wasn't
working as intended.

* `FUN` argument of `calcLoopEnrichment` now accepts
environmental variables and uses flexible argument names
for `fg` and `bg`.

# mariner 0.99.0

Bug fixes and improvements:

* Improve dispatch speed of `mergePairs()` by removing
S4 method dispatch on all arguments to just `x` and
`radius`.

* Fix bug in `mergePairs()` where all pairs are altered
during mean of mode transformation. Now original pairs
are preserved when accessed with `clusters()`.

* Set replace method for `counts<-` accessor for
`InteractionMatrix` objects. Helpful for converting
`DelayedMatrix` to `matrix`.

* Update `pixelsToMatrix` to preserve metadata columns
and include some additional tests.

* Add `plotMatrix()` function for plotting matrix data
as a heatmap. Useful for visualizing `DelayedMatrices`
from `pullHicMatrices()` and `aggHicMatrices()`.
Compatible with `plotgardener` package.

    * Allow `plotMatrix()` to accept `na.color`

* Bug fix in `mergePairs()` that allows columns
named "radius" and/or "method".

* Swap "binSize" and "files" argument order
in `pullHicPixels` and `pullHicMatrices`

* Allow `pullHicPixels` to overwrite existing
HDF5 files.

* Validity checks and functions to access/update
the HDF5 paths for `InteractionMatrix` objects,
even when those paths have been broken.

* Add temporary `plotBullseye` function.

* Selection functions for selecting indices of a
matrix:
    
    * `selectCenterPixel`
    * `selectRadius`
    * `selectSubmatrix`
    * `selectCoordinates`
    * `selectBlock`
    * `selectTopLeft`
    * `selectTopRight`
    * `selectBottomRight`
    * `selectBottomLeft`
    * `selectCorners`
    * `selectRows`
    * `selectCols`
    * `selectInner`
    * `selectOuter`
    
* `calcLoopEnrichment` function for flexibly
calculating enrichment of interactions compared
to their local background.

* `adjustEnrichment` and `plotEnrichment` for
adjusting the loop enrichment to remove the
effect of loop size on enrichment and visualize
this correction across chosen parameters.


# mariner 0.2.0

Methods for pulling Hi-C pixels and matrices from
.hic files and storing them on-disk with HDF5Array
and DelayedArray.

## Overview of functionality

New or updated functions:

* `snapToBins()`
    * "snaps" ranges in `GInteractions` objects to
    their nearest bin boundary. Allows spanning
    of multiple bins.

* `pullHicPixels()` extracts contact frequency
from `.hic` files and returns an `InteractionMatrix`
object containing a matrix of Hi-C interactions (rows)
and samples (columns).
    * Includes `counts()` accessor for matrix.
    * Custom `show()` method.
    * `rbind()` and `cbind()` methods.

* `pullHicMatrices()` extracts submatrices of
contact frequency from `.hic` files and returns
an `InteractionArray` object containing a
4-dimensional array of Hi-C submatrices, rownames,
and colnames.
    * Includes `counts()` accessor for submatrices.
    * Custom `show()` method.
    * `rbind()` and `cbind()` methods.

* `pixelsToMatrices()` takes `GInteractions` containing
single pixels (i.e., each range represents one `binSize`)
and expands ranges such that there is a `buffer` of pixels
around each range.

* `changePixelRes()` takes a `GInteractions` object
containing pixels of interest and is resized to the `from`
resolution/binSize (if its not already). Then count
matrices are extracted for each interaction and `.hic` file
using the new `to` resolution. Count matrices are
aggregated by interactions with the supplied `aggFUN` and
a new pixel is selected with the supplied `selectFUN`.
Allows block processing for large datasets. The object
returned is a `GInteractions` object with the updated
pixel ranges along with a column containing the aggregated
min/max value for that pixel.

* `calcLoopEnrichment()` pulls Hi-C pixels and calculates
the enrichment over background returning a `DelayedMatrix`
of enrichment scores where rows are interactions and
columns are Hi-C files.

* Accessors for GInteractions objects such as `seqnames1()`,
`start1()`, `end1()`, `seqnames2()`, `start2()`, `end2()`.

# mariner 0.1.0

First pre-release of mariner functionality focused
on manipulating, clustering, and merging paired
interactions.

## Overview of functionality

* Conversion of paired-range data to `GInteractions`
with `as_ginteractions`/`makeGInteractionsFromDataFrame`

* Functions for manipulating `GInteractions` and `GRanges`
objects with `binPairs`, `binRanges`, `shiftRanges`.

* Functions for clustering and merging lists of
`GInteractions` objects with `mergePairs`.

* Extensions to `GInteractions` class with
`MergedGInteractions`, and `DelegatingGInteractions`.

* Accessor functions for `MergedGInteractions`:
    * `aggPairMcols`
        * Aggregate metadata columns of clustered
        interactions.
    * `getPairClusters`
        * Return interactions for each cluster of
        interactions.
    * `selectionMethod`
        * Method used to select pairs from each
        cluster of interactions.
    * `sources`
        * List of names (or indices) used as input
        for clustering and merging.
    * `subsetBySource`
        * Return interactions unique to each source
        or combination of sources.
        
