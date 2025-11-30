# Package index

## About mariner

- [`mariner`](http://ericscottdavis.com/mariner/reference/mariner-package.md)
  [`mariner-package`](http://ericscottdavis.com/mariner/reference/mariner-package.md)
  : Mariner: Explore the Hi-Cs

## Manipulating Paired Ranges

Convert, bin, and shift paired genomic ranges

- [`as_ginteractions()`](http://ericscottdavis.com/mariner/reference/as_ginteractions.md)
  [`makeGInteractionsFromDataFrame()`](http://ericscottdavis.com/mariner/reference/as_ginteractions.md)
  : Convert DataFrames to GInteraction objects
- [`shiftRanges()`](http://ericscottdavis.com/mariner/reference/shiftRanges.md)
  : Flexibly shifting GRanges according to strand
- [`binRanges()`](http://ericscottdavis.com/mariner/reference/binRanges.md)
  : Flexibly bin ranges
- [`assignToBins()`](http://ericscottdavis.com/mariner/reference/assignToBins.md)
  : Flexibly bin paired ranges
- [`snapToBins()`](http://ericscottdavis.com/mariner/reference/snapToBins.md)
  : Snap GRanges or GInteractions to nearest bins
- [`seqnames1()`](http://ericscottdavis.com/mariner/reference/GInteractions-accessors.md)
  [`seqnames2()`](http://ericscottdavis.com/mariner/reference/GInteractions-accessors.md)
  [`start1()`](http://ericscottdavis.com/mariner/reference/GInteractions-accessors.md)
  [`end1()`](http://ericscottdavis.com/mariner/reference/GInteractions-accessors.md)
  [`start2()`](http://ericscottdavis.com/mariner/reference/GInteractions-accessors.md)
  [`end2()`](http://ericscottdavis.com/mariner/reference/GInteractions-accessors.md)
  : Access each portion of a GInteractions-like object
- [`makeRandomGRanges()`](http://ericscottdavis.com/mariner/reference/makeRandomGInteractions.md)
  [`makeRandomGInteractions()`](http://ericscottdavis.com/mariner/reference/makeRandomGInteractions.md)
  : Creating random GRanges & GInteractions

## Clustering & Merging Interactions

Group nearby interactions and select one as representative

- [`mergePairs()`](http://ericscottdavis.com/mariner/reference/mergePairs.md)
  : Merge sets of paired interactions
- [`MergedGInteractions-class`](http://ericscottdavis.com/mariner/reference/MergedGInteractions-class.md)
  [`MergedGInteractions`](http://ericscottdavis.com/mariner/reference/MergedGInteractions-class.md)
  : MergedGInteractions Class
- [`aggMetadata()`](http://ericscottdavis.com/mariner/reference/aggMetadata.md)
  : Aggregate the metadata columns of merged pairs
- [`clusters()`](http://ericscottdavis.com/mariner/reference/clusters.md)
  : Get clustered pairs from MergedGInteractions object
- [`sources()`](http://ericscottdavis.com/mariner/reference/sources.md)
  : Accessor for sources
- [`sets()`](http://ericscottdavis.com/mariner/reference/sets.md) : Get
  each set from a MergedGInteractions object

## Extracting & Aggregating Interactions

Pull Hi-C pixels or matrices, then aggregate by files or interactions

- [`pullHicPixels()`](http://ericscottdavis.com/mariner/reference/pullHicPixels.md)
  : Pull contact frequency from \`.hic\` files
- [`InteractionMatrix()`](http://ericscottdavis.com/mariner/reference/InteractionMatrix-class.md)
  [`show(`*`<InteractionMatrix>`*`)`](http://ericscottdavis.com/mariner/reference/InteractionMatrix-class.md)
  [`rbind(`*`<InteractionMatrix>`*`)`](http://ericscottdavis.com/mariner/reference/InteractionMatrix-class.md)
  [`cbind(`*`<InteractionMatrix>`*`)`](http://ericscottdavis.com/mariner/reference/InteractionMatrix-class.md)
  : InteractionMatrix Class
- [`pixelsToMatrices()`](http://ericscottdavis.com/mariner/reference/pixelsToMatrices.md)
  : Expand pixels to submatrices
- [`removeShortPairs()`](http://ericscottdavis.com/mariner/reference/removeShortPairs.md)
  : Remove interactions that would cross the Hi-C diagonal or a
  specified distance from the diagonal.
- [`pullHicMatrices()`](http://ericscottdavis.com/mariner/reference/pullHicMatrices.md)
  : Pull submatrices from \`.hic\` files
- [`InteractionArray()`](http://ericscottdavis.com/mariner/reference/InteractionArray-class.md)
  [`show(`*`<InteractionArray>`*`)`](http://ericscottdavis.com/mariner/reference/InteractionArray-class.md)
  [`rbind(`*`<InteractionArray>`*`)`](http://ericscottdavis.com/mariner/reference/InteractionArray-class.md)
  [`cbind(`*`<InteractionArray>`*`)`](http://ericscottdavis.com/mariner/reference/InteractionArray-class.md)
  : InteractionArray Class
- [`counts(`*`<InteractionArray>`*`)`](http://ericscottdavis.com/mariner/reference/counts.md)
  [`counts(`*`<InteractionMatrix>`*`)`](http://ericscottdavis.com/mariner/reference/counts.md)
  [`` `counts<-`( ``*`<InteractionMatrix>`*`)`](http://ericscottdavis.com/mariner/reference/counts.md)
  : Access count matrices from InteractionArray or InteractionMatrix
- [`path(`*`<InteractionMatrix>`*`)`](http://ericscottdavis.com/mariner/reference/path.md)
  [`` `path<-`( ``*`<InteractionMatrix>`*`)`](http://ericscottdavis.com/mariner/reference/path.md)
  : Accessor for h5File path from an InteractionMatrix
- [`show(`*`<InteractionJaggedArray>`*`)`](http://ericscottdavis.com/mariner/reference/InteractionJaggedArray-class.md)
  [`dim(`*`<InteractionJaggedArray>`*`)`](http://ericscottdavis.com/mariner/reference/InteractionJaggedArray-class.md)
  [`interactions(`*`<InteractionJaggedArray>`*`)`](http://ericscottdavis.com/mariner/reference/InteractionJaggedArray-class.md)
  [`metadata(`*`<InteractionJaggedArray>`*`)`](http://ericscottdavis.com/mariner/reference/InteractionJaggedArray-class.md)
  [`colData(`*`<InteractionJaggedArray>`*`)`](http://ericscottdavis.com/mariner/reference/InteractionJaggedArray-class.md)
  [`counts(`*`<InteractionJaggedArray>`*`)`](http://ericscottdavis.com/mariner/reference/InteractionJaggedArray-class.md)
  [`path(`*`<InteractionJaggedArray>`*`)`](http://ericscottdavis.com/mariner/reference/InteractionJaggedArray-class.md)
  [`length(`*`<InteractionJaggedArray>`*`)`](http://ericscottdavis.com/mariner/reference/InteractionJaggedArray-class.md)
  [`` `[`( ``*`<InteractionJaggedArray>`*`,`*`<ANY>`*`,`*`<ANY>`*`,`*`<ANY>`*`)`](http://ericscottdavis.com/mariner/reference/InteractionJaggedArray-class.md)
  : InteractionJaggedArray Class
- [`show(`*`<JaggedArray>`*`)`](http://ericscottdavis.com/mariner/reference/JaggedArray-class.md)
  [`` `[`( ``*`<JaggedArray>`*`,`*`<ANY>`*`,`*`<ANY>`*`,`*`<ANY>`*`)`](http://ericscottdavis.com/mariner/reference/JaggedArray-class.md)
  [`as.list(`*`<JaggedArray>`*`)`](http://ericscottdavis.com/mariner/reference/JaggedArray-class.md)
  [`path(`*`<JaggedArray>`*`)`](http://ericscottdavis.com/mariner/reference/JaggedArray-class.md)
  [`dim(`*`<JaggedArray>`*`)`](http://ericscottdavis.com/mariner/reference/JaggedArray-class.md)
  : JaggedArray Class
- [`findOverlaps(`*`<InteractionJaggedArray>`*`,`*`<InteractionJaggedArray>`*`)`](http://ericscottdavis.com/mariner/reference/InteractionJaggedArray-overlaps.md)
  [`findOverlaps(`*`<InteractionJaggedArray>`*`,`*`<Vector>`*`)`](http://ericscottdavis.com/mariner/reference/InteractionJaggedArray-overlaps.md)
  [`findOverlaps(`*`<InteractionJaggedArray>`*`,`*`<missing>`*`)`](http://ericscottdavis.com/mariner/reference/InteractionJaggedArray-overlaps.md)
  [`countOverlaps(`*`<InteractionJaggedArray>`*`,`*`<InteractionJaggedArray>`*`)`](http://ericscottdavis.com/mariner/reference/InteractionJaggedArray-overlaps.md)
  [`countOverlaps(`*`<InteractionJaggedArray>`*`,`*`<Vector>`*`)`](http://ericscottdavis.com/mariner/reference/InteractionJaggedArray-overlaps.md)
  [`countOverlaps(`*`<InteractionJaggedArray>`*`,`*`<missing>`*`)`](http://ericscottdavis.com/mariner/reference/InteractionJaggedArray-overlaps.md)
  [`overlapsAny(`*`<InteractionJaggedArray>`*`,`*`<InteractionJaggedArray>`*`)`](http://ericscottdavis.com/mariner/reference/InteractionJaggedArray-overlaps.md)
  [`overlapsAny(`*`<InteractionJaggedArray>`*`,`*`<Vector>`*`)`](http://ericscottdavis.com/mariner/reference/InteractionJaggedArray-overlaps.md)
  [`overlapsAny(`*`<InteractionJaggedArray>`*`,`*`<missing>`*`)`](http://ericscottdavis.com/mariner/reference/InteractionJaggedArray-overlaps.md)
  [`subsetByOverlaps(`*`<InteractionJaggedArray>`*`,`*`<InteractionJaggedArray>`*`)`](http://ericscottdavis.com/mariner/reference/InteractionJaggedArray-overlaps.md)
  [`subsetByOverlaps(`*`<InteractionJaggedArray>`*`,`*`<Vector>`*`)`](http://ericscottdavis.com/mariner/reference/InteractionJaggedArray-overlaps.md)
  [`subsetByOverlaps(`*`<InteractionJaggedArray>`*`,`*`<missing>`*`)`](http://ericscottdavis.com/mariner/reference/InteractionJaggedArray-overlaps.md)
  : Overlap methods for InteractionJaggedArray
- [`regularize()`](http://ericscottdavis.com/mariner/reference/regularize.md)
  : Regularize JaggedArray or InteractionJaggedArray objects
- [`aggHicMatrices()`](http://ericscottdavis.com/mariner/reference/aggHicMatrices.md)
  : Aggregate count matrices from InteractionArray objects
- [`hdf5BlockApply()`](http://ericscottdavis.com/mariner/reference/hdf5BlockApply.md)
  : HDF5-backed blockApply
- [`plotMatrix()`](http://ericscottdavis.com/mariner/reference/plotMatrix.md)
  : Plot matrix

## Pileup analysis

Wraps several mariner functions to facilitate pileup analysis of Hi-C
pixels, domains or boundaries.

- [`pileupPixels()`](http://ericscottdavis.com/mariner/reference/pileupPixels.md)
  : Pileup Hi-C pixels
- [`pileupDomains()`](http://ericscottdavis.com/mariner/reference/pileupDomains.md)
  : Pileup Hi-C domains
- [`pileupBoundaries()`](http://ericscottdavis.com/mariner/reference/pileupBoundaries.md)
  : Pileup Hi-C contacts around boundary regions

## Calculating Loop Enrichment

Determine loop enrichment to local background with selection functions
to flexibility select foreground and background.

- [`calcLoopEnrichment()`](http://ericscottdavis.com/mariner/reference/calcLoopEnrichment.md)
  : Calculate loop enrichment over background.
- [`plotEnrichment()`](http://ericscottdavis.com/mariner/reference/adjustEnrichment.md)
  [`adjustEnrichment()`](http://ericscottdavis.com/mariner/reference/adjustEnrichment.md)
  : Adjust loop enrichment to remove distance- dependent effect.

### Selection functions

- [`MatrixSelection-class`](http://ericscottdavis.com/mariner/reference/MatrixSelection-class.md)
  [`MatrixSelection`](http://ericscottdavis.com/mariner/reference/MatrixSelection-class.md)
  : MatrixSelection Class
- [`selectRadius()`](http://ericscottdavis.com/mariner/reference/selection-functions.md)
  [`selectCenterPixel()`](http://ericscottdavis.com/mariner/reference/selection-functions.md)
  [`selectSubmatrix()`](http://ericscottdavis.com/mariner/reference/selection-functions.md)
  [`selectCoordinates()`](http://ericscottdavis.com/mariner/reference/selection-functions.md)
  [`selectBlock()`](http://ericscottdavis.com/mariner/reference/selection-functions.md)
  [`selectTopLeft()`](http://ericscottdavis.com/mariner/reference/selection-functions.md)
  [`selectTopRight()`](http://ericscottdavis.com/mariner/reference/selection-functions.md)
  [`selectBottomRight()`](http://ericscottdavis.com/mariner/reference/selection-functions.md)
  [`selectBottomLeft()`](http://ericscottdavis.com/mariner/reference/selection-functions.md)
  [`selectCorners()`](http://ericscottdavis.com/mariner/reference/selection-functions.md)
  [`selectRows()`](http://ericscottdavis.com/mariner/reference/selection-functions.md)
  [`selectCols()`](http://ericscottdavis.com/mariner/reference/selection-functions.md)
  [`selectInner()`](http://ericscottdavis.com/mariner/reference/selection-functions.md)
  [`selectOuter()`](http://ericscottdavis.com/mariner/reference/selection-functions.md)
  [`show(`*`<MatrixSelection>`*`)`](http://ericscottdavis.com/mariner/reference/selection-functions.md)
  : Visualize selection for a MatrixSelection object
- [`selectionMethod()`](http://ericscottdavis.com/mariner/reference/selectionMethod.md)
  : Get selectionMethod from MergedGInteractions object
- [`defaultBuffer()`](http://ericscottdavis.com/mariner/reference/defaultBuffer.md)
  : Return default buffer If InteractionArray is supplied, it uses the
  dimensions of counts matrices to set the buffer dimensions.

## Miscellaneous

- [`changePixelRes()`](http://ericscottdavis.com/mariner/reference/changePixelRes.md)
  : Change pixels from one resolution to another selecting the new pixel
  using Hi-C data.
- [`selectPixel()`](http://ericscottdavis.com/mariner/reference/selectPixel.md)
  : Get the pixel representing the strongest or weakest interaction in
  an InteractionArray
