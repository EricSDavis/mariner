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
`MergedGInteractions`, `BinnedGInteractions`, and
`DelegatingGInteractions`.

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
        
