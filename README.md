# Mariner: Explore the Hi-C's <img src="man/figures/mariner.png" id="mariner_logo" align="right" width="125"/>

`mariner` is an R package with a collection of useful tools for performing Hi-C data analysis. The tools fall into the following broad categories:

-   Manipulating/Merging Anchors

    A collection of functions for converting to `GInteractions` objects and manipulating anchors.
    
-   Generating Pairs from Ranges

    Methods to easily generate combinations of paired interactions.
  
-   Extracting/Aggregating Interactions

    Functions that are optimized to quickly extract and aggregate data from `.hic` files.

-   Visualization

    Additional visualization functions that are compatible with the `plotgardener` package.

## Installation

This package can be installed via github:

```{r}
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("EricSDavis/mariner")
```
