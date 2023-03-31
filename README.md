# Mariner: Explore the Hi-Cs <img src="man/figures/mariner.png" id="mariner_logo" align="right" width="150px" style="padding-left:20px; background-color:white"/>

<!-- badges: start -->
![GitHub R package version](https://img.shields.io/github/r-package/v/EricSDavis/mariner?style=plastic)
[![DOI](https://zenodo.org/badge/475953890.svg)](https://zenodo.org/badge/latestdoi/475953890)
<!-- badges: end -->

## Why mariner?

Disruption or aberrant formation of chromatin interactions can result in
developmental abnormalities and disease. Therefore, deriving biological
insights from 3D chromatin structure experiments, such as Hi-C or Micro-C,
is essential for understanding and correcting human disease.

`mariner` is an R/Bioconductor package for exploring Hi-C data. It enables
users to flexibly manipulate, extract, and aggregate chromatin interaction
data quickly and efficiently.

**One ecosystem:**
`mariner` extends common Bioconductor classes, leveraging the thousands of
existing tools for analyzing and visualizing genomic data.

**Modular design:**
`mariner's` functions can be combined and chained in various ways to produce
custom workflows.

**Fast and efficient:**
`mariner` leverages HDF5 to store large results and uses block processing
to minimize hardware requirements.

## Key features

**Extracting & Aggregating Interactions**

Pull Hi-C pixels or matrices, then aggregate by files or interactions

**Clustering & Merging Interactions**

Group nearby interactions and select one as representative

**Manipulating Pairs**

Convert, bin, and shift paired genomic ranges

**Calculating Loop Enrichment**

Determine loop enrichment to local background with
_selection functions_ to flexibiliy select foreground
and background.

## Installation

This package can be installed via github:

```{r}
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("EricSDavis/mariner")
```
