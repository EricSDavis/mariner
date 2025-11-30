# Calculate loop enrichment over background.

Pulls Hi-C pixels and calculates the enrichment of the selected
foreground (\`fg\`) over the selected background (\`bg\`).

## Usage

``` r
calcLoopEnrichment(
  x,
  files,
  fg = selectCenterPixel(mhDist = 1, buffer = defaultBuffer()),
  bg = selectTopLeft(n = 4, buffer = defaultBuffer()) + selectBottomRight(n = 4, buffer =
    defaultBuffer()),
  FUN = function(fg, bg) median(fg + 1)/median(bg + 1),
  nBlocks = 5,
  verbose = TRUE,
  BPPARAM = bpparam(),
  ...
)

# S4 method for class 'GInteractions,character'
calcLoopEnrichment(
  x,
  files,
  fg = selectCenterPixel(mhDist = 1, buffer = defaultBuffer()),
  bg = selectTopLeft(n = 4, buffer = defaultBuffer()) + selectBottomRight(n = 4, buffer =
    defaultBuffer()),
  FUN = function(fg, bg) median(fg + 1)/median(bg + 1),
  nBlocks = 5,
  verbose = TRUE,
  BPPARAM = bpparam(),
  ...
)

# S4 method for class 'InteractionArray,missing'
calcLoopEnrichment(
  x,
  files,
  fg = selectCenterPixel(mhDist = 1, buffer = defaultBuffer()),
  bg = selectTopLeft(n = 4, buffer = defaultBuffer()) + selectBottomRight(n = 4, buffer =
    defaultBuffer()),
  FUN = function(fg, bg) median(fg + 1)/median(bg + 1),
  nBlocks = 5,
  verbose = TRUE,
  BPPARAM = bpparam(),
  ...
)
```

## Arguments

- x:

  GInteractions object or an InteractionArray object.

- files:

  Character file paths to \`.hic\` files. Required only if GInteractions
  object is supplied for x.

- fg:

  MatrixSelection object of matrix indices for the foreground.

- bg:

  MatrixSelection object of matrix indices for the background.

- FUN:

  Function with at least two parameters (i.e., \`fg\`, \`bg\`) defining
  how enrichment should be calculated. Must produce a single value
  (numeric of length one). The first and second parameters must
  represent fg and bg, respectively.

- nBlocks:

  Number of blocks for block-processing arrays. Default is 5. Increase
  this for large datasets. To read and process all data at once, set
  this value to 1.

- verbose:

  Boolean (TRUE or FALSE) describing whether to report block-processing
  progress.

- BPPARAM:

  Parallelization params (passed to \`BiocParallel::bplapply()\`).
  Default is the result of \`BiocParallel::bpparams()\`. Parallel
  processing is not available when \`by=interactions\`.

- ...:

  Additional arguments passed to \`pullHicMatrices\`. See
  ?\[\`pullHicMatrices\`\].

## Value

A DelayedMatrix of enrichment scores where rows are interactions (i.e.
loops) and columns are Hi-C files.

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

## Read in loops as GInteractions object
loops <-
    WT_5kbLoops.txt() |>
    setNames("WT") |>
    read.table(header=TRUE) |>
    as_ginteractions(keep.extra.columns=FALSE)
#> see ?marinerData and browseVignettes('marinerData') for documentation
#> loading from cache

## Removes the "chr" prefix for compatibility
## with the preprocessed hic files
GenomeInfoDb::seqlevelsStyle(loops) <- 'ENSEMBL'

## Expand binSize of loops
loops <- assignToBins(x=loops, binSize=100e3)

## Calculate loop enrichment
calcLoopEnrichment(x=loops[1:10],
                   files=hicFiles)
#> '0' = foreground;
#> 'X' = background;
#> '*' = both;
#> '-' = unselected
#>                                  
#>  X  X  X  X  -  -  -  -  -  -  - 
#>  X  X  X  X  -  -  -  -  -  -  - 
#>  X  X  X  X  -  -  -  -  -  -  - 
#>  X  X  X  X  -  -  -  -  -  -  - 
#>  -  -  -  -  -  0  -  -  -  -  - 
#>  -  -  -  -  0  0  0  -  -  -  - 
#>  -  -  -  -  -  0  -  -  -  -  - 
#>  -  -  -  -  -  -  -  X  X  X  X 
#>  -  -  -  -  -  -  -  X  X  X  X 
#>  -  -  -  -  -  -  -  X  X  X  X 
#>  -  -  -  -  -  -  -  X  X  X  X 
#> <10 x 2> DelayedMatrix object of type "double":
#>              FS        WT
#>  [1,] 1.0000000 1.6666667
#>  [2,] 1.0000000 1.4000000
#>  [3,] 1.5789474 2.1666667
#>  [4,] 1.3333333 1.0000000
#>  [5,] 0.7368421 0.7777778
#>  [6,] 1.3333333 1.0000000
#>  [7,] 0.7500000 1.0000000
#>  [8,] 1.7647059 1.4482759
#>  [9,] 0.8181818 0.9090909
#> [10,] 0.3333333 0.8000000

## Customize different foreground/background
## with selection functions
buffer <- 10 # choose pixel radius around center
fg <- selectCenterPixel(mhDist=seq(0,4), buffer=buffer)
bg <- selectCorners(n=6, buffer=buffer) +
    selectOuter(n=2, buffer=buffer)

## Calculate loop enrichment
calcLoopEnrichment(x=loops[1:10],
                   files=hicFiles,
                   fg=fg,
                   bg=bg)
#> '0' = foreground;
#> 'X' = background;
#> '*' = both;
#> '-' = unselected
#>                                                                
#>  X  X  X  X  X  X  X  X  X  X  X  X  X  X  X  X  X  X  X  X  X 
#>  X  X  X  X  X  X  X  X  X  X  X  X  X  X  X  X  X  X  X  X  X 
#>  X  X  X  X  X  X  -  -  -  -  -  -  -  -  -  X  X  X  X  X  X 
#>  X  X  X  X  X  X  -  -  -  -  -  -  -  -  -  X  X  X  X  X  X 
#>  X  X  X  X  X  X  -  -  -  -  -  -  -  -  -  X  X  X  X  X  X 
#>  X  X  X  X  X  X  -  -  -  -  -  -  -  -  -  X  X  X  X  X  X 
#>  X  X  -  -  -  -  -  -  -  -  0  -  -  -  -  -  -  -  -  X  X 
#>  X  X  -  -  -  -  -  -  -  0  0  0  -  -  -  -  -  -  -  X  X 
#>  X  X  -  -  -  -  -  -  0  0  0  0  0  -  -  -  -  -  -  X  X 
#>  X  X  -  -  -  -  -  0  0  0  0  0  0  0  -  -  -  -  -  X  X 
#>  X  X  -  -  -  -  0  0  0  0  0  0  0  0  0  -  -  -  -  X  X 
#>  X  X  -  -  -  -  -  0  0  0  0  0  0  0  -  -  -  -  -  X  X 
#>  X  X  -  -  -  -  -  -  0  0  0  0  0  -  -  -  -  -  -  X  X 
#>  X  X  -  -  -  -  -  -  -  0  0  0  -  -  -  -  -  -  -  X  X 
#>  X  X  -  -  -  -  -  -  -  -  0  -  -  -  -  -  -  -  -  X  X 
#>  X  X  X  X  X  X  -  -  -  -  -  -  -  -  -  X  X  X  X  X  X 
#>  X  X  X  X  X  X  -  -  -  -  -  -  -  -  -  X  X  X  X  X  X 
#>  X  X  X  X  X  X  -  -  -  -  -  -  -  -  -  X  X  X  X  X  X 
#>  X  X  X  X  X  X  -  -  -  -  -  -  -  -  -  X  X  X  X  X  X 
#>  X  X  X  X  X  X  X  X  X  X  X  X  X  X  X  X  X  X  X  X  X 
#>  X  X  X  X  X  X  X  X  X  X  X  X  X  X  X  X  X  X  X  X  X 
#> <10 x 2> DelayedMatrix object of type "double":
#>             FS       WT
#>  [1,] 2.500000 2.000000
#>  [2,] 2.666667 3.500000
#>  [3,] 2.000000 2.500000
#>  [4,] 2.000000 1.500000
#>  [5,] 2.666667 3.000000
#>  [6,] 3.000000 2.500000
#>  [7,] 1.500000 1.500000
#>  [8,] 6.000000 4.500000
#>  [9,] 2.333333 1.250000
#> [10,] 1.000000 1.000000

## Extract count matrices first
mats <- assignToBins(loops[1:10],100e3) |>
  pixelsToMatrices(buffer=10) |>
    pullHicMatrices(
    files=hicFiles,
    binSize=100e3)

## Calculate loop enrichment from count matrices 
calcLoopEnrichment(x = mats)
#> '0' = foreground;
#> 'X' = background;
#> '*' = both;
#> '-' = unselected
#>                                                                
#>  X  X  X  X  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#>  X  X  X  X  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#>  X  X  X  X  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#>  X  X  X  X  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  0  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  0  0  0  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  0  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  X  X  X  X 
#>  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  X  X  X  X 
#>  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  X  X  X  X 
#>  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  X  X  X  X 
#> <10 x 2> DelayedMatrix object of type "double":
#>              FS        WT
#>  [1,] 0.8333333 1.1111111
#>  [2,] 0.7142857 1.1666667
#>  [3,] 1.7647059 2.0000000
#>  [4,] 1.0000000 0.7500000
#>  [5,] 0.7000000 0.7000000
#>  [6,] 1.5000000 1.3333333
#>  [7,] 1.0000000 1.2000000
#>  [8,] 2.0689655 2.2105263
#>  [9,] 0.6206897 0.7142857
#> [10,] 0.5000000 1.0000000
```
