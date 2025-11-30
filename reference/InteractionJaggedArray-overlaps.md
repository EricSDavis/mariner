# Overlap methods for InteractionJaggedArray

Overlap methods for InteractionJaggedArray

## Usage

``` r
# S4 method for class 'InteractionJaggedArray,InteractionJaggedArray'
findOverlaps(
  query,
  subject,
  maxgap = -1L,
  minoverlap = 0L,
  type = c("any", "start", "end", "within", "equal"),
  select = c("all", "first", "last", "arbitrary"),
  ignore.strand = TRUE,
  ...,
  use.region = "both"
)

# S4 method for class 'InteractionJaggedArray,Vector'
findOverlaps(
  query,
  subject,
  maxgap = -1L,
  minoverlap = 0L,
  type = c("any", "start", "end", "within", "equal"),
  select = c("all", "first", "last", "arbitrary"),
  ignore.strand = TRUE,
  ...,
  use.region = "both"
)

# S4 method for class 'InteractionJaggedArray,missing'
findOverlaps(
  query,
  subject,
  maxgap = -1L,
  minoverlap = 0L,
  type = c("any", "start", "end", "within", "equal"),
  select = c("all", "first", "last", "arbitrary"),
  ignore.strand = TRUE,
  ...,
  use.region = "both"
)

# S4 method for class 'InteractionJaggedArray,InteractionJaggedArray'
countOverlaps(
  query,
  subject,
  maxgap = -1L,
  minoverlap = 0L,
  type = c("any", "start", "end", "within", "equal"),
  select = c("all", "first", "last", "arbitrary"),
  ignore.strand = TRUE,
  ...,
  use.region = "both"
)

# S4 method for class 'InteractionJaggedArray,Vector'
countOverlaps(
  query,
  subject,
  maxgap = -1L,
  minoverlap = 0L,
  type = c("any", "start", "end", "within", "equal"),
  select = c("all", "first", "last", "arbitrary"),
  ignore.strand = TRUE,
  ...,
  use.region = "both"
)

# S4 method for class 'InteractionJaggedArray,missing'
countOverlaps(
  query,
  subject,
  maxgap = -1L,
  minoverlap = 0L,
  type = c("any", "start", "end", "within", "equal"),
  select = c("all", "first", "last", "arbitrary"),
  ignore.strand = TRUE,
  ...,
  use.region = "both"
)

# S4 method for class 'InteractionJaggedArray,InteractionJaggedArray'
overlapsAny(
  query,
  subject,
  maxgap = -1L,
  minoverlap = 0L,
  type = c("any", "start", "end", "within", "equal"),
  ...,
  use.region = "both"
)

# S4 method for class 'InteractionJaggedArray,Vector'
overlapsAny(
  query,
  subject,
  maxgap = -1L,
  minoverlap = 0L,
  type = c("any", "start", "end", "within", "equal"),
  ...,
  use.region = "both"
)

# S4 method for class 'InteractionJaggedArray,missing'
overlapsAny(
  query,
  subject,
  maxgap = -1L,
  minoverlap = 0L,
  type = c("any", "start", "end", "within", "equal"),
  ...,
  use.region = "both"
)

# S4 method for class 'InteractionJaggedArray,InteractionJaggedArray'
subsetByOverlaps(
  x,
  ranges,
  maxgap = -1L,
  minoverlap = 0L,
  type = c("any", "start", "end", "within", "equal"),
  invert = FALSE,
  ...,
  use.region = "both"
)

# S4 method for class 'InteractionJaggedArray,Vector'
subsetByOverlaps(
  x,
  ranges,
  maxgap = -1L,
  minoverlap = 0L,
  type = c("any", "start", "end", "within", "equal"),
  invert = FALSE,
  ...,
  use.region = "both"
)

# S4 method for class 'InteractionJaggedArray,missing'
subsetByOverlaps(
  x,
  ranges,
  maxgap = -1L,
  minoverlap = 0L,
  type = c("any", "start", "end", "within", "equal"),
  invert = FALSE,
  ...,
  use.region = "both"
)
```

## Arguments

- query, subject, x, ranges:

  An InteractionJaggedArray, Vector, GInteractions or InteractionSet
  object, depending on the specified method. At least one of these must
  be a \`subject\` can be missing if query is an InteractionJaggedArray
  object.

- maxgap, minoverlap, type, select:

  see ?\`findOverlaps\` in the GenomicRanges package.

- ignore.strand:

  see ?\`findOverlaps\` in InteractionSet package for more information.

- ...:

  see ?\`findOverlaps\` in InteractionSet package for more information

- use.region:

  see ?\`findOverlaps\` in InteractionSet package for more information.

- invert:

  Boolean (TRUE/FALSE) to invert selection. Default is TRUE.

## Value

\`findOverlaps\` returns a Hits object. \`countOverlaps\` returns an
integer vector of overlaps for each interaction in \`query\`.
\`overlapsAny\` returns a logical vector of overlaps for each
interaction in \`query\`. \`subsetByOverlaps\` returns overlapping
interactions as an InteractionJaggedArray.

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

## Create test interactions
gi <- read.table(text="
            1 51000000 51300000 1 51000000 51500000
            2 52000000 52300000 3 52000000 52500000
            1 150000000 150500000 1 150000000 150300000
            2 52000000 52300000 2 52000000 52800000") |>
    as_ginteractions()

## InteractionJaggedArray object
iarr <- pullHicMatrices(gi, hicFiles, 100e03, half="both")

## Shift first two ranges out of range
gi2 <- c(assignToBins(gi[1:2], binSize=100e3, pos1=-200e3), gi[3:4])

## Find overlaps
findOverlaps(iarr, gi2)
#> Hits object with 2 hits and 0 metadata columns:
#>       queryHits subjectHits
#>       <integer>   <integer>
#>   [1]         3           3
#>   [2]         4           4
#>   -------
#>   queryLength: 4 / subjectLength: 4
countOverlaps(iarr, gi2)
#> [1] 0 0 1 1
countOverlaps(iarr, gi2, maxgap=100e3)
#> [1] 1 1 1 1
overlapsAny(iarr, gi2)
#> [1] FALSE FALSE  TRUE  TRUE
subsetByOverlaps(iarr, gi2)
#> class: InteractionJaggedArray
#> dim: 2 interaction(s), 2 file(s), variable count matrix(es)
#> metadata(3): binSize, norm, matrix
#> colData: FS, WT
#> colData names(2): files, fileNames
#> HDF5: /tmp/Rtmp9FSwoA/file7b572bccf993.h5
subsetByOverlaps(iarr, gi2, invert=TRUE)
#> class: InteractionArray 
#> dim: 2 interaction(s), 2 file(s), 3x5 count matrix(es)
#> metadata(3): binSize norm matrix
#> assays(1): counts
#> rownames: NULL
#> rowData names(0):
#> colnames(2): FS WT
#> colData names(2): files fileNames
#> type: GInteractions
#> regions: 7
```
