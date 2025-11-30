# Expand pixels to submatrices

Pixels are defined as paired-ranges with starts & ends equal to their
\`binSize\`. This function takes GInteractions fitting this description
and expands the ranges such that there is a \`buffer\` of pixels around
each range.

## Usage

``` r
pixelsToMatrices(x, buffer)

# S4 method for class 'GInteractions,numeric'
pixelsToMatrices(x, buffer)
```

## Arguments

- x:

  GInteractions object.

- buffer:

  Number (length one numeric vector) of pixels around the pixels in
  \`x\`.

## Value

\`x\` with updated ranges.

## Details

For example, a buffer of 3 would return a GInteractions object with 3
pixels surrounding the original pixel ranges.

After using \`pullHicMatrices()\`, the result will return a matrix of
row and column dimensions of buffer\*2+1.

Note, this function does not handle out-of-bound ranges.

## Examples

``` r
## Define example 100bp pixel
library(InteractionSet)
pixel <- GInteractions(
    anchor1=GRanges("chr1:500-600"),
    anchor2=GRanges("chr1:2000-2100")
)

## Expand pixel to matrix with
## 3 pixels surrounding the center
## pixel
region <- pixelsToMatrices(x=pixel, buffer=3)
region
#> GInteractions object with 1 interaction and 0 metadata columns:
#>       seqnames1   ranges1     seqnames2   ranges2
#>           <Rle> <IRanges>         <Rle> <IRanges>
#>   [1]      chr1   200-900 ---      chr1 1700-2400
#>   -------
#>   regions: 2 ranges and 0 metadata columns
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
