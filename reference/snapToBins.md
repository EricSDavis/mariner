# Snap GRanges or GInteractions to nearest bins

Snap GRanges or GInteractions to nearest bins

Snap paired-objects to nearest bins

## Usage

``` r
snapToBins(x, binSize)

# S4 method for class 'GRanges,numeric'
snapToBins(x, binSize)

# S4 method for class 'GInteractions,numeric'
snapToBins(x, binSize)
```

## Arguments

- x:

  \`GInteractions\` object.

- binSize:

  Integer (numeric) describing the new size of each range.

## Value

GRanges object snapped to the nearest \`binSize\`.

Input object snapped to the nearest \`binSize\`.

## Examples

``` r
library(GenomicRanges)
## Example GRanges object
x <- GRanges(seqnames = c("chr1"),
             ranges = IRanges(start = c(1, 1, 25, 19, 21),
                              end = c(15, 11, 31, 31, 39)))

snapToBins(x, binSize = 5)
#> GRanges object with 5 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chr1      0-15      *
#>   [2]     chr1      0-10      *
#>   [3]     chr1     25-30      *
#>   [4]     chr1     20-30      *
#>   [5]     chr1     20-40      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
snapToBins(x, binSize = 10)
#> GRanges object with 5 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chr1      0-20      *
#>   [2]     chr1      0-10      *
#>   [3]     chr1     20-30      *
#>   [4]     chr1     20-30      *
#>   [5]     chr1     20-40      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
snapToBins(x, binSize = 20)
#> GRanges object with 5 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chr1      0-20      *
#>   [2]     chr1      0-20      *
#>   [3]     chr1     20-40      *
#>   [4]     chr1     20-40      *
#>   [5]     chr1     20-40      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

library(InteractionSet)
## Sample GInteractions object
x <- GInteractions(anchor1 = c(GRanges("chr1:1-15"),
                               GRanges("chr1:1-11")),
                   anchor2 = c(GRanges("chr1:25-31"),
                               GRanges("chr1:19-31")))

snapToBins(x, binSize = 5)
#> GInteractions object with 2 interactions and 0 metadata columns:
#>       seqnames1   ranges1     seqnames2   ranges2
#>           <Rle> <IRanges>         <Rle> <IRanges>
#>   [1]      chr1      0-15 ---      chr1     25-30
#>   [2]      chr1      0-10 ---      chr1     20-30
#>   -------
#>   regions: 4 ranges and 0 metadata columns
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
snapToBins(x, binSize = 10)
#> GInteractions object with 2 interactions and 0 metadata columns:
#>       seqnames1   ranges1     seqnames2   ranges2
#>           <Rle> <IRanges>         <Rle> <IRanges>
#>   [1]      chr1      0-20 ---      chr1     20-30
#>   [2]      chr1      0-10 ---      chr1     20-30
#>   -------
#>   regions: 3 ranges and 0 metadata columns
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
snapToBins(x, binSize = 20)
#> GInteractions object with 2 interactions and 0 metadata columns:
#>       seqnames1   ranges1     seqnames2   ranges2
#>           <Rle> <IRanges>         <Rle> <IRanges>
#>   [1]      chr1      0-20 ---      chr1     20-40
#>   [2]      chr1      0-20 ---      chr1     20-40
#>   -------
#>   regions: 2 ranges and 0 metadata columns
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
