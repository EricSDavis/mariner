# Flexibly shifting GRanges according to strand

Flexibly shifting GRanges according to strand

## Usage

``` r
shiftRanges(x, pos)

# S4 method for class 'GRanges,character_OR_numeric'
shiftRanges(x, pos)
```

## Arguments

- x:

  GRanges object

- pos:

  Position within anchors to resize the bin. Can be a character or
  integer vector of length 1 or \`length(x)\` designating the position
  for each element in bedpe. Character options are "start", "end" and
  "center". Integers are referenced from the start position for '+' and
  '\*' strands and from the end position for the '-' strand.

## Value

GRanges object with a single position range that has been shifted
appropriately.

## Examples

``` r
library(GenomicRanges)

## Create example GRanges
gr1 <- GRanges(seqnames = "chr1",
               ranges = IRanges::IRanges(start = rep(5000,3),
                                         end = rep(6000,3)),
               strand = c('+', '-', '*'))

gr2 <- gr1 |> promoters(upstream = 2000, downstream = 200)

## Shifting anchors by keyword
shiftRanges(gr1, 'start')
#> GRanges object with 3 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chr1      5000      +
#>   [2]     chr1      6000      -
#>   [3]     chr1      5000      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
shiftRanges(gr1, 'end')
#> GRanges object with 3 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chr1      6000      +
#>   [2]     chr1      5000      -
#>   [3]     chr1      6000      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
shiftRanges(gr1, 'center')
#> GRanges object with 3 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chr1      5500      +
#>   [2]     chr1      5500      -
#>   [3]     chr1      5500      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

## Shifting anchors by position
shiftRanges(gr1, 100)
#> GRanges object with 3 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chr1      5100      +
#>   [2]     chr1      5900      -
#>   [3]     chr1      5100      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
shiftRanges(gr1, c(100, 200, 300))
#> GRanges object with 3 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chr1      5100      +
#>   [2]     chr1      5800      -
#>   [3]     chr1      5300      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

## Shifting back to TSS
shiftRanges(gr2, 2000)
#> GRanges object with 3 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chr1      5000      +
#>   [2]     chr1      6000      -
#>   [3]     chr1      5000      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
