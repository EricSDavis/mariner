# Flexibly bin ranges

Flexibly bin ranges

## Usage

``` r
binRanges(x, binSize, pos = "center")

# S4 method for class 'GRanges,numeric,character_OR_numeric_OR_missing'
binRanges(x, binSize, pos = "center")
```

## Arguments

- x:

  \`GRanges\` object

- binSize:

  Integer (numeric) describing the new size of each range.

- pos:

  Position within range to resize the bin. Can be a character or integer
  vector of length 1 or \`length(x)\` designating the position for each
  element in \`x\`. Character options are "start", "end" and "center".
  Integers are referenced from the start position for '+' and '\*'
  strands and from the end position for the '-' strand.

## Value

\`GRanges\` object that has been shifted by \`pos\` and assigned to bins
of \`binSize\`.

## Examples

``` r
library(GenomicRanges)

## Create example GRanges
gr1 <- GRanges(seqnames = "chr1",
               ranges = IRanges::IRanges(start = rep(5000,3),
                                         end = rep(6000,3)),
               strand = c('+', '-', '*'))

gr2 <- gr1 |> promoters(upstream = 2000, downstream = 200)

## Binning the results
binRanges(x = gr1, binSize = 1000, pos = 'start')
#> GRanges object with 3 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chr1 5000-6000      +
#>   [2]     chr1 6000-7000      -
#>   [3]     chr1 5000-6000      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
binRanges(x = gr1, binSize = 1000, pos = 'end')
#> GRanges object with 3 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chr1 6000-7000      +
#>   [2]     chr1 5000-6000      -
#>   [3]     chr1 6000-7000      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
binRanges(x = gr1, binSize = 1000, pos = 'center')
#> GRanges object with 3 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chr1 5000-6000      +
#>   [2]     chr1 5000-6000      -
#>   [3]     chr1 5000-6000      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

## Bin after shifting back to TSS
binRanges(x = gr2, binSize = 1000, pos = 2000)
#> GRanges object with 3 ranges and 0 metadata columns:
#>       seqnames    ranges strand
#>          <Rle> <IRanges>  <Rle>
#>   [1]     chr1 5000-6000      +
#>   [2]     chr1 6000-7000      -
#>   [3]     chr1 5000-6000      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
