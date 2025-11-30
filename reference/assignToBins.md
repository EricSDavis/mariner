# Flexibly bin paired ranges

Paired range objects (like \`GInteractions\` or BEDPE-formatted
\`data.frame\`-like objects) can be binned separately for each set of
ranges.

## Usage

``` r
assignToBins(x, binSize, pos1 = "center", pos2 = "center", ...)

# S4 method for class 'DF_OR_df_OR_dt,numeric,character_OR_numeric_OR_missing,character_OR_numeric_OR_missing'
assignToBins(x, binSize, pos1, pos2)

# S4 method for class 'GInteractions,numeric,character_OR_numeric_OR_missing,character_OR_numeric_OR_missing'
assignToBins(x, binSize, pos1, pos2)
```

## Arguments

- x:

  \`GInteractions\` or \`data.frame\`-like object with paired
  interactions.

- binSize:

  Integer (numeric) vector describing the new size of each pair of
  ranges. Accepts up to 2 values for adjusting each pair.

- pos1, pos2:

  Position within anchors to resize the bin. Can be a character or
  integer vector of length 1 or \`length(x)\` designating the position
  for each element in \`x\`. Character options are "start", "end" and
  "center". Integers are referenced from the start position for '+' and
  '\*' strands and from the end position for the '-' strand.

- ...:

  Additional arguments.

## Value

GInteractions-like object binned to \`binSize\` by \`pos1\` and
\`pos2\`.

## Examples

``` r
## Construct interactions as data.frame
df1 <-
    data.frame(chr1 = "chr1", x1 = 10000, x2 = 20000,
               chr2 = "chr1", y1 = 30000, y2 = 40000)

## Assign each range to 20-kb bins from the start positions
assignToBins(x = df1,
         binSize = 20000,
         pos1 = 'start',
         pos2 = 'start')
#> GInteractions object with 1 interaction and 0 metadata columns:
#>       seqnames1   ranges1     seqnames2     ranges2
#>           <Rle> <IRanges>         <Rle>   <IRanges>
#>   [1]      chr1   0-20000 ---      chr1 20000-40000
#>   -------
#>   regions: 2 ranges and 0 metadata columns
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

## Construct GInteractions
library(InteractionSet)
gi1 <-
    data.frame(chr1 = "chr1", x1 = 10000, x2 = 20000,
               chr2 = "chr1", y1 = 30000, y2 = 40000) |>
    as_ginteractions()

## Assign each range to 20-kb bins from the start positions
assignToBins(x = gi1,
         binSize = 20000,
         pos1 = 'start',
         pos2 = 'start')
#> GInteractions object with 1 interaction and 0 metadata columns:
#>       seqnames1   ranges1     seqnames2     ranges2
#>           <Rle> <IRanges>         <Rle>   <IRanges>
#>   [1]      chr1   0-20000 ---      chr1 20000-40000
#>   -------
#>   regions: 2 ranges and 0 metadata columns
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
