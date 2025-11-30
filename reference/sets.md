# Get each set from a MergedGInteractions object

Returns the subset of MergedGInteractions that belong to each input
source object (see these with \`sources(x)\`). If the source pairs all
come from the same object, their corresponding merged pair is returned.
However, if at least one source pair comes from a different object, then
that merged pair is not returned.

## Usage

``` r
sets(x, include, exclude)

# S4 method for class 'MergedGInteractions,missing,missing'
sets(x)

# S4 method for class 'MergedGInteractions,character_OR_missing,missing'
sets(x, include)

# S4 method for class 'MergedGInteractions,missing,character_OR_missing'
sets(x, exclude)

# S4 method for class 'MergedGInteractions,character_OR_missing,character_OR_missing'
sets(x, include, exclude)
```

## Arguments

- x:

  MergedGInteractions object.

- include:

  (Optional) A character vector of sources in which a pair must be
  present. For a list of available sources use \`sources(x)\`.

- exclude:

  (Optional) A character vector of sources in which a pair must be
  absent. For a list of available sources use \`sources(x)\`.

## Value

A list of subsetted \`MergedGInteractions\` objects or a
\`MergedGInteractions\` object (if \`include\` and/or \`exclude\` are
used).

## Details

Optional \`include\` and \`exclude\` parameters modulate the behavior of
\`sets\` to return different subsets of originating pairs. For example,
\`include\` requires that the returned pairs be present in specific
sources, while \`exclude\` requires that returned pairs be absent from
specific sources. Sources not listed in either \`include\` or
\`exclude\` are ignored (they may or may not) be present in the returned
\`MergedGInteractions\` object. \`include\` and \`exclude\` can be used
indepedently or in combination to return every possible set. If any of
the same sources are used in both \`include\` and \`exclude\` the
function will return a 0-length MergedGInteractions object.

## Examples

``` r
## Load required packages
library(GenomicRanges)
library(InteractionSet)

## Define example anchor regions
gr1 <-
    GRanges(seqnames = "chr1",
            ranges = IRanges(start = c(30,40,40,70,80),
                             end = c(40,50,50,80,90)))
gr2 <-
    GRanges(seqnames = "chr1",
            ranges = IRanges(start = c(30,30,50,10,30),
                             end = c(40,40,60,20,40)))

## Form GInteractions and split into two files
giList <- split(x = GInteractions(gr1, gr2),
                f = c(rep(1,3), rep(2,2)))

## Merge pairs
x <- mergePairs(x = giList, radius = 20)

sets(x)
#> $`1`
#> MergedGInteractions object with 1 interaction and 0 metadata columns:
#>       seqnames1   ranges1     seqnames2   ranges2
#>           <Rle> <IRanges>         <Rle> <IRanges>
#>   [1]      chr1     40-50 ---      chr1     30-40
#>   -------
#>   regions: 5 ranges and 0 metadata columns
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
#> 
#> $`2`
#> MergedGInteractions object with 2 interactions and 0 metadata columns:
#>       seqnames1   ranges1     seqnames2   ranges2
#>           <Rle> <IRanges>         <Rle> <IRanges>
#>   [1]      chr1     70-80 ---      chr1     10-20
#>   [2]      chr1     80-90 ---      chr1     30-40
#>   -------
#>   regions: 5 ranges and 0 metadata columns
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
#> 
#> $`1_2`
#> MergedGInteractions object with 0 interactions and 0 metadata columns:
#>    seqnames1   ranges1     seqnames2   ranges2
#>        <Rle> <IRanges>         <Rle> <IRanges>
#>   -------
#>   regions: 5 ranges and 0 metadata columns
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
#> 
```
