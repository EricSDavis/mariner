# Remove interactions that would cross the Hi-C diagonal or a specified distance from the diagonal.

Removes short interactions with some padding from the diagonal. If you
are resizing the regions with a function like \`pixelsToMatrices()\`,
make sure this function is run afterwards.

## Usage

``` r
removeShortPairs(x, padding = 0)

# S4 method for class 'GInteractions'
removeShortPairs(x, padding = 0)
```

## Arguments

- x:

  A GInteractions object.

- padding:

  Minimum distance away from the diagonal.

## Value

A GInteractions object with the short pairs removed.

## Details

Note this is only applies to intrachromosomal pairs, as pair distance is
meaningless for interchromosomal pairs. Therefore, all interchromosomal
pairs are kept.

## Examples

``` r
## Example GInteractions object
gi <- as_ginteractions(read.table(
    text="
        seqnames1 start1 end1 seqnames2 start2 end2 keep
        chr1 300 400 chr1 300 400 'no'
        chr1 100 200 chr1 300 400 'yes'
        chr1 300 400 chr1 100 200 'yes'
        chr1 300 400 chr2 300 400 'yes'
        chr1 250 350 chr1 300 400 'only_with_padding_50'
        chr1 300 400 chr1 250 350 'only_with_padding_50'
        ",
    header=TRUE
))

## Remove pairs that would cross the diagonal
removeShortPairs(gi)
#> GInteractions object with 3 interactions and 1 metadata column:
#>       seqnames1   ranges1     seqnames2   ranges2 |        keep
#>           <Rle> <IRanges>         <Rle> <IRanges> | <character>
#>   [1]      chr1   100-200 ---      chr1   300-400 |         yes
#>   [2]      chr1   300-400 ---      chr1   100-200 |         yes
#>   [3]      chr1   300-400 ---      chr2   300-400 |         yes
#>   -------
#>   regions: 4 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths

## Add 50bp of padding
removeShortPairs(gi, padding=50)
#> GInteractions object with 3 interactions and 1 metadata column:
#>       seqnames1   ranges1     seqnames2   ranges2 |        keep
#>           <Rle> <IRanges>         <Rle> <IRanges> | <character>
#>   [1]      chr1   100-200 ---      chr1   300-400 |         yes
#>   [2]      chr1   300-400 ---      chr1   100-200 |         yes
#>   [3]      chr1   300-400 ---      chr2   300-400 |         yes
#>   -------
#>   regions: 4 ranges and 0 metadata columns
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
```
