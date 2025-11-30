# Convert DataFrames to GInteraction objects

\`as_ginteractions\` takes a paired-interaction (i.e. BEDPE) formatted
data-frame-like object and converts it to a GInteractions object. For
convenience, \`makeGInteractionsFromDataFrame\` can be used as an alias.

## Usage

``` r
as_ginteractions(
  df,
  keep.extra.columns = TRUE,
  starts.in.df.are.0based = FALSE,
  ...
)

makeGInteractionsFromDataFrame(
  df,
  keep.extra.columns = TRUE,
  starts.in.df.are.0based = FALSE,
  ...
)

# S4 method for class 'DF_OR_df_OR_dt,logical_OR_missing,logical_OR_missing'
makeGInteractionsFromDataFrame(df, keep.extra.columns, starts.in.df.are.0based)

# S4 method for class 'DF_OR_df_OR_dt,logical_OR_missing,logical_OR_missing'
as_ginteractions(df, keep.extra.columns, starts.in.df.are.0based)
```

## Arguments

- df:

  A data.table, data.frame, or DataFrame object. Assumes that the first
  6 colummns are in the format chr1, start1, end1 and chr2, start2,
  end2, representing each pair of interactions.

- keep.extra.columns:

  TRUE or FALSE (the default). If TRUE, the columns in df that are not
  used to form the genomic ranges of the returned GRanges object are
  then returned as metadata columns on the object. Otherwise, they are
  ignored. If df has a width column, then it's always ignored.

- starts.in.df.are.0based:

  TRUE or FALSE (the default). If TRUE, then the start positions of the
  genomic ranges in df are considered to be 0-based and are converted to
  1-based in the returned GRanges object. This feature is intended to
  make it more convenient to handle input that contains data obtained
  from resources using the "0-based start" convention. A notorious
  example of such resource is the UCSC Table Browser
  (http://genome.ucsc.edu/cgi-bin/hgTables).

- ...:

  Additional arguments.

## Value

GInteraction object

## Examples

``` r
## data.frame
df <- data.frame(chr1 = "chr1", x1 = 10000, x2 = 20000,
                 chr2 = "chr1", y1 = 30000, y2 = 40000)
makeGInteractionsFromDataFrame(df)
#> GInteractions object with 1 interaction and 0 metadata columns:
#>       seqnames1     ranges1     seqnames2     ranges2
#>           <Rle>   <IRanges>         <Rle>   <IRanges>
#>   [1]      chr1 10000-20000 ---      chr1 30000-40000
#>   -------
#>   regions: 2 ranges and 0 metadata columns
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

## data.frame
df <- data.frame(chr1 = "chr1", x1 = 10000, x2 = 20000,
                 chr2 = "chr1", y1 = 30000, y2 = 40000)
as_ginteractions(df)
#> GInteractions object with 1 interaction and 0 metadata columns:
#>       seqnames1     ranges1     seqnames2     ranges2
#>           <Rle>   <IRanges>         <Rle>   <IRanges>
#>   [1]      chr1 10000-20000 ---      chr1 30000-40000
#>   -------
#>   regions: 2 ranges and 0 metadata columns
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

## data.table
library(data.table)
df <- data.table::data.table(chr1 = "chr1", x1 = 10000, x2 = 20000,
                             chr2 = "chr1", y1 = 30000, y2 = 40000)
as_ginteractions(df)
#> GInteractions object with 1 interaction and 0 metadata columns:
#>       seqnames1     ranges1     seqnames2     ranges2
#>           <Rle>   <IRanges>         <Rle>   <IRanges>
#>   [1]      chr1 10000-20000 ---      chr1 30000-40000
#>   -------
#>   regions: 2 ranges and 0 metadata columns
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

## DataFrame
library(S4Vectors)
df <- DataFrame(chr1 = "chr1", x1 = 10000, x2 = 20000,
                chr2 = "chr1", y1 = 30000, y2 = 40000)
as_ginteractions(df)
#> GInteractions object with 1 interaction and 0 metadata columns:
#>       seqnames1     ranges1     seqnames2     ranges2
#>           <Rle>   <IRanges>         <Rle>   <IRanges>
#>   [1]      chr1 10000-20000 ---      chr1 30000-40000
#>   -------
#>   regions: 2 ranges and 0 metadata columns
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

## Alias
df <- data.frame(chr1 = "chr1", x1 = 10000, x2 = 20000,
                 chr2 = "chr1", y1 = 30000, y2 = 40000,
                 pval = 0.05, dist = 10000)
makeGInteractionsFromDataFrame(df)
#> GInteractions object with 1 interaction and 2 metadata columns:
#>       seqnames1     ranges1     seqnames2     ranges2 |      pval      dist
#>           <Rle>   <IRanges>         <Rle>   <IRanges> | <numeric> <numeric>
#>   [1]      chr1 10000-20000 ---      chr1 30000-40000 |      0.05     10000
#>   -------
#>   regions: 2 ranges and 0 metadata columns
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

## Additional metadata
df <- data.frame(chr1 = "chr1", x1 = 10000, x2 = 20000,
                 chr2 = "chr1", y1 = 30000, y2 = 40000,
                 pval = 0.05, dist = 10000)
as_ginteractions(df)
#> GInteractions object with 1 interaction and 2 metadata columns:
#>       seqnames1     ranges1     seqnames2     ranges2 |      pval      dist
#>           <Rle>   <IRanges>         <Rle>   <IRanges> | <numeric> <numeric>
#>   [1]      chr1 10000-20000 ---      chr1 30000-40000 |      0.05     10000
#>   -------
#>   regions: 2 ranges and 0 metadata columns
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

## Remove additional metadata
as_ginteractions(df, keep.extra.columns = FALSE)
#> GInteractions object with 1 interaction and 0 metadata columns:
#>       seqnames1     ranges1     seqnames2     ranges2
#>           <Rle>   <IRanges>         <Rle>   <IRanges>
#>   [1]      chr1 10000-20000 ---      chr1 30000-40000
#>   -------
#>   regions: 2 ranges and 0 metadata columns
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

## Add 1 to starts (for 0-based programs)
as_ginteractions(df, starts.in.df.are.0based = TRUE)
#> GInteractions object with 1 interaction and 2 metadata columns:
#>       seqnames1     ranges1     seqnames2     ranges2 |      pval      dist
#>           <Rle>   <IRanges>         <Rle>   <IRanges> | <numeric> <numeric>
#>   [1]      chr1 10001-20000 ---      chr1 30001-40000 |      0.05     10000
#>   -------
#>   regions: 2 ranges and 0 metadata columns
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
