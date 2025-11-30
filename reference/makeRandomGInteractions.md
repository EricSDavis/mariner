# Creating random GRanges & GInteractions

Creating random GRanges & GInteractions

## Usage

``` r
makeRandomGRanges(seqinfo, n = 100, ...)

makeRandomGInteractions(seqinfo, n = 100, interchromosomal = TRUE, ...)

# S4 method for class 'Seqinfo'
makeRandomGRanges(seqinfo, n, .rows = NULL)

# S4 method for class 'Seqinfo'
makeRandomGInteractions(seqinfo, n, interchromosomal)
```

## Arguments

- seqinfo:

  A Seqinfo object containing the chromosome names, lengths, and genome
  build.

- n:

  Integer describing the number of random sequences to generate

- ...:

  Additional arguments.

- interchromosomal:

  Boolean (TRUE/FALSE) indicating whether interchromosomal interactions
  should be allowed. Default is TRUE.

- .rows:

  (internal use only) vector of row positions to sample from seqinfo.

## Value

A GRanges or GInteractions object with ranges selected randomly with
replacement on the provided seqinfo.

## Examples

``` r
## Define Seqinfo containing chromosome info
if (require(TxDb.Hsapiens.UCSC.hg38.knownGene)) {
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    si <- seqinfo(txdb)
    si <- GenomeInfoDb::keepStandardChromosomes(si)
} else {
    si <- Seqinfo(
        seqnames=c("chr1", "chr2"),
        seqlengths=rep(200e6, 2),
        genome="hg38"
    )
}
#> Loading required package: TxDb.Hsapiens.UCSC.hg38.knownGene
#> Loading required package: GenomicFeatures
#> Loading required package: AnnotationDbi

## Make some GRanges
set.seed(123)
makeRandomGRanges(si, 100)
#> GRanges object with 100 ranges and 0 metadata columns:
#>         seqnames              ranges strand
#>            <Rle>           <IRanges>  <Rle>
#>     [1]    chr15   13828782-70429770      *
#>     [2]    chr15   53296438-78742890      *
#>     [3]    chr15   15454818-60964584      *
#>     [4]    chr15   21625558-57530611      *
#>     [5]    chr15    8277964-64393835      *
#>     ...      ...                 ...    ...
#>    [96]     chr2   80405238-95890355      *
#>    [97]     chr1 170893955-179388667      *
#>    [98]    chr16    5063321-35660766      *
#>    [99]    chr16   10589506-65167455      *
#>   [100]    chr16   24806742-37689512      *
#>   -------
#>   seqinfo: 25 sequences (1 circular) from hg38 genome

## Make some GInteractions
set.seed(123)
makeRandomGInteractions(si, n=100)
#> GInteractions object with 100 interactions and 0 metadata columns:
#>         seqnames1             ranges1     seqnames2            ranges2
#>             <Rle>           <IRanges>         <Rle>          <IRanges>
#>     [1]     chr15   13828782-70429770 ---      chrY    4147731-8428764
#>     [2]     chr15   53296438-78742890 ---      chrY  15005282-45321438
#>     [3]     chr15   15454818-60964584 ---      chr1  54781435-60189865
#>     [4]     chr15   21625558-57530611 ---      chr1 39402618-181970367
#>     [5]     chr15    8277964-64393835 ---      chr1 21648261-133668401
#>     ...       ...                 ... ...       ...                ...
#>    [96]      chr2   80405238-95890355 ---     chr20  13604328-59298515
#>    [97]      chr1 170893955-179388667 ---     chr20  35007959-62780398
#>    [98]     chr16    5063321-35660766 ---     chr20  13802693-47656725
#>    [99]     chr16   10589506-65167455 ---     chr19  35598384-44849872
#>   [100]     chr16   24806742-37689512 ---     chr19  31855464-51306138
#>   -------
#>   regions: 200 ranges and 0 metadata columns
#>   seqinfo: 25 sequences (1 circular) from hg38 genome

## Make some GInteractions only on same chromosome
set.seed(123)
makeRandomGInteractions(si, n=100, interchromosomal=FALSE)
#> GInteractions object with 100 interactions and 0 metadata columns:
#>         seqnames1             ranges1     seqnames2             ranges2
#>             <Rle>           <IRanges>         <Rle>           <IRanges>
#>     [1]     chr15   13828782-70429770 ---     chr15   14465720-14982368
#>     [2]     chr15   53296438-78742890 ---     chr15   51381722-77619820
#>     [3]     chr15   15454818-60964584 ---     chr15   29294834-45801465
#>     [4]     chr15   21625558-57530611 ---     chr15   53327816-96005291
#>     [5]     chr15    8277964-64393835 ---     chr15   60670359-68571354
#>     ...       ...                 ... ...       ...                 ...
#>    [96]      chr2   80405238-95890355 ---      chr2 118379329-165665755
#>    [97]      chr1 170893955-179388667 ---      chr1 112541676-185146454
#>    [98]     chr16    5063321-35660766 ---     chr16   46674131-69857380
#>    [99]     chr16   10589506-65167455 ---     chr16   35836834-39690234
#>   [100]     chr16   24806742-37689512 ---     chr16   11783461-46043479
#>   -------
#>   regions: 200 ranges and 0 metadata columns
#>   seqinfo: 25 sequences (1 circular) from hg38 genome

## Use specific binSizes
n <- 100
binOptions <- seq(5e3, 200e3, by=5e3)
si <- Seqinfo(seqnames="chr1", seqlengths=200e6, genome="hg38")
set.seed(123)
bins <- sample(binOptions, n, replace=TRUE)
makeRandomGInteractions(si, n) |>
    resize(bins) |>
    trim()
#> GInteractions object with 100 interactions and 0 metadata columns:
#>         seqnames1             ranges1     seqnames2             ranges2
#>             <Rle>           <IRanges>         <Rle>           <IRanges>
#>     [1]      chr1 121243273-121398272 ---      chr1   14211768-14366767
#>     [2]      chr1     7057416-7132415 ---      chr1   13550004-13625003
#>     [3]      chr1   41013509-41083508 ---      chr1     3618070-3688069
#>     [4]      chr1   70620055-70635054 ---      chr1   86204803-86219802
#>     [5]      chr1 102324459-102509458 ---      chr1   14961015-15146014
#>     ...       ...                 ... ...       ...                 ...
#>    [96]      chr1   32331617-32506616 ---      chr1   97867288-98042287
#>    [97]      chr1   44551938-44751937 ---      chr1     2859232-3059231
#>    [98]      chr1   37802441-37952440 ---      chr1   67576720-67726719
#>    [99]      chr1       643242-718241 ---      chr1 108206447-108281446
#>   [100]      chr1   68355278-68475277 ---      chr1 120060909-120180908
#>   -------
#>   regions: 400 ranges and 0 metadata columns
#>   seqinfo: 1 sequence from hg38 genome
```
