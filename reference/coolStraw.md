# Equivalent to \`strawr::straw\` for \`.cool\` and \`.mcool\` files

Reads the cooler file, finds the appropriate matrix and slice of data,
and outputs as data.frame in sparse upper triangular format. Currently
only supporting "observed" matrixes.

## Usage

``` r
coolStraw(norm, fname, chr1loc, chr2loc, binsize)
```

## Arguments

- norm:

  Normalization to apply. Must be one of the normalizations in the given
  file. Use \`readCoolNormTypes(fname)\` for accepted normalization
  types

- fname:

  path to .cool or .mcool file

- chr1loc:

  first chromosome location, in the format "chr1:start1:end1"

- chr2loc:

  second chromosome location, in the format "chr2:start2:end2"

- binsize:

  The bin size in basepairs. Must be one of the resolutions in the given
  file. Use \`readCoolBpResolutions(fname)\` for accepted binsizes

## Value

data.frame of a sparse matrix of data from cool file. x,y,counts

## Details

Normalizations are applied following the convention of the format they
come from. \`"BALANCE"\` reads cooler's own \`weight\` column, which
holds multiplicative biases, so counts are multiplied by the weight of
each bin. Any other normalization is read from a like-named bin column
and divided out instead, matching the divisive convention of juicer
vectors carried into cooler files by converters such as hic2cool, and
therefore matching \`strawr::straw()\`.
