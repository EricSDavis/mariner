# Function for reading available normalizations from .cool or .mcool file

Function for reading available normalizations from .cool or .mcool file

## Usage

``` r
readCoolNormTypes(fname, resolution)
```

## Arguments

- fname:

  path to .cool or .mcool file

- resolution:

  optional, specify which resolution(s) to read normalization types
  from. Default is all resolutions in \`fname\`

## Value

Vector of available normalizations

## Details

The "BALANCE" normalization refers to applying the pre-calculated matrix
balancing weights in the \`weight\` dataset of \`fname\`, typically
present in files created using cooler. VC is vanilla coverage, VC_SQRT
is square root of vanilla coverage, and KR is Knight-Ruiz normalization.
