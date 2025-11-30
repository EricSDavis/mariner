# InteractionArray Class

The \`InteractionArray\` class extends \`InteractionSet\` to provide an
interface for accessing submatrices pulled from Hi-C data.

## Usage

``` r
InteractionArray(assays, interactions, ...)

# S4 method for class 'ANY,GInteractions'
InteractionArray(assays, interactions, ...)

# S4 method for class 'missing,missing'
InteractionArray(assays, interactions, ...)

# S4 method for class 'InteractionArray'
show(object)

# S4 method for class 'InteractionArray'
rbind(..., deparse.level = 1)

# S4 method for class 'InteractionArray'
cbind(..., deparse.level = 1)
```

## Arguments

- assays, interactions:

  See
  `?`[`InteractionSet`](https://rdrr.io/pkg/InteractionSet/man/InteractionSet-class.html)

- ...:

  InteractionArray objects to be combined column-wise. All objects must
  be the same class.

- object:

  InteractionArray object.

- deparse.level:

  An integer scalar; see \`?base::cbind\` for a description of this
  argument.

## Value

An InteractionArray (see description)

## Details

This class is constructed with the \`pullHicMatrices()\` function when
all paired ranges have equal dimensions.

## See also

\[InteractionSet::InteractionSet\]

## Examples

``` r
InteractionArray()
#> class: InteractionArray 
#> dim: 0 interaction(s), 0 file(s), 0x0 count matrix(es)
#> metadata(0):
#> assays(0):
#> rownames: NULL
#> rowData names(0):
#> colnames: NULL
#> colData names(0):
#> type: GInteractions
#> regions: 0
```
