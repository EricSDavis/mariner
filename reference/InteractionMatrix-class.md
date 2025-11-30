# InteractionMatrix Class

The \`InteractionMatrix\` class extends the \`InteractionSet\` to
provide an interface for accessing the count matrix pulled from Hi-C
data.

## Usage

``` r
InteractionMatrix(assays, interactions, ...)

# S4 method for class 'ANY,GInteractions'
InteractionMatrix(assays, interactions, ...)

# S4 method for class 'missing,missing'
InteractionMatrix(assays, interactions, ...)

# S4 method for class 'InteractionMatrix'
show(object)

# S4 method for class 'InteractionMatrix'
rbind(..., deparse.level = 1)

# S4 method for class 'InteractionMatrix'
cbind(..., deparse.level = 1)
```

## Arguments

- assays, interactions:

  See
  `?`[`InteractionSet`](https://rdrr.io/pkg/InteractionSet/man/InteractionSet-class.html)

- ...:

  InteractionMatrix objects to be combined column-wise. All objects must
  be the same class.

- object:

  InteractionMatrix object.

- deparse.level:

  An integer scalar; see \`?base::cbind\` for a description of this
  argument.

## Value

An InteractionMatrix (see description)

## Details

This class is constructed with the \`pullHicPixels()\` function when all
paired ranges define a single pixel.

## See also

\[InteractionSet::InteractionSet\]

## Examples

``` r
InteractionMatrix()
#> class: InteractionMatrix 
#> dim: count matrix with 0 interactions and 0 file(s)
#> metadata(0):
#> assays(0):
#> rownames: NULL
#> rowData names(0):
#> colnames: NULL
#> colData names(0):
#> type: GInteractions
#> regions: 0
```
