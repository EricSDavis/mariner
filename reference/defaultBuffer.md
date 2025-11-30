# Return default buffer If InteractionArray is supplied, it uses the dimensions of counts matrices to set the buffer dimensions.

Return default buffer If InteractionArray is supplied, it uses the
dimensions of counts matrices to set the buffer dimensions.

## Usage

``` r
defaultBuffer(x)
```

## Arguments

- x:

  InteractionArray

## Value

5 (set default), the buffer of the provided InteractionArray, or an
error message if the InteractionArray is not odd and square (no buffer)
