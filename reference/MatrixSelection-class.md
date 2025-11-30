# MatrixSelection Class

An object containing the selected indices of a matrix.

## Value

A MatrixSelection object (see description)

## Slots

- `x`:

  Vector of selected indices from a matrix of \`dim = buffer\*2+1\`.

- `buffer`:

  Integer indicating the buffer size, or number of pixels around a
  matrix.

## Examples

``` r
selectCenterPixel(0, 5)
#> '0' = selected; '-' = unselected
#>                                  
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  0  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
```
