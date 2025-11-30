# Visualize selection for a MatrixSelection object

Note: that buffer must be the same as the selection functions to work
appropriately

For \`selectCoordinates\`, \`rowInd\` and \`colInd\` are paired such
that the selected position in the matrix is \`c(rowInd\[1:i\],
colInd\[1:j\])\` for \`i\` rows and \`j\` columns.

## Usage

``` r
selectRadius(x, buffer, invert = FALSE)

selectCenterPixel(mhDist, buffer, invert = FALSE)

selectSubmatrix(m, invert = FALSE)

selectCoordinates(rowInd, colInd, buffer, invert = FALSE)

selectBlock(rowInd, colInd, buffer, invert = FALSE)

selectTopLeft(n, buffer, inset = 0, invert = FALSE)

selectTopRight(n, buffer, inset = 0, invert = FALSE)

selectBottomRight(n, buffer, inset = 0, invert = FALSE)

selectBottomLeft(n, buffer, inset = 0, invert = FALSE)

selectCorners(n, buffer, inset = 0, invert = FALSE)

selectRows(rows, buffer, invert = FALSE)

selectCols(cols, buffer, invert = FALSE)

selectInner(n, buffer, invert = FALSE)

selectOuter(n, buffer, invert = FALSE)

# S4 method for class 'MatrixSelection'
show(object)

# S4 method for class 'numeric'
selectRadius(x, buffer, invert = FALSE)

# S4 method for class 'numeric'
selectCenterPixel(mhDist, buffer, invert = FALSE)

# S4 method for class 'matrix'
selectSubmatrix(m, invert = FALSE)

# S4 method for class 'numeric'
selectCoordinates(rowInd, colInd, buffer, invert = FALSE)

# S4 method for class 'numeric'
selectBlock(rowInd, colInd, buffer, invert = FALSE)

# S4 method for class 'numeric'
selectTopLeft(n, buffer, inset = 0, invert = FALSE)

# S4 method for class 'numeric'
selectTopRight(n, buffer, inset = 0, invert = FALSE)

# S4 method for class 'numeric'
selectBottomRight(n, buffer, inset = 0, invert = FALSE)

# S4 method for class 'numeric'
selectBottomLeft(n, buffer, inset = 0, invert = FALSE)

# S4 method for class 'numeric'
selectCorners(n, buffer, inset = 0, invert = FALSE)

# S4 method for class 'numeric'
selectRows(rows, buffer, invert = FALSE)

# S4 method for class 'numeric'
selectCols(cols, buffer, invert = FALSE)

# S4 method for class 'numeric'
selectInner(n, buffer, invert = FALSE)

# S4 method for class 'numeric'
selectOuter(n, buffer, invert = FALSE)
```

## Arguments

- x:

  Integer vector of manhattan distances to select.

- buffer:

  Integer describing the number of pixels surrounding the central pixel.

- invert:

  Boolean indicating whether to invert the selection.

- mhDist:

  Integer vector of manhattan distances to select along with center
  pixel.

- m:

  matrix with 1's indicating selected positions and 0's indicated
  unselected positions.

- rowInd:

  Integer describing the row indices.

- colInd:

  Integer describing the column indices.

- n:

  Integer describing the number of outer pixels to select. Must be
  length of one.

- inset:

  Integer describing the number of pixels to inset the selection from
  the outer edge of the matrix. Default of 0 uses no inset.

- rows:

  Integer describing which rows to select.

- cols:

  Integer describing which cols to select.

- object:

  A MatrixSelection object.

## Value

A text-based visualization of the select matrix indices.

Numeric vector of matrix indices (byRow).

Numeric vector of matrix indices (byRow).

Numeric vector of matrix indices (byRow).

Numeric vector of matrix indices (byRow).

Numeric vector of matrix indices (byRow).

Numeric vector of matrix indices (byRow).

Numeric vector of matrix indices (byRow).

Numeric vector of matrix indices (byRow).

Numeric vector of matrix indices (byRow).

Numeric vector of matrix indices (byRow).

Numeric vector of matrix indices (byRow).

Numeric vector of matrix indices (byRow).

Numeric vector of matrix indices (byRow).

Numeric vector of matrix indices (byRow).

## Examples

``` r
res <- selectCenterPixel(0, 3)
show(res)
#> '0' = selected; '-' = unselected
#>                      
#>  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  - 
#>  -  -  -  0  -  -  - 
#>  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  - 
selectRadius(x=c(2,3,4), buffer=5, invert=FALSE)
#> '0' = selected; '-' = unselected
#>                                  
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  0  -  -  -  -  - 
#>  -  -  -  -  0  0  0  -  -  -  - 
#>  -  -  -  0  0  0  0  0  -  -  - 
#>  -  -  0  0  0  -  0  0  0  -  - 
#>  -  0  0  0  -  -  -  0  0  0  - 
#>  -  -  0  0  0  -  0  0  0  -  - 
#>  -  -  -  0  0  0  0  0  -  -  - 
#>  -  -  -  -  0  0  0  -  -  -  - 
#>  -  -  -  -  -  0  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
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
selectSubmatrix(m = matrix(rep(c(1,0,1), 3), nrow=3, ncol=3))
#> '0' = selected; '- ' = unselected
#>          
#>  0  0  0 
#>  -  -  - 
#>  0  0  0 
#> [1] 1 3 4 6 7 9
selectCoordinates(rowInd=1:3, colInd=1:3, buffer=5)
#> '0' = selected; '-' = unselected
#>                                  
#>  0  -  -  -  -  -  -  -  -  -  - 
#>  -  0  -  -  -  -  -  -  -  -  - 
#>  -  -  0  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
selectBlock(rowInd=1:3, colInd=1:3, buffer=5)
#> '0' = selected; '-' = unselected
#>                                  
#>  0  0  0  -  -  -  -  -  -  -  - 
#>  0  0  0  -  -  -  -  -  -  -  - 
#>  0  0  0  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
selectTopLeft(n=3, buffer=5, inset=1, invert=FALSE)
#> '0' = selected; '-' = unselected
#>                                  
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  0  0  0  -  -  -  -  -  -  - 
#>  -  0  0  0  -  -  -  -  -  -  - 
#>  -  0  0  0  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
selectTopRight(n=3, buffer=5, inset=1, invert=FALSE)
#> '0' = selected; '-' = unselected
#>                                  
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  0  0  0  - 
#>  -  -  -  -  -  -  -  0  0  0  - 
#>  -  -  -  -  -  -  -  0  0  0  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
selectBottomRight(n=3, buffer=5, inset=1, invert=FALSE)
#> '0' = selected; '-' = unselected
#>                                  
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  0  0  0  - 
#>  -  -  -  -  -  -  -  0  0  0  - 
#>  -  -  -  -  -  -  -  0  0  0  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
selectBottomLeft(n=3, buffer=5, inset=1, invert=FALSE)
#> '0' = selected; '-' = unselected
#>                                  
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  0  0  0  -  -  -  -  -  -  - 
#>  -  0  0  0  -  -  -  -  -  -  - 
#>  -  0  0  0  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
selectCorners(n=3, buffer=5, inset=1, invert=FALSE)
#> '0' = selected; '-' = unselected
#>                                  
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  0  0  0  -  -  -  0  0  0  - 
#>  -  0  0  0  -  -  -  0  0  0  - 
#>  -  0  0  0  -  -  -  0  0  0  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  0  0  0  -  -  -  0  0  0  - 
#>  -  0  0  0  -  -  -  0  0  0  - 
#>  -  0  0  0  -  -  -  0  0  0  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
selectRows(rows=1:3, buffer=5, invert=FALSE)
#> '0' = selected; '-' = unselected
#>                                  
#>  0  0  0  0  0  0  0  0  0  0  0 
#>  0  0  0  0  0  0  0  0  0  0  0 
#>  0  0  0  0  0  0  0  0  0  0  0 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
selectCols(cols=1:3, buffer=5, invert=FALSE)
#> '0' = selected; '-' = unselected
#>                                  
#>  0  0  0  -  -  -  -  -  -  -  - 
#>  0  0  0  -  -  -  -  -  -  -  - 
#>  0  0  0  -  -  -  -  -  -  -  - 
#>  0  0  0  -  -  -  -  -  -  -  - 
#>  0  0  0  -  -  -  -  -  -  -  - 
#>  0  0  0  -  -  -  -  -  -  -  - 
#>  0  0  0  -  -  -  -  -  -  -  - 
#>  0  0  0  -  -  -  -  -  -  -  - 
#>  0  0  0  -  -  -  -  -  -  -  - 
#>  0  0  0  -  -  -  -  -  -  -  - 
#>  0  0  0  -  -  -  -  -  -  -  - 
selectInner(n=1, buffer=5, invert=FALSE)
#> '0' = selected; '-' = unselected
#>                                  
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  0  0  0  -  -  -  - 
#>  -  -  -  -  0  0  0  -  -  -  - 
#>  -  -  -  -  0  0  0  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
#>  -  -  -  -  -  -  -  -  -  -  - 
selectOuter(n=1, buffer=5, invert=FALSE)
#> '0' = selected; '-' = unselected
#>                                  
#>  0  0  0  0  0  0  0  0  0  0  0 
#>  0  -  -  -  -  -  -  -  -  -  0 
#>  0  -  -  -  -  -  -  -  -  -  0 
#>  0  -  -  -  -  -  -  -  -  -  0 
#>  0  -  -  -  -  -  -  -  -  -  0 
#>  0  -  -  -  -  -  -  -  -  -  0 
#>  0  -  -  -  -  -  -  -  -  -  0 
#>  0  -  -  -  -  -  -  -  -  -  0 
#>  0  -  -  -  -  -  -  -  -  -  0 
#>  0  -  -  -  -  -  -  -  -  -  0 
#>  0  0  0  0  0  0  0  0  0  0  0 
```
