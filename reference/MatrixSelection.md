# Names S3 method for autocomplete

Names S3 method for autocomplete

Extract \`\$\` operator for MatrixSelection

Extract \`\[\[\` operator for MatrixSelection

Concatenate MatrixSelection objects

Concatenate MatrixSelection objects

## Usage

``` r
# S3 method for class 'MatrixSelection'
names(x)

# S4 method for class 'MatrixSelection'
x$name

# S4 method for class 'MatrixSelection,ANY,ANY'
x[[i]]

# S4 method for class 'MatrixSelection,MatrixSelection'
e1 + e2

# S4 method for class 'MatrixSelection,MatrixSelection'
e1 - e2
```

## Arguments

- x:

  A \`MatrixSelection\` object.

- name:

  Name of slot.

- i:

  A character or numeric to extract.

## Value

different data or metadata concerning a MatrixSelection object
