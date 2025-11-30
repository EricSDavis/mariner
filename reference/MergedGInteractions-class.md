# MergedGInteractions Class

The \`MergedGInteractions\` class extends the \`GInteractions\` to
contain additional information about the pairs being merged.

## Value

A MergedGInteractions object (see description)

## Details

The \`MergedGInteractions\` class uses a delegate object during
initialization to assign its \`GInteractions\` slots. In addition to
containing information from all pairs, it also behaves as a
\`GInteractions\` object. \`mergePairs()\` builds this object.

## Slots

- `delegate`:

  A \`GInteractions\` object used to initialize
  \`GInteractions\`-specific slots. This is the mergedPairs set of
  interactions.

- `ids`:

  An integer vector of ids linking indices in the \`delegate\` slot all
  pairs (\`allPairs\` slot). These indices are parallel to \`delegate\`.

- `allPairs`:

  A \`data.table\` containing all input pairs combined. Also contains
  all metadata for each pair and 1) the source of the file, 2) an id, 3)
  which chromosome pair it belongs to (i.e. \`grp\`), and 4) the
  assigned cluster from \`dbscan\` (i.e. \`clst\`).

- `selectionMethod`:

  Character describing which method was used to select the final pair
  from the cluster of merged pairs.

## See also

\[InteractionSet::GInteractions\]

## Examples

``` r
## Load required packages
library(data.table, include.only="fread")

## Load marinerData
if (!require("marinerData", quietly = TRUE))
    BiocManager::install("marinerData")

## Reference BEDPE files (loops called with SIP)
bedpeFiles <- c(
    marinerData::FS_5kbLoops.txt(),
    marinerData::WT_5kbLoops.txt()
)
#> see ?marinerData and browseVignettes('marinerData') for documentation
#> downloading 1 resources
#> retrieving 1 resource
#> 
#> loading from cache
#> see ?marinerData and browseVignettes('marinerData') for documentation
#> loading from cache
names(bedpeFiles) <- c("FS", "WT")

## Read in bedpeFiles as a list of GInteractions
## Use only first 1000 rows for fast example
giList <-
    lapply(bedpeFiles, fread, nrows=1000) |>
    lapply(as_ginteractions)

## Cluster & merge pairs
x <- mergePairs(x = giList,
                radius = 10e03,
                column = "APScoreAvg")

class(x)
#> [1] "MergedGInteractions"
#> attr(,"package")
#> [1] "mariner"
```
