#' Internal mergePairs function
#' @inheritParams mergePairs
#' @noRd
.mergePairs <- function(x, binSize, column, distMethod, minPts) {
    return("hello!")
}

#' Merge sets of paired interactions
#'
#' Sets of paired range objects (like `GInteractions`
#' or BEDPE-formatted `data.frame`-like objects) are
#' merged by genomic distance with dbscan.
#'
#' Interactions are clustered into groups using the
#' `binSize`, `dist_method`, and `minPts` parameters
#' passed to `dbscan()`. A column and function provided
#' by the user is then used to select which interaction
#' should be chosen to represent the group. The resulting
#' object reports the source (which file or list
#' name/index) of the chosen interaction.
#'
#' @param x Either a list of `GInteractions` objects
#'  or a list of file paths of BEDPE-formatted data
#'  to be merged.
#' @param binSize Integer (numeric) describing the
#'  resolution (range widths) of the paired data.
#'  Used to determine the epsilon value for `dbscan()`
#'  (i.e. `eps = binSize*2`).
#' @param column Integer or character denoting the
#'  column to be used to select among clustered
#'  interactions.
#' @param distMethod Character - distance measure
#'  passed to `dist()`. For available methods see
#'  `?dist()`.
#' @param minPts Numeric (integer) of the minimum
#'  number of interactions to form a `dbscan` cluster.
#'
#' @return Returns a merged `GInteractions` object.
#'
#' @rdname mergePairs
#' @export
setMethod("mergePairs",
          signature(x = 'list',
                    binSize = 'numeric',
                    column = 'character_OR_numeric',
                    distMethod = 'character_OR_missing',
                    minPts = 'numeric'),
          definition = .mergePairs)

#' @rdname mergePairs
#' @export
setMethod("mergePairs",
          signature(x = 'character',
                    binSize = 'numeric',
                    column = 'character_OR_numeric',
                    distMethod = 'character_OR_missing',
                    minPts = 'numeric'),
          definition = .mergePairs)
