#' Check list type input for BEDPE and validate
#' @inheritParams mergePairs
#' @importFrom glue glue
#' @importFrom rlang abort
#' @return NULL (if no errors)
#' @noRd
.checkListFormat <- function(x) {

    ## Ensure the list is not empty
    if (length(x) == 0) {
        abort("List must have length > 0.")
    }

    ## Test for each type of list input
    if (is(x[[1]], 'data.frame')) {
        listType <- 'data.frame'
    } else if (is(x[[1]], 'DFrame')) {
        listType <- class(x[[1]])
    } else if (is(x[[1]], 'GInteractions')) {
        listType <- class(x[[1]])
    } else {
        listType <- unique(unlist(lapply(x, \(u) class(u))))
        msg <- c("Incorrect list input type.",
                 'i' = glue("Input must be data.frame-like ",
                            "object or a GInteractions object."),
                 'x' = glue("Your list is type ",
                            "{paste(listType, collapse=', ')}."))
        abort(msg)
    }

    ## Ensure all elements of the list
    ## are the same type.
    allTypesEqual <- all(unlist(lapply(x, \(u) is(u, listType))))

    ## If not equal, find the types and display
    ## them in an error message.
    if (!allTypesEqual) {
        types <- unique(unlist(lapply(x, \(u) class(u))))
        msg <- c("Incorrect list input type.",
                 'i' = "All objects in the list must be the same type.",
                 'x' = glue("Your list contains the ",
                            "following types: ",
                            "{paste(types, collapse=', ')}"))
        abort(msg)
    }
}

#' Read in Bedpe from list of objects
#'
#' Objects in the list can be `data.frame`-like
#' or GInteractions types.
#'
#' @inheritParams mergePairs
#' @importFrom data.table as.data.table
#' @return A list of data.tables
#' @noRd
.readBedpeFromList <- function(x) {

    ## Transform multiple times to make
    ## the structure consistent.
    lapply(x, \(u){
        u |>
            as.data.table() |>
            makeGInteractionsFromDataFrame() |>
            as.data.table()
    })
}

#' Internal mergePairs function
#' @inheritParams mergePairs
#' @noRd
.mergePairs <- function(x, binSize, column, distMethod, minPts) {

    return(x)

}

#' Internal mergePairs for list type
#' @inheritParams mergePairs
#' @noRd
.mergePairsList <- function(x, binSize, column,
                            distMethod, minPts) {

    ## Check that list is formatted correctly
    .checkListFormat(x)

    ## Read in Bedpe
    bedpe <- .readBedpeFromList(x)

    ## Add new column for source
    if (!is.null(names(bedpe))) {
        bedpe <- Map(cbind, bedpe, source = names(bedpe))
    } else {
        bedpe <- Map(cbind, bedpe, source = seq_along(bedpe))
    }

    ## Pass to internal merging function
    .mergePairs(x = bedpe,
                binSize = binSize,
                column = column,
                distMethod = distMethod,
                minPts = minPts)
}

#' Internal mergePairs for character type
#' @inheritParams mergePairs
#' @importFrom data.table fread
#' @noRd
.mergePairsCharacter <- function(x, binSize, column,
                                 distMethod, minPts) {

    ## Read in files as list of bedpe
    bedpe <- lapply(x, fread)

    ## Add new column for source
    bedpe <- Map(cbind, bedpe, source = basename(x))

    ## Pass to internal merging function
    .mergePairs(x = bedpe,
                binSize = binSize,
                column = column,
                distMethod = distMethod,
                minPts = minPts)
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
#' @param x Either a list of `GInteractions` objects,
#'  a list of `data.frame`-like objects, or a list of
#'  file paths of BEDPE-formatted data to be merged.
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
                    binSize = 'numeric_OR_missing',
                    column = 'character_OR_numeric',
                    distMethod = 'character_OR_missing',
                    minPts = 'numeric_OR_missing'),
          definition = .mergePairsList)

#' @rdname mergePairs
#' @export
setMethod("mergePairs",
          signature(x = 'character',
                    binSize = 'numeric_OR_missing',
                    column = 'character_OR_numeric',
                    distMethod = 'character_OR_missing',
                    minPts = 'numeric_OR_missing'),
          definition = .mergePairsCharacter)
