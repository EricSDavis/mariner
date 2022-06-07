#' Internal read in paired files
#' @importFrom data.table fread as.data.table
#' @importFrom glue glue
#' @importFrom rlang abort
#' @noRd
# .readFiles <- function(x)  {
#
#     ## Handle list of BEDPE files
#     if (is(x, "character")) {
#         return(lapply(x, fread))
#     }
#
#     ## Handle list of GInterctions objects
#     if (is(x, "list")) {
#
#         if (length(x) == 0) {
#             abort(glue("List must have length > 0."))
#         }
#
#         ## Are all list elements GInteractions objects?
#         GInteractionsList <-
#             lapply(x, \(y) is(y, "GInteractions")) |>
#             unlist() |>
#             all()
#
#         if (!GInteractionsList) {
#             msg <- c(glue("All list elements must be ",
#                           "GInteractions objects."))
#             abort(msg)
#         }
#
#         return(lapply(x, as.data.table))
#     }
# }

#' Internal mergePairs function
#' @inheritParams mergePairs
#' @noRd
.mergePairs <- function(x, binSize, column, distMethod, minPts) {

    ## Read in list of BEDPE
    bedpe <- .readFiles(x)

    ## Add new column for source
    bedpe <- Map(cbind, bedpe, source = basename(x))
}

#' Check list type input for BEDPE and validate
#' @inheritParams mergePairs
#' @importFrom glue glue
#' @importFrom rlang abort
#' @return Character; either "data.frame" or "GInteractions"
#' @noRd
.checkListType <- function(x) {

    ## Ensure the list is not empty
    if (length(x) == 0) {
        abort("List must have length > 0.")
    }

    ## Test for each type of list input
    if (is(x[[1]], 'data.frame')) {

        listType <- 'data.frame'
        arg <- listType # needed because data.table class is length 2

    } else if (is(x[[1]], 'DFrame')) {

        listType <- 'data.frame'
        arg <- class(x[[1]])

    } else if (is(x[[1]], 'GInteractions')) {

        listType <- "GInteractions"
        arg <- listType

    } else {

        listType <-
            lapply(x, \(u) class(u)) |>
            unlist() |>
            unique()
        msg <- c("Incorrect list input type.",
                 'i' = glue("Input must be data.frame-like ",
                            "object or a GInteractions object."),
                 'x' = glue("Your list is type ",
                            "{paste(listType, collapse=', ')}."))
        abort(msg)

    }

    ## Ensure all elements of the list
    ## are the same type.
    allTypesEqual <-
        lapply(x, \(u) is(u, arg)) |>
        unlist() |>
        all()

    ## If not equal, find the types and display
    ## them in an error message.
    if (!allTypesEqual) {
        types <-
            unique(unlist(lapply(x, \(u) class(u)))) |>
            unlist() |>
            unique()
        msg <- c("Incorrect list input type.",
                 'i' = "All objects in the list must be the same type.",
                 'x' = glue("Your list contains the ",
                            "following types: ",
                            "{paste(types, collapse=', ')}"))
        abort(msg)
    }

    return(listType)
}

#' Internal mergePairs for list type
#' @inheritParams mergePairs
#' @noRd
.mergePairsList <- function(x, binSize, column,
                            distMethod, minPts) {

    ## Determine the type of supplied list
    listType <- .checkListType(x)
}

#' Internal mergePairs for character type
#' @inheritParams mergePairs
#' @noRd
.mergePairsCharacter <- function(x, binSize, column,
                                 distMethod, minPts) {

    ## Determine the type of supplied list
    listType <- .checkListType(x)
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
                    binSize = 'numeric',
                    column = 'character_OR_numeric',
                    distMethod = 'character_OR_missing',
                    minPts = 'numeric'),
          definition = .mergePairsList)

#' @rdname mergePairs
#' @export
setMethod("mergePairs",
          signature(x = 'character',
                    binSize = 'numeric',
                    column = 'character_OR_numeric',
                    distMethod = 'character_OR_missing',
                    minPts = 'numeric'),
          definition = .mergePairsCharacter)
