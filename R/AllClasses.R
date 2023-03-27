## Generic class unions --------------------------------------------------------

#' Class union for "DataFrame-like" objects
#' @importClassesFrom data.table data.table
#' @importClassesFrom S4Vectors DFrame
#' @import methods
#' @noRd
setClassUnion("DF_OR_df_OR_dt", c("DFrame", "data.frame", "data.table"))

#' Class unions for general types
#' @noRd
setClassUnion("character_OR_missing", c("character", "missing"))
#' @noRd
setClassUnion("logical_OR_missing", c("logical", "missing"))
#' @noRd
setClassUnion("numeric_OR_missing", c("numeric", "missing"))
#' @noRd
setClassUnion("integer_OR_NULL", c("integer", "NULL"))
#' @noRd
setClassUnion("character_OR_numeric", c("character", "numeric"))
#' @noRd
setClassUnion("character_OR_numeric_OR_missing",
              c("character", "numeric", "missing"))
#' @noRd
setClassUnion("character_OR_function_OR_list",
              c("character", "function", "list"))

#' Class union for DelayedMatrix and matrix
#' @importClassesFrom DelayedArray DelayedMatrix
#' @noRd
setClassUnion("DelayedMatrix_OR_matrix",
              c("DelayedMatrix", "matrix"))

#' Class union for GInteractionsList
#' @importClassesFrom S4Vectors SimpleList
#' @importClassesFrom InteractionSet GInteractions
#' @noRd
setClassUnion("list_OR_SimpleList_OR_GInteractions",
              c("list", "SimpleList", "GInteractions"))

#' Class union for GInteractions or InteractionSet
#' @importClassesFrom InteractionSet InteractionSet GInteractions
#' @noRd
setClassUnion("GInteractions_OR_InteractionSet",
              c("GInteractions", "InteractionSet"))


## DelegatingGInteractions class -----------------------------------------------

#' Virtual class for delegating GInteractions
#'
#' Uses a delegate `GInteractions` object during
#' initialization to assign its `GInteractions` slots.
#'
#' @slot delegate A `GInteractions` object used to initialize
#'  `GInteractions`-specific slots.
#'
#' @seealso [InteractionSet::GInteractions]
#'
#' @rdname DelegatingGInteractions-class
#' @return DelegatingGInteractions virtual class
#' @export
setClass(
    Class = "DelegatingGInteractions",
    contains = "VIRTUAL",
    slots = list(
        delegate = "GInteractions"
    )
)

#' Internal method for updating GInteractions objects
#' @param .Object GInteractions object
#' @param ... Additional arguments
#' @param delegate GInteractins object to use for
#'  updating slots of `.Object`.
#' @return An updated GInteractions object.
#' @noRd
.updateGInteractions <- function(.Object, ..., delegate = GInteractions()) {

    ## Use GInteractions object to set other attributes
    .Object@anchor1 <- anchorIds(delegate)$first
    .Object@anchor2 <- anchorIds(delegate)$second
    .Object@regions <- regions(delegate)
    .Object@NAMES <- names(delegate)
    .Object@elementMetadata <- elementMetadata(delegate)
    .Object@metadata <- metadata(delegate)

    .Object
}

#' Initialization method for DelegatingGInteractions
#' @param .Object GInteractions object
#' @param ... Additional arguments
#' @param delegate GInteractions object
#' @importFrom InteractionSet anchorIds
#' @importFrom InteractionSet regions
#' @importFrom S4Vectors elementMetadata
#' @importFrom S4Vectors metadata
#' @return DelegatingGInteractions object
setMethod(
    f = "initialize",
    signature = "DelegatingGInteractions",
    definition = function(.Object, ..., delegate = GInteractions()) {

        ## Set delegate (for children of virtual class)
        .Object@delegate <- delegate

        ## Update GInteractions
        .Object <- .updateGInteractions(.Object = .Object, ...,
                                        delegate = delegate)

        ## Use default class generator function
        .Object <- callNextMethod(.Object, ...)
        .Object
    })

## MergedGInteractions Class ---------------------------------------------------

#' MergedGInteractions Class
#'
#' The `MergedGInteractions` class extends the
#' `GInteractions` to contain additional information
#' about the pairs being merged.
#'
#' The `MergedGInteractions` class uses a delegate object
#' during initialization to assign its `GInteractions` slots.
#' In addition to containing information from all pairs, it
#' also behaves as a `GInteractions` object. `mergePairs()`
#' builds this object.
#'
#' @slot delegate A `GInteractions` object used to initialize
#'  `GInteractions`-specific slots. This is the mergedPairs
#'  set of interactions.
#' @slot ids An integer vector of ids linking indices in the
#'  `delegate` slot all pairs (`allPairs` slot). These indices
#'  are parallel to `delegate`.
#' @slot allPairs A `data.table` containing all input pairs
#'  combined. Also contains all metadata for each pair and
#'  1) the source of the file, 2) an id, 3) which chromosome
#'  pair it belongs to (i.e. `grp`), and 4) the assigned
#'  cluster from `dbscan` (i.e. `clst`).
#' @slot selectionMethod Character describing which method
#'  was used to select the final pair from the cluster of
#'  merged pairs.
#'
#' @seealso [InteractionSet::GInteractions]
#'
#' @examples
#' ## Load required packages
#' library(data.table, include.only="fread")
#'
#' ## Reference BEDPE files (loops called with SIP)
#' bedpeFiles <-
#'     system.file("extdata", package = "mariner") |>
#'     list.files(pattern = "Loops.txt", full.names = TRUE)
#'
#' giList <-
#'     lapply(bedpeFiles, fread) |>
#'     lapply(as_ginteractions)
#'
#' x <- mergePairs(x = giList,
#'                 radius = 10e03,
#'                 column = "APScoreAvg")
#'
#' class(x)
#' @rdname MergedGInteractions-class
#' @export
MergedGInteractions <- setClass(
    Class = "MergedGInteractions",
    contains = c("GInteractions", "DelegatingGInteractions"),
    slots = list(
        delegate = "GInteractions",
        ids = "integer",
        allPairs = "data.table",
        selectionMethod = "character"
    )
)

#' Initialize MergedGInteractions
#' @param .Object GInteractions object
#' @param ... Additional arguments
#' @param delegate GInteractions object
#' @return MergedGInteractions
setMethod(
    f = "initialize",
    signature = "MergedGInteractions",
    definition = function(.Object, ..., delegate = GInteractions()) {

        ## Pass GInteractions construction to virtual class
        .Object <- callNextMethod(.Object, ..., delegate = delegate)
        .Object
    })

#' Method for parallel slots for MergedGInteractions
#' @param x See
#'   \code{? \link[S4Vectors]{parallel_slot_names}}
#' @importFrom S4Vectors parallel_slot_names
#' @return MergedGInteractions with parallel slots
setMethod("parallel_slot_names", "MergedGInteractions", function(x) {
    c("ids", callNextMethod())
})

## InteractionArray Class ------------------------------------------------------

#' InteractionArray Class
#'
#' The `InteractionArray` class extends
#' `InteractionSet` to provide an interface
#' for accessing submatrices pulled from
#' Hi-C data.
#'
#' This class is constructed with the
#' `pullHicMatrices()` function when all
#' paired ranges have equal dimensions.
#'
#' @returns An InteractionArray (see description)
#'
#' @seealso [InteractionSet::InteractionSet]
#'
#' @examples
#' InteractionArray()
#'
#' @rdname InteractionArray-class
#' @export
setClass(
    Class = "InteractionArray",
    contains = "InteractionSet"
)

## InteractionMatrix Class -----------------------------------------------------

#' InteractionMatrix Class
#'
#' The `InteractionMatrix` class extends the
#' `InteractionSet` to provide an interface
#' for accessing the count matrix pulled from
#' Hi-C data.
#'
#' This class is constructed with the
#' `pullHicPixels()` function when all
#' paired ranges define a single pixel.
#'
#' @returns An InteractionMatrix (see description)
#'
#' @seealso [InteractionSet::InteractionSet]
#'
#' @examples
#' InteractionMatrix()
#'
#' @rdname InteractionMatrix-class
#' @export
setClass(
    Class = "InteractionMatrix",
    contains = "InteractionSet"
)

#' Set validity for InteractionMatrix class
#' @name InteractionMatrix-class
setValidity("InteractionMatrix", function(object) {

    ## If there are assays...
    if (length(object@assays) == 0) {
        return(TRUE)
    }

    ## And they are HDF5...
    if (!is(object@assays@data@listData$counts, "HDF5Array")) {
        return(TRUE)
    }

    ## Check for .h5 file
    h5File <- object@assays@data@listData$counts@seed@filepath
    if (!file.exists(h5File)) {
        return(".h5 file doesn't exist. Update path to .h5 file.")
    }

    return(TRUE)

})

## MatrixSelection Class -------------------------------------------------------

#' MatrixSelection Class
#'
#' An object containing the selected indices
#' of a matrix.
#'
#' @slot x Vector of selected indices from
#'  a matrix of `dim = buffer*2+1`.
#' @slot buffer Integer indicating the
#'  buffer size, or number of pixels
#'  around a matrix.
#'
#' @examples
#' selectCenterPixel(0, 5)
#'
#' @rdname MatrixSelection-class
#' @export
MatrixSelection <- setClass(
    Class = "MatrixSelection",
    slots = list(
        x = "numeric",
        buffer = "numeric"
    )
)

