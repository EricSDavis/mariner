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
setClassUnion("character_OR_numeric", c("character", "numeric"))
#' @noRd
setClassUnion("character_OR_numeric_OR_missing",
              c("character", "numeric", "missing"))

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
#' @slot anchor1 `anchorIds(delegate)$first`
#' @slot anchor2 `anchorIds(delegate)$second`
#' @slot regions `regions(delegate)`
#' @slot NAMES `names(delegate)`
#' @slot elementMetadata `elementMetadata(delegate)`
#' @slot metadata `metadata(delegate)`
#'
#' @seealso [InteractionSet::GInteractions]
#'
#' @examples
#' ## Reference BEDPE files (loops called with SIP)
#' bedpeFiles <-
#'     system.file("extdata", package = "mariner") |>
#'     list.files(pattern = "Loops.txt", full.names = TRUE)
#'
#' x <- mergePairs(x = bedpeFiles,
#'                 binSize = 5e03,
#'                 radius = 2,
#'                 column = "APScoreAvg")
#' class(x)
#' @rdname MergedGInteractions-class
#' @export
MergedGInteractions <- setClass(
    Class = "MergedGInteractions",
    contains = "GInteractions",
    slots = list(
        delegate = "GInteractions",
        ids = "integer",
        allPairs = "data.table",
        selectionMethod = "character"
    )
)

#' Initialization method for MergedGInteractions
#' @importFrom InteractionSet anchorIds
#' @importFrom InteractionSet regions
#' @importFrom S4Vectors elementMetadata
#' @importFrom S4Vectors metadata
setMethod(
    f = "initialize",
    signature = "MergedGInteractions",
    definition = function(.Object, ..., delegate = GInteractions()) {

        ## Use GInteractions object to set other attributes
        .Object@delegate <- delegate
        .Object@anchor1 <- anchorIds(delegate)$first
        .Object@anchor2 <- anchorIds(delegate)$second
        .Object@regions <- regions(delegate)
        .Object@NAMES <- names(delegate)
        .Object@elementMetadata <- elementMetadata(delegate)
        .Object@metadata <- metadata(delegate)

        ## Use default class generator function
        .Object <- callNextMethod(.Object, ...)
        .Object
    })

#' Method for parallel slots for MergedGInteractions
#' @importFrom S4Vectors parallel_slot_names
setMethod("parallel_slot_names", "MergedGInteractions", function(x) {
    c("ids", callNextMethod())
})
