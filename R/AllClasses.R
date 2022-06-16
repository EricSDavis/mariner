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
