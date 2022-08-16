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
#' @noRd
setClassUnion("character_OR_function_OR_list",
              c("character", "function", "list"))

#' Class union for GInteractionsList
#' @importClassesFrom InteractionSet GInteractions
#' @noRd
setClassUnion("list_OR_GInteractions", c("list", "GInteractions"))


#' Virtual class for delegating GInteractions
#'
#' Uses a delegate `GInteractions` object during
#' initialization to assign its `GInteractions` slots.
#'
#' @slot delegate A `GInteractions` object used to initialize
#'  `GInteractions`-specific slots.
#' @slot anchor1 `anchorIds(delegate)$first`
#' @slot anchor2 `anchorIds(delegate)$second`
#' @slot regions `regions(delegate)`
#' @slot NAMES `names(delegate)`
#' @slot elementMetadata `elementMetadata(delegate)`
#' @slot metadata `metadata(delegate)`
#'
#' @seealso [InteractionSet::GInteractions]
#'
#' @rdname DelegatingGInteractions-class
#' @export
setClass(
    Class = "DelegatingGInteractions",
    contains = "VIRTUAL",
    slots = list(
        delegate = "GInteractions"
    )
)

#' Initialization method for DelegatingGInteractions
#' @importFrom InteractionSet anchorIds
#' @importFrom InteractionSet regions
#' @importFrom S4Vectors elementMetadata
#' @importFrom S4Vectors metadata
setMethod(
    f = "initialize",
    signature = "DelegatingGInteractions",
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

#' BinnedGInteractions Class
#'
#' The `BinnedGInteractions` class extends
#' `GInteractions` to enforce consistent
#' range widths.
#'
#' @slot delegate A `GInteractions` object used to initialize
#'  `GInteractions`-specific slots.
#' @slot firstBinSize Integer describing first pair bin width.
#' @slot secondBinSize Integer describing second pair bin width.
#' @slot pairBinsEqual Logical. TRUE if `firstBinSize` ==
#'  `secondBinSize` and FALSE if `firstBinSize` != `secondBinSize`.
#'
#' @seealso [InteractionSet::GInteractions]
#'
#' @rdname BinnedGInteractions-class
#' @export
setClass(
    Class = "BinnedGInteractions",
    contains = c("GInteractions", "DelegatingGInteractions"),
    slots = list(
        delegate = "GInteractions",
        firstBinSize = "integer",
        secondBinSize = "integer",
        pairBinsEqual = "logical"
    )
)

#' Constructor for BinnedGInteractions Class
#' @rdname BinnedGInteractions-class
#' @export
BinnedGInteractions <- function(delegate) {
    new("BinnedGInteractions", delegate = delegate)
}

#' Initialize BinnedGInteractions
#' @importFrom rlang abort
#' @importFrom glue glue
setMethod(
    f = "initialize",
    signature = "BinnedGInteractions",
    definition = function(.Object, ..., delegate = GInteractions()) {

        ## Get the widths of each pair
        w <- width(delegate)
        fbs <- w$first
        sbs <- w$second

        ## Check & set length of first bin
        if (length(unique(fbs)) == 1) {
            .Object@firstBinSize <- fbs[1]
        } else {
            msg <- c("First pair of ranges must all have equal widths.",
                     'i' = glue("Use `binPairs()` to bin paired ranges."))
            abort(msg)
        }

        ## Check & set length of second bin
        if (length(unique(sbs)) == 1) {
            .Object@secondBinSize <- sbs[1]
        } else {
            msg <- c("Second pair of ranges must all have equal widths.",
                     'i' = glue("Use `binPairs()` to bin paired ranges."))
            abort(msg)
        }

        ## Set pairBinsEqual flag
        if (.Object@firstBinSize == .Object@secondBinSize) {
            .Object@pairBinsEqual <- TRUE
        } else {
            .Object@pairBinsEqual <- FALSE
        }

        ## Pass GInteractions construction to virtual class
        .Object <- callNextMethod(.Object, ..., delegate = delegate)
        validObject(.Object)
        .Object
    })

#' Validate BinnedGInteractions
#' @name BinnedGInteractions-class
setValidity("BinnedGInteractions", function(object) {

    ## Get the widths of each pair
    w <- width(object)
    fbs <- w$first
    sbs <- w$second

    if (length(unique(fbs)) != 1)
        "Widths of the first anchor are not the same."

    if (length(unique(sbs)) != 1)
        "Widths of the second anchor are not the same."

    if (fbs[1] != object@firstBinSize)
        "Incorrect `firstBinSize`."

    if (sbs[1] != object@secondBinSize)
        "Incorrect `secondBinSize`."

    if (object@pairBinsEqual) {
        if (object@firstBinSize != object@secondBinSize) {
            "Incorrect `pairBinsEqual` flag."
        }
    }

    if (!object@pairBinsEqual) {
        if (object@firstBinSize == object@secondBinSize) {
            "Incorrect `pairBinsEqual` flag."
        }
    }

    return(TRUE)

})


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
    contains = c("GInteractions", "DelegatingGInteractions"),
    slots = list(
        delegate = "GInteractions",
        ids = "integer",
        allPairs = "data.table",
        selectionMethod = "character"
    )
)

#' Initialize MergedGInteractions
setMethod(
    f = "initialize",
    signature = "MergedGInteractions",
    definition = function(.Object, ..., delegate = GInteractions()) {

        ## Pass GInteractions construction to virtual class
        .Object <- callNextMethod(.Object, ..., delegate = delegate)
        .Object
    })

#' Method for parallel slots for MergedGInteractions
#' @importFrom S4Vectors parallel_slot_names
setMethod("parallel_slot_names", "MergedGInteractions", function(x) {
    c("ids", callNextMethod())
})
