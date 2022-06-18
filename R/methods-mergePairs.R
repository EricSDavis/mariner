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

#' Find clusters by manhattan distance with DBSCAN
#' @param x `data.table` with at least 'start1' and
#'  'start2' columns. Data should be binned to
#'  the `binSize` parameter.
#' @param radius manhattan distance used to define
#'  a cluster of pairs.
#' @param binSize Integer (numeric) describing the
#'  resolution (range widths) of the paired data.
#'
#' @importFrom dbscan dbscan
#' @noRd
.findClusters <- function(x, radius, binSize) {
    d <- dist(x/binSize, method = 'manhattan')
    dbscan(d, eps = radius, minPts = 2)$cluster
}

#' Calculate new paired ranges with the mean of modes method
#' @param start1 vector of start1 positions
#' @param start2 vector of start2 positions
#' @param binSize Integer (numeric) describing the
#'  resolution (range widths) of the paired data.
#' @noRd
.newPairRanges <- function(start1, start2, binSize) {
    s1 <- .meanOfModes(start1, binSize)
    s2 <- .meanOfModes(start2, binSize)
    list(start1=s1, end1=s1+binSize, start2=s2, end2=s2+binSize)
}

#' Internal function for renaming mergePairs columns
#' @returns character vector of column names
#' @noRd
.renameCols <- function(x) {
    x |>
        gsub("^src_1$", "src", x = _) |>
        gsub("^id_1$", "id", x = _) |>
        gsub("^grp_1$", "grp", x = _) |>
        gsub("^clst_1$", "clst", x = _)
}

#' Internal mergePairs function
#' @inheritParams mergePairs
#' @importFrom rlang inform
#' @importFrom data.table as.data.table uniqueN `:=`
#' @noRd
.mergePairs <- function(x, binSize, radius, column) {

    ## Concatenate bedpe
    bedpe <- do.call(rbind, x)

    ## Check column name exists in bedpe
    if (!missing(column)) {
        if (!column %in% colnames(bedpe)) {
            abort(glue("Column '{column}' does not exist."))
        }
    }

    ## Convert to a GInteractions object for binning
    bedpe <- makeGInteractionsFromDataFrame(bedpe)

    ## Ensure interactions are binned
    binned <- .checkBinnedPairs(bedpe, binSize = binSize)
    if (!binned) {
        bedpe <- binPairs(x = bedpe,
                          binSize = binSize)
        msg <- c("Pairs are not binned to `binSize`.",
                 'i' = glue("Binning with binSize={binSize}, ",
                            "pos1='center' and pos2='center'."),
                 'i' = glue("Use `binPairs()` for different binning."))
        inform(msg)
    }

    ## Convert to data.table format
    dt <- as.data.table(bedpe)

    ## Annotate id & groups (essentially pairs of chromosomes)
    dt[,id := .I]
    dt[,grp := .GRP,by = .(seqnames1, seqnames2)]

    ## Initialize progress bar
    pb <- progress::progress_bar$new(
        format = "  :step [:bar] :percent elapsed: :elapsedfull",
        clear = FALSE,
        total = uniqueN(dt$grp)
    )
    pb$tick(0)

    ## Find clusters by manhattan distance with dbscan
    dt[,clst := {
        pb$tick(tokens=list(step="Merging pairs"));
        .findClusters(x = .SD[,c('start1','start2')],
                      radius = radius,
                      binSize = binSize)},
       by = grp]

    if (missing(column)) {
        ## Fast mean of modes
        selectionMethod <- "Mean of modes"

        ## Define start & end columns to update
        columnsToUpdate <- c("start1", "end1", "start2", "end2")

        ## Update with mean of modes for each group and cluster
        dt[clst != 0,
           (columnsToUpdate) :=
               .newPairRanges(start1 = .SD[['start1']],
                              start2 = .SD[['start2']],
                              binSize = binSize),
           by = .(grp, clst)]

        ## Select the first of the duplicated pairs
        single <- dt[clst == 0]
        selected <- dt[dt[clst != 0,.I[1], by = .(grp, clst)]$V1]
        mergedPairs <- rbind(single, selected)

    } else {
        selectionMethod <- glue("Select by column '{column}'")

        ## Edit column if it matches src, id, grp, or clst
        column <- gsub("(^src$|^id$|^grp$|^clst$)", "\\1_1", column)

        ## Fast find by column
        single <- dt[clst == 0]
        selected <- dt[dt[clst != 0,
                          .I[which.max(.SD[[column]])],
                          by=.(grp,clst)]$V1]
        mergedPairs <- rbind(single, selected)
    }

    ## Save and remove src, id, grp, and clst
    ids <- mergedPairs$id
    mergedPairs <-
        mergedPairs[,-c("src", "id", "grp", "clst")] |>
        as_ginteractions()

    ## Return to original column names
    colnames(mcols(mergedPairs)) <- .renameCols(colnames(mcols(mergedPairs)))

    ## Build MergedGInteractions object
    obj <-
        MergedGInteractions(delegate = mergedPairs,
                            ids = ids,
                            allPairs = dt,
                            selectionMethod = selectionMethod)

    return(obj)

}

#' Internal mergePairs for list type
#' @inheritParams mergePairs
#' @noRd
.mergePairsList <- function(x, binSize, radius, column) {

    ## Check that list is formatted correctly
    .checkListFormat(x)

    ## Read in Bedpe
    bedpe <- .readBedpeFromList(x)

    ## Rename columns to prevent overwriting
    newColumnNames <-
        lapply(bedpe, colnames) |>
        lapply(\(x) gsub("(^src$|^id$|^grp$|^clst$)", "\\1_1", x))
    bedpe <- Map(setNames, bedpe, newColumnNames)

    ## Add new column for source
    if (!is.null(names(bedpe))) {
        bedpe <- Map(cbind, bedpe, src = names(bedpe))
    } else {
        bedpe <- Map(cbind, bedpe, src = seq_along(bedpe))
    }

    ## Pass to internal merging function
    .mergePairs(x = bedpe,
                binSize = binSize,
                radius = radius,
                column = column)
}

#' Internal mergePairs for character type
#' @inheritParams mergePairs
#' @importFrom data.table fread
#' @noRd
.mergePairsCharacter <- function(x, binSize, radius, column) {

    ## Read in files as list of bedpe
    bedpe <- lapply(x, fread)

    ## Rename columns to prevent overwriting
    newColumnNames <-
        lapply(bedpe, colnames) |>
        lapply(\(x) gsub("(^src$|^id$|^grp$|^clst$)", "\\1_1", x))
    bedpe <- Map(setNames, bedpe, newColumnNames)

    ## Add new column for source
    bedpe <- Map(cbind, bedpe, src = basename(x))

    ## Pass to internal merging function
    .mergePairs(x = bedpe,
                binSize = binSize,
                radius = radius,
                column = column)
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
#' @param radius Character - distance measure
#'  passed to `dist()`. For available methods see
#'  `?dist()`.
#' @param column Integer or character denoting the
#'  column to be used to select among clustered
#'  interactions.
#'
#' @return Returns a merged `GInteractions` object.
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
#' x
#'
#' @rdname mergePairs
#' @export
setMethod("mergePairs",
          signature(x = 'list',
                    binSize = 'numeric_OR_missing',
                    radius = 'numeric_OR_missing',
                    column = 'character_OR_missing'),
          definition = .mergePairsList)

#' @rdname mergePairs
#' @export
setMethod("mergePairs",
          signature(x = 'character',
                    binSize = 'numeric_OR_missing',
                    radius = 'numeric_OR_missing',
                    column = 'character_OR_missing'),
          definition = .mergePairsCharacter)
