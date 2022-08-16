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

    ## Ensure every element of list
    ## is GInteractions
    if (is(x[[1]], 'GInteractions')) {
        listType <- class(x[[1]])
    } else {
        listType <- unique(unlist(lapply(x, \(u) class(u))))
        msg <- c("Incorrect list input type.",
                 'i' = glue("Input must be a list ",
                            "of GInteractions objects."),
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

#' Internal function for reading in lists of GInteractions
#' or GInteractions objects.
#' @inheritParams mergePairs
#' @importFrom data.table data.table as.data.table
#' @returns A concatenated data.table from the list
#' @noRd
.readGInteractionsList <- function(x, column) {
    ## Handle list or single GInteractions object
    if (is(x, "list")) {

        ## Check that list is formatted correctly
        .checkListFormat(x)

        ## Convert to list of data.tables
        bedpe <- lapply(x, as.data.table)

    } else {

        ## Convert to list of data.tables
        bedpe <- list(as.data.table(x))

    }

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

    ## Concatenate bedpe
    dt <- do.call(rbind, bedpe)

    ## Check column name exists in bedpe
    if (!missing(column)) {
        if (!column %in% colnames(dt)) {
            abort(glue("Column '{column}' does not exist."))
        }
    }

    ## Return concatenated data.table
    return(dt)
}

#' Find clusters by manhattan distance with DBSCAN
#' @param x `data.table` with at least 'start1' and
#'  'start2' columns.
#' @inheritParams mergePairs
#'
#' @importFrom dbscan dbscan
#' @noRd
.findClusters <- function(x, radius, method) {
    d <- dist(x, method = method)
    dbscan(d, eps = radius, minPts = 2)$cluster
}

#' Cluster pairs with DBSCAN
#'
#' @param x concatenated data.table
#' @inheritParams mergePairs
#' @returns Returns data.table with cluster information
#' @importFrom data.table as.data.table uniqueN `:=`
#' @noRd
.clusterPairs <- function(x, radius, method, pos) {

    ## Parse pos parameter
    pos <- match.arg(pos, choices = c("start", "end", "center"))

    ## Annotate id & groups (essentially pairs of chromosomes)
    x[,id := .I]
    x[,grp := .GRP,by = .(seqnames1, seqnames2)]

    ## Initialize progress bar
    pb <- progress::progress_bar$new(
        format = "  :step [:bar] :percent elapsed: :elapsedfull",
        clear = FALSE,
        total = uniqueN(x$grp)
    )
    pb$tick(0)

    ## Find clusters by distance with dbscan
    if (identical(pos, 'start')) {

        x[,clst := {
            pb$tick(tokens=list(step="Clustering pairs"));
            .findClusters(x = .SD[,c('start1','start2')],
                          radius = radius,
                          method = method)},
           by = grp]

    }

    if (identical(pos, 'end')) {

        x[,clst := {
            pb$tick(tokens=list(step="Clustering pairs"));
            .findClusters(x = .SD[,c('end1','end2')],
                          radius = radius,
                          method = method)},
           by = grp]
    }

    if (identical(pos, 'center')) {

        ## Calculate midpoint between starts & ends
        x$mid1 <- rowMeans(x[,c('start1', 'end1')])
        x$mid2 <- rowMeans(x[,c('start2', 'end2')])

        ## Cluster
        x[,clst := {
            pb$tick(tokens=list(step="Clustering pairs"));
            .findClusters(x = .SD[,c('mid1','mid2')],
                          radius = radius,
                          method = method)},
          by = grp]

        ## Remove calculated midpoints
        x[, "mid1" := NULL]
        x[, "mid2" := NULL]
    }

    ## Separate unique pairs by denoting as negative integers
    x[clst == 0, clst := seq(1, length(clst))*-1, by = grp]

    ## Return data.table with cluster information
    return(x)

}

#' Calculate new paired ranges with the mean of modes method
#' @param start1 vector of start1 positions
#' @param start2 vector of start2 positions
#' @param end1 vector of end1 positions
#' @param end2 vector of end2 positions
#'
#' Calculates the mean of modes. Results are rounded
#' to the nearest whole number.
#'
#' @returns list of updated start/end positions
#' @noRd
.newPairRanges <- function(start1, start2, end1, end2) {
    list(start1=round(mean(.modes(start1))),
         end1=round(mean(.modes(end1))),
         start2=round(mean(.modes(start2))),
         end2=round(mean(.modes(end2))))
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
.mergePairs <- function(x, radius, method, column, selectMax, pos) {

    ## Parse list or GInteractions input and return
    ## as a concatenated data.table
    dt <- .readGInteractionsList(x, column)

    ## Perform clustering
    dt <- .clusterPairs(x=dt, radius=radius, method=method, pos=pos)

    ## Perform selection
    if (missing(column)) {
        ## Fast mean of modes
        selectionMethod <- "Mean of modes"

        ## Define start & end columns to update
        columnsToUpdate <- c("start1", "end1", "start2", "end2")

        ## Update with mean of modes for each group and cluster
        dt[clst > 0,
           (columnsToUpdate) :=
               .newPairRanges(start1 = .SD[['start1']],
                              start2 = .SD[['start2']],
                              end1 = .SD[['end1']],
                              end2 = .SD[['end2']]),
           by = .(grp, clst)]

        ## Select the first of the duplicated pairs
        single <- dt[clst < 0]
        selected <- dt[dt[clst > 0,.I[1], by = .(grp, clst)]$V1]
        mergedPairs <- rbind(single, selected)

    } else {
        selectionMethod <- glue("Select by column '{column}'")

        ## Edit column if it matches src, id, grp, or clst
        column <- gsub("(^src$|^id$|^grp$|^clst$)", "\\1_1", column)

        ## Select max or min
        fun <- ifelse(selectMax, `which.max`, `which.min`)

        ## Fast find by column
        single <- dt[clst < 0]
        selected <- dt[dt[clst > 0,
                          .I[fun(.SD[[column]])],
                          by=.(grp,clst)]$V1]
        mergedPairs <- rbind(single, selected)

        ## Remove src
        mergedPairs[, "src" := NULL]
    }

    ## Save and remove id, grp, and clst (src removed above)
    ids <- mergedPairs$id
    mergedPairs <-
        mergedPairs[,c("id", "grp", "clst") := NULL] |>
        as_ginteractions()

    ## Return to original column names
    colnames(mcols(mergedPairs)) <- .renameCols(colnames(mcols(mergedPairs)))

    ## Remove metadata if using mean of modes
    ## since the choice is arbitrary.
    if (missing(column)) mcols(mergedPairs) <- NULL

    ## Build MergedGInteractions object
    obj <-
        MergedGInteractions(delegate = mergedPairs,
                            ids = ids,
                            allPairs = dt,
                            selectionMethod = selectionMethod)

    return(obj)

}


#' Merge sets of paired interactions
#'
#' Sets of paired range objects (like `GInteractions`
#' or BEDPE-formatted `data.frame`-like objects) are
#' first clustered by genomic distance with `dbscan`,
#' then a representative interaction is selected for
#' each cluster.
#'
#' Interactions are clustered into groups using the
#' provided base pair `radius`, and distance `method`
#' with `dbscan()`. Representative interactions are
#' selected for each group by one of two methods.
#' If `column` and `selectMax` arguments are provided,
#' the representative interaction with the maximum
#' (or minimum) value in `column` is returned for each
#' cluster. If these parameters are missing, new ranges
#' for each pair are returned by calculating the median
#' of modes for each cluster.
#'
#' @param x List of `GInteractions` or `data.frame`-like
#'  objects.
#' @param radius Numeric describing the distance in base
#'  pairs used to define a cluster or pairs.
#' @param method Character describing the distance measure
#'  to be used. This must be one of "euclidean", "maximum",
#'  "manhattan", "canberra", "binary" or "minkowski".
#'  Any unambiguous substring can be given.
#'  Default is "manhattan".
#' @param column Character denoting the column to be used
#'  to select among clustered interactions.
#' @param selectMax Logical. TRUE (default) uses `which.max()`
#'  to select the interaction pair. FALSE uses `which.min()`.
#'  Only applicable when `column` is specified.
#' @param pos Positions used for clustering pairs. Must be
#'  one of "start", "end" or "center". Default is "center".
#'
#' @return Returns a `MergedGInteractions` object.
#'
#' @examples
#' ## Reference BEDPE files (loops called with SIP)
#' bedpeFiles <-
#'     system.file("extdata", package = "mariner") |>
#'     list.files(pattern = "Loops.txt", full.names = TRUE)
#'
#' x <- mergePairs(x = bedpeFiles,
#'                 radius = 10e03,
#'                 column = "APScoreAvg")
#' x
#'
#' @rdname mergePairs
#' @export
setMethod("mergePairs",
          signature(x = 'list_OR_GInteractions',
                    radius = 'numeric',
                    method = 'character_OR_missing',
                    column = 'character_OR_missing',
                    selectMax = 'logical_OR_missing',
                    pos = "character_OR_missing"),
          definition = .mergePairs)
