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
#' @importClassesFrom S4Vectors SimpleList List
#' @returns A concatenated data.table from the list
#' @noRd
.readGInteractionsList <- function(x, column) {
    ## Handle list or single GInteractions object
    if (is(x, "list")) {

        ## Check that list is formatted correctly
        .checkListFormat(x)

        ## Convert to list of data.tables
        bedpe <- lapply(x, as.data.table)

    } else if (is(x, "SimpleList")) {

        ## Convert to "list"
        x <- as.list(x)

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
    if (!is.null(column)) {
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
#' @importFrom dbscan dbscan
#' @importFrom stats dist
#' @return vector of DBSCAN cluster designation.
#' @noRd
.findClusters <- function(x, radius, method) {
    d <- dist(x, method = method)
    dbscan(d, eps = radius, minPts = 2)$cluster
}

#' Cluster pairs with DBSCAN
#' @param x concatenated data.table
#' @inheritParams mergePairs
#' @importFrom data.table as.data.table uniqueN `:=`
#' @importFrom progress progress_bar
#' @returns Returns data.table with cluster information
#' @noRd
.clusterPairs <- function(x, radius, method, pos) {

    ## Suppress NSE notes in R CMD check
    id <- grp <- clst <- seqnames1 <- seqnames2 <- NULL

    ## Parse pos parameter
    pos <- match.arg(pos, choices = c("start", "end", "center"))

    ## Annotate id & groups (essentially pairs of chromosomes)
    x[,id := .I]
    x[,grp := .GRP,by = .(seqnames1, seqnames2)]

    ## Rename "radius" and "method" columns if they exist
    colnames(x) <- gsub("(^radius$|^method$)", "\\1_1", colnames(x))

    ## Initialize progress bar
    pb <- progress_bar$new(
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

        ## Rename "mid1" and "mid2" columns if they exist
        colnames(x) <- gsub("(^mid1$|^mid2$)", "\\1_1", colnames(x))

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

        ## Restore column names
        colnames(x) <- gsub("(^(mid1)_1$|^(mid2)_1$)", "\\1", colnames(x))
    }

    ## Restore column names
    colnames(x) <- gsub("(^(radius)_1$|^(method)_1$)", "\\1", colnames(x))

    ## Separate unique pairs by denoting as negative integers
    x[clst == 0, clst := seq(1, length(clst))*-1, by = grp]

    ## Return data.table with cluster information
    return(x)

}

#' Calculate new paired ranges with the mean of modes method
#'
#' Calculates the mean of modes. Results are rounded
#' to the nearest whole number.
#'
#' @param start1 vector of start1 positions
#' @param start2 vector of start2 positions
#' @param end1 vector of end1 positions
#' @param end2 vector of end2 positions
#'
#' @returns list of updated start/end positions.
#' @noRd
.newPairRanges <- function(start1, start2, end1, end2) {
    list(start1=round(mean(.modes(start1))),
         end1=round(mean(.modes(end1))),
         start2=round(mean(.modes(start2))),
         end2=round(mean(.modes(end2))))
}

#' Internal function for renaming mergePairs columns
#' @param x vector of column names
#' @returns character vector of column names
#' @noRd
.renameCols <- function(x) {
    x |>
        gsub("^src_1$", "src", x = _) |>
        gsub("^id_1$", "id", x = _) |>
        gsub("^grp_1$", "grp", x = _) |>
        gsub("^clst_1$", "clst", x = _) |>
        gsub("^mid1_1$", "mid1", x = _) |>
        gsub("^mid2_1$", "mid2", x = _) |>
        gsub("^radius_1$", "radius", x = _) |>
        gsub("^method_1$", "method", x = _)
}


#' Internal mergePairs function
#' @inheritParams mergePairs
#' @importFrom rlang inform
#' @importFrom data.table as.data.table uniqueN `:=` copy
#' @importFrom S4Vectors mcols
#' @noRd
.mergePairs <- function(x, radius, method, column, selectMax, pos) {

    ## Suppress NSE notes in R CMD check
    clst <- grp <- NULL

    ## Check argument types
    ## column is checked later (can be string or NULL)
    .checkTypes(c(
        method="string",
        selectMax="boolean",
        pos="string"
    ))

    ## Parse list or GInteractions input and return
    ## as a concatenated data.table
    dt <- .readGInteractionsList(x, column)

    ## Perform clustering
    dt <- .clusterPairs(x=dt, radius=radius, method=method, pos=pos)

    ## Perform selection
    if (is.null(column)) {
        ## Fast mean of modes
        selectionMethod <- "Mean of modes"

        ## Define start & end columns to update
        columnsToUpdate <- c("start1", "end1", "start2", "end2")

        ## Make copy of pairs to retain original ranges
        dt2 <- copy(dt)

        ## Update with mean of modes for each group and cluster
        dt2[clst > 0,
           (columnsToUpdate) :=
               .newPairRanges(start1 = .SD[['start1']],
                              start2 = .SD[['start2']],
                              end1 = .SD[['end1']],
                              end2 = .SD[['end2']]),
           by = .(grp, clst)]

        ## Select the first of the duplicated pairs
        single <- dt2[clst < 0]
        selected <- dt2[dt2[clst > 0,.I[1], by = .(grp, clst)]$V1]
        mergedPairs <- rbind(single, selected)

    } else {
        selectionMethod <- glue("Select by column '{column}'")

        ## Edit column if it matches src, id, grp, or clst
        pattern <- "(^src$|^id$|^grp$|^clst$|^radius$|^method$)"
        column <- gsub(pattern, "\\1_1", column)

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
    if (is.null(column)) mcols(mergedPairs) <- NULL

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
#' Sets of paired range objects (i.e., `GInteractions`)
#' are first clustered by genomic distance with `dbscan`,
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
#' ## Load marinerData
#' if (!require("marinerData", quietly = TRUE))
#'     install.packages("marinerData")
#'
#' bedpeFiles <- c(
#'     marinerData::FS_5kbLoops.txt(),
#'     marinerData::WT_5kbLoops.txt()
#' )
#' names(bedpeFiles) <- c("FS", "WT")
#'
#' ## Read in bedpeFiles as a list of GInteractions
#' ## Use only first 1000 rows for fast example
#' giList <-
#'     lapply(bedpeFiles, read.table, header=TRUE, nrows=1000) |>
#'     lapply(as_ginteractions)
#'
#' ## Cluster & merge pairs
#' x <- mergePairs(x = giList,
#'                 radius = 10e03,
#'                 column = "APScoreAvg")
#' x
#'
#' @rdname mergePairs
#' @export
setMethod("mergePairs",
          signature(x = 'list_OR_SimpleList_OR_GInteractions',
                    radius = 'numeric'),
          definition = .mergePairs)
