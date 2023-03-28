## Accessors -------------------------------------------------------------------

#' Internal function for mapping ids
#' for merged objects
#'
#' Map the ids in allPairs to the ids
#' used in the MergedGInteractions object.
#'
#' @param x MergedGInteractions object.
#' @return A key mapping ids by group and
#'  cluster. Contains a "i.id" column with
#'  the new id grouped by grp and clst.
#' @noRd
.mapIds <- function(x) {

    ## Suppress NSE notes in R CMD check
    id <- grp <- clst <- NULL

    ## Pull out all pairs
    ap <- x@allPairs

    ## Create a key using ids from merged pairs
    keys <- ap[id %in% x@ids, .(id, grp, clst)]

    ## Join keys by group and cluster
    res <- ap[keys, on = .(grp, clst)]

    ## Return mapped result
    return(res)
}

#' Internal aggPairMcols function
#' @inheritParams aggPairMcols
#' @importFrom purrr pmap
#' @importFrom glue glue glue_collapse single_quote
#' @importFrom S4Vectors mcols
#' @importFrom rlang abort
#' @rawNamespace import(data.table,
#'  except = c(between, shift, first, second, indices))
#' @return `x` with aggregated metadata columns
#' @noRd
.aggPairMcols <- function(x, columns, funs){

    ## Suppress NSE notes in R CMD check
    i.id <- NULL

    ## Check columns exist
    if (!any(columns %in% colnames(x@allPairs))) {
        abort(glue("Column(s) ",
                   "{glue_collapse(single_quote(columns), sep=',')} ",
                   "do(es) not exist."))
    }

    ## Map ids joined by group and cluster
    res <- .mapIds(x)

    ## Capture functions as a vector
    funs <- c(funs)

    ## Define column names
    colNames <-
        pmap(.l = list(columns, funs, seq_along(funs)),
             .f = \(col, f, i) {
                 if (is(f, "function"))
                     glue("fun{i}.{col}")
                 else
                     glue("{f}.{col}")
             }) |>
        unlist()

    ## Apply functions to columns (and retain ids)
    aggregated <-
        res[, pmap(.l = list(columns, funs),
                   .f = \(col, f) {
                       do.call(f, list(.SD[[col]]))
                   }), by = i.id]

    ## Check that aggregation was appropriate
    if (nrow(aggregated) != length(x)) {
        msg <- c("Improper aggregation function.",
                 "i" = "`funs` must return a single value.")
        abort(msg)
    }

    ## Add id's
    aggregatedWithIds <-
        cbind(aggregated, id = res[match(aggregated$i.id, res$i.id)]$id)

    ## Match to ids in x
    aggregatedColumns <-
        aggregatedWithIds[match(x@ids, aggregatedWithIds$id)]

    ## Remove excess columns
    aggregatedColumns[, c("i.id", "id") := NULL]

    ## Name columns
    colnames(aggregatedColumns) <- colNames

    ## Add result to x
    mcols(x) <- cbind(mcols(x), aggregatedColumns)

    return(x)
}

#' Aggregate the metadata columns of merged pairs
#'
#' @param x MergedGInteractions object.
#' @param columns Character vector of columns to
#'  aggregate.
#' @param funs Character vector of functions to apply to
#'  `columns`.
#'
#' @examples
#' ## Load required packages
#' library(data.table, include.only = "fread")
#' library(purrr, include.only = "pmap")
#'
#' ## Reference BEDPE files
#' bedpeFiles <-
#'     system.file("extdata", package = "mariner") |>
#'     list.files(pattern = "Loops.txt", full.names = TRUE)
#'
#' ## Read in bedpeFiles as a list of GInteractions
#' giList <-
#'     lapply(bedpeFiles, fread) |>
#'     lapply(as_ginteractions) |>
#'     setNames(gsub("^.*extdata/(.{2}).*$", "\\1", bedpeFiles))
#'
#' ## Add names describing the source and loop
#' giList <-
#'     pmap(.l = list(giList, names(giList)),
#'          .f = \(x, n) {
#'              x$name <- paste0(n, "_loop_", seq_len(length(x)))
#'              return(x)
#'          })
#'
#' ## Merge pairs
#' x <- mergePairs(x = giList,
#'                 radius = 5e03)
#'
#' ## List loop names
#' aggPairMcols(x, columns = "name", fun = "list")
#'
#' ## Aggregate values
#' aggPairMcols(x, columns = c("APScoreAvg"), fun = "mean")
#' aggPairMcols(x, columns = c("APScoreAvg", "avg"), fun = "mean")
#' aggPairMcols(x, columns = c("APScoreAvg"), fun = c("mean", "median"))
#'
#' ## Custom functions
#' aggPairMcols(x, columns = c("APScoreAvg"), fun = \(x) {
#'     ifelse(is.na(sd(x)), 0, sd(x))
#' })
#'
#' @return `x` with aggregated metadata columns
#'
#' @rdname aggPairMcols
#' @export
setMethod("aggPairMcols", signature(x = "MergedGInteractions",
                                    columns = "character",
                                    funs = "character_OR_function_OR_list"),
          definition = .aggPairMcols)

#' Accessor for sources
#'
#' Access the names or source files of
#' a `MergedGInteractions` object.
#'
#' @inheritParams subsetBySource
#' @return A character vector of names or source
#'  files of a `MergedGInteractions` object.
#' @examples
#' ## Load required packages
#' library(data.table, include.only="fread")
#'
#' loopFiles <-
#'     list.files(path = system.file("extdata", package = "mariner"),
#'                pattern = "Loops.txt",
#'                full.names = TRUE)
#'
#' giList <-
#'     lapply(loopFiles, fread) |>
#'     lapply(as_ginteractions)
#'
#' x <- mergePairs(x = giList,
#'                 radius = 10e03)
#'
#' sources(x)
#'
#' @rdname sources
#' @export
setMethod("sources", "MergedGInteractions",
          function(x) as.character(unique(x@allPairs$src)))

#' Create logial matrix of sources for each merged pair
#' @inheritParams subsetBySource
#' @importFrom data.table data.table as.data.table
#' @importFrom stats setNames
#' @return A logical matrix with columns for each source
#'  and rows for each merged pair indicating that pair's
#'  presence or abscence from the source (file or object).
#' @noRd
.makeSrcMatrix <- function(x) {

    ## Suppress NSE notes in R CMD check
    i.id <- NULL

    ## Map ids joined by group and cluster
    res <- .mapIds(x)

    ## Define the available sources
    srcs <- sources(x)

    ## Create source matrix
    srcMat <-
        lapply(srcs, \(x) x == res$src) |>
        setNames(srcs) |>
        as.data.table()

    ## Bind identifiers to the source matrix
    srcMat <- cbind(res, srcMat)

    ## Group by identifiers (x@ids)
    srcMat <- srcMat[,lapply(.SD, any), by=.(i.id), .SDcols=srcs]

    ## Put them in the correct order
    stopifnot(all(srcMat[rank(x@ids)]$i.id == x@ids))
    srcMat <- srcMat[rank(x@ids)]

    ## Drop id column and convert to pure matrix
    srcMat[, i.id := NULL]
    m <- as.matrix(srcMat)

    ## Return the source matrix
    return(m)
}

#' Check source names
#' @inheritParams subsetBySource
#' @param vec Input character vector (i.e. include, exclude or both)
#' @importFrom glue glue glue_collapse double_quote
#' @importFrom rlang abort
#' @return Nothing or Error message
#' @noRd
.checkSourceNames <- function(x, vec) {
    src <- sources(x)
    y <- !(vec %in% src)
    if (any(y)) {
        errorVals <-
            vec[which(y)] |>
            unique() |>
            double_quote() |>
            glue_collapse(sep = ', ')
        msg <- c(glue("{errorVals} not source option(s) ",
                      "for parameters `include` or `exclude`."),
                 "i" = glue("Use `sources(x)` to see available options."))
        abort(msg)
    }
}

#' Find subset pairs using include and exclude sources
#' @inheritParams subsetBySource
#' @importFrom data.table data.table
#' @return A `MergedGInteractions` object of subset pairs.
#' @noRd
.subsetBySourceIncludeExclude <- function(x, include, exclude) {

    ## Check source names
    .checkSourceNames(x, vec = c(include, exclude))

    ## Calculate the src matrix
    srcMat <- .makeSrcMatrix(x)

    ## Include and exclude
    inBool <- apply(srcMat[,include, drop=FALSE], 1, all)
    exBool <- apply(!srcMat[,exclude, drop=FALSE], 1, all)

    ## Use intersecting conditions to return mergedPairs
    if (any(inBool & exBool)) {
        x[inBool & exBool]
    } else {
        x[0]
    }
}

#' Find subset pairs using include sources
#' @inheritParams subsetBySource
#' @importFrom data.table data.table
#' @return A `MergedGInteractions` object of subset pairs.
#' @noRd
.subsetBySourceInclude <- function(x, include) {

    ## Check source names
    .checkSourceNames(x, vec = include)

    ## Calculate the src matrix
    srcMat <- .makeSrcMatrix(x)

    ## Include
    bool <- apply(srcMat[,include, drop=FALSE], 1, all)

    ## Use intersecting conditions
    if (any(bool)) {
        x[bool]
    } else {
        x[0]
    }
}

#' Find subset pairs using exclude sources
#' @inheritParams subsetBySource
#' @importFrom data.table data.table
#' @return A `MergedGInteractions` object of subset pairs.
#' @noRd
.subsetBySourceExclude <- function(x, exclude) {

    ## Check source names
    .checkSourceNames(x, vec = exclude)

    ## Calculate the src matrix
    srcMat <- .makeSrcMatrix(x)

    ## Exclude
    bool <- apply(!srcMat[,exclude, drop=FALSE], 1, all)

    ## Use intersecting conditions
    if (any(bool)) {
        x[bool]
    } else {
        x[0]
    }
}

#' Internal find de novo pairs
#' @inheritParams selectionMethod
#' @importFrom data.table data.table
#' @return A `MergedGInteractions` object of de novo pairs.
#' @noRd
.subsetBySource <- function(x) {

    ## Suppress NSE notes in R CMD check
    id <- grp <- clst <- src <- nDistinct <- NULL

    ## Pull out all pairs
    ap <- x@allPairs

    ## Find ids of de novo pairs
    ## All clst < 0 are new
    single <- ap[clst < 0, .(grp, clst, id, src, nDistinct = 1)]

    ## Only new if from the same source
    double <- ap[clst > 0,
                 .(id, src, nDistinct = length(unique(src))),
                 by = .(grp, clst)]

    ## Combine
    united <- rbind(single, double)

    ## Extract ids for each source
    ids <-
        split(united[nDistinct == 1], by = "src") |>
        lapply(`[[`, "id")

    ## Select merged pairs that are "de novo"
    ## There will be extra ids when using mean of modes.
    novo <- lapply(ids, \(y) x[x@ids %in% y])

    return(novo)
}

#' Subset MergedGInteractions by source
#'
#' Returns the subset of MergedGInteractions that belong
#' to each input source object (see these with `sources(x)`).
#' If the source pairs all come from the same object, their
#' corresponding merged pair is returned. However, if at least
#' one source pair comes from a different object, then that
#' merged pair is not returned.
#'
#' Optional `include` and `exclude` parameters modulate
#' the behaveior of `subsetBySource` to return different
#' subsets of originating pairs. For example, `include` requires
#' that the returned pairs be present in specific sources,
#' while `exclude` requires that returned pairs be absent
#' from specific sources. Sources not listed in either
#' `include` or `exclude` are ignored (they may or may not)
#' be present in the returned `MergedGInteractions` object.
#' `include` and `exclude` can be used indepedently or in
#' combination to return every possible set. If any of the
#' same sources are used in both `include` and `exclude` the
#' function will return a 0-length MergedGInteractions object.
#'
#' @inheritParams selectionMethod
#' @param include (Optional) A character vector of sources
#'  in which a pair must be present. For a list of available
#'  sources use `sources(x)`.
#' @param exclude (Optional) A character vector of sources
#'  in which a pair must be absent. For a list of available
#'  sources use `sources(x)`.
#'
#' @return A list of subsetted `MergedGInteractions` objects
#'  or a `MergedGInteractions` object (if `include` and/or
#'  `exclude` are used).
#'
#' @examples
#' ## Load required packages
#' library(GenomicRanges)
#' library(InteractionSet)
#'
#' ## Define example anchor regions
#' gr1 <-
#'     GRanges(seqnames = "chr1",
#'             ranges = IRanges(start = c(30,40,40,70,80),
#'                              end = c(40,50,50,80,90)))
#' gr2 <-
#'     GRanges(seqnames = "chr1",
#'             ranges = IRanges(start = c(30,30,50,10,30),
#'                              end = c(40,40,60,20,40)))
#'
#' ## Form GInteractions and split into two files
#' giList <- split(x = GInteractions(gr1, gr2),
#'                 f = c(rep(1,3), rep(2,2)))
#'
#' ## Merge pairs
#' x <- mergePairs(x = giList, radius = 20)
#'
#' subsetBySource(x)
#'
#' @rdname subsetBySource
#' @export
setMethod("subsetBySource",
          signature(x = "MergedGInteractions",
                    include = "missing",
                    exclude = "missing"),
          definition = .subsetBySource)

#' @rdname subsetBySource
#' @export
setMethod("subsetBySource",
          signature(x = "MergedGInteractions",
                    include = "character_OR_missing",
                    exclude = "missing"),
          definition = .subsetBySourceInclude)

#' @rdname subsetBySource
#' @export
setMethod("subsetBySource",
          signature(x = "MergedGInteractions",
                    include = "missing",
                    exclude = "character_OR_missing"),
          definition = .subsetBySourceExclude)

#' @rdname subsetBySource
#' @export
setMethod("subsetBySource",
          signature(x = "MergedGInteractions",
                    include = "character_OR_missing",
                    exclude = "character_OR_missing"),
          definition = .subsetBySourceIncludeExclude)



#' Internal accessor function for getPairClusters
#' @inheritParams getPairClusters
#' @importFrom data.table data.table
#' @noRd
.getPairClusters <- function(x) {

    ## Suppress NSE notes in R CMD check
    i.id <- NULL

    ## Map ids joined by group and cluster
    res <- .mapIds(x)

    ## Save the split vector & drop extra columns
    s <- res$i.id
    res <- res[, -c("id", "grp", "clst", "i.id")]

    ## Split and order result by input vector (slowest step)
    res <-
        split(res, s)[rank(x@ids)] |>
        `names<-`(value = names(x))

    return(res)
}

#' Get clustered pairs from MergedGInteractions object
#'
#' Returns the clustered pairs associated with each
#' range in the `MergedGInteractions` object. Order always
#' follows the indices of the `MergedGInteractions`
#' object.
#'
#' @inheritParams selectionMethod
#'
#' @return  A list of data.tables cooresponding to each pair
#'  in `x`.
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
#' getPairClusters(x[1:3])
#' getPairClusters(x[3:1])
#' getPairClusters(x[c(3, 1, 2)])
#' getPairClusters(x)
#'
#' @rdname getPairClusters
#' @export
setMethod("getPairClusters", signature(x = "MergedGInteractions"),
          definition = .getPairClusters)

#' Get selectionMethod from MergedGInteractions object
#'
#' @param x MergedGInteractions object.
#' @param ... Additional arguments.
#'
#' @return A character vector describing which selection
#'  method was used for merging.
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
#' selectionMethod(x)
#'
#' @rdname selectionMethod
#' @export
setMethod("selectionMethod", "MergedGInteractions", function(x, ...) {
    x@selectionMethod
})
