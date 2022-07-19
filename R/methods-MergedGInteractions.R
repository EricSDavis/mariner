## Accessors -------------------------------------------------------------------

#' Internal aggPairMcols function
#' @inheritParams aggPairMcols
#' @importFrom purrr pmap
#' @importFrom glue glue glue_collapse single_quote
#' @importFrom S4Vectors mcols
#' @importFrom rlang abort
#' @return `x` with aggregated metadata columns
#' @noRd
.aggPairMcols <- function(x, columns, funs){

    ## Pull out all pairs
    ap <- x@allPairs

    ## Check columns exist
    if (!any(columns %in% colnames(x@allPairs))) {
        abort(glue("Column(s) ",
                   "{glue_collapse(single_quote(columns), sep=',')} ",
                   "do(es) not exist."))
    }

    ## Create a key using ids from merged pairs
    keys <- ap[id %in% x@ids, .(id, grp, clst)]

    ## Join keys by group and cluster
    res <- ap[keys, on = .(grp, clst)]

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

    ## Apply functions to columns
    aggregatedColumns <-
        res[, pmap(.l = list(columns, funs),
                   .f = \(col, f) {
                       do.call(f, list(.SD[[col]]))
                   }), by = i.id][, i.id := NULL] |>
        `colnames<-`(value = c(colNames))

    if (nrow(aggregatedColumns) != length(x)) {
        msg <- c("Improper aggregation function.",
                 "i" = "`funs` must return a single value.")
        abort(msg)
    }

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
#' @examples
#' loopFiles <- list.files(path = system.file("extdata", package = "mariner"),
#'                         pattern = "Loops.txt",
#'                         full.names = TRUE)
#' x <- mergePairs(x = loopFiles)
#' sources(x)
#'
#' @rdname sources
#' @export
setMethod("sources", "MergedGInteractions",
          function(x) as.character(unique(x@allPairs$src)))

#' Create logial matrix of sources for each merged pair
#' @inheritParams subsetBySource
#' @importFrom data.table data.table
#' @return A logical matrix with columns for each source
#'  and rows for each merged pair indicating that pair's
#'  presence or abscence from the source (file or object).
#' @noRd
.makeSrcMatrix <- function(x) {

    ## Pull out all pairs
    ap <- x@allPairs

    ## Create a key using ids from merged pairs
    keys <- ap[id %in% x@ids, .(id, grp, clst)]

    ## Join keys by group and cluster
    res <- ap[keys, on = .(grp, clst)]

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
#' @return A `MergedGInteractions` object of de novo pairs.
#' @noRd
.subsetBySource <- function(x) {

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
#' ## Form GInteractions and convert to data.table
#' dt <- GInteractions(gr1, gr2) |> as.data.table()
#'
#' ## Split into two files
#' dts <- split(dt, c(rep(1,3), rep(2, 2)))
#'
#' ## Merge pairs
#' x <- mergePairs(dts, binSize = 10, radius = 2)
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



#' Internal accessor function for allPairs
#' @inheritParams allPairs
#' @importFrom S4Vectors isEmpty
#' @importFrom glue glue
#' @importFrom rlang abort
#' @noRd
.allPairs <- function(x) {

    ## Pull out all pairs
    ap <- x@allPairs

    ## Create a key using ids from merged pairs
    keys <- ap[id %in% x@ids, .(id, grp, clst)]

    ## Join keys by group and cluster
    res <- ap[keys, on = .(grp, clst)]

    ## Save the split vector & drop extra columns
    s <- res$i.id
    res <- res[, -c("id", "grp", "clst", "i.id")]

    ## Split and order result by input vector (slowest step)
    res <-
        split(res, s)[rank(x@ids)] |>
        `names<-`(values = names(x))

    return(res)
}

#' Get all pairs from MergedGInteractions object
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
#' ## Reference BEDPE files (loops called with SIP)
#' bedpeFiles <-
#'     system.file("extdata", package = "mariner") |>
#'     list.files(pattern = "Loops.txt", full.names = TRUE)
#'
#' x <- mergePairs(x = bedpeFiles,
#'                 binSize = 5e03,
#'                 radius = 2,
#'                 column = "APScoreAvg")
#' allPairs(x[1:3])
#' allPairs(x[3:1])
#' allPairs(x[c(3, 1, 2)])
#' allPairs(x)
#'
#' @rdname allPairs
#' @export
setMethod("allPairs", signature(x = "MergedGInteractions"),
          definition = .allPairs)

#' Get selectionMethod from MergedGInteractions object
#'
#' @param x MergedGInteractions object.
#' @param ... Additional arguments.
#'
#' @return A character vector describing which selection
#'  method was used for merging.
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
#' selectionMethod(x)
#'
#' @rdname selectionMethod
#' @export
setMethod("selectionMethod", "MergedGInteractions", function(x, ...) {
    x@selectionMethod
})
