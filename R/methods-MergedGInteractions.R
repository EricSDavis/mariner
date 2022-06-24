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

#' Internal find de novo pairs
#' @inheritParams selectionMethod
#' @return A `MergedGInteractions` object of de novo pairs.
#' @noRd
.deNovo <- function(x) {

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

    ## Select merged pairs that are de novo
    ## There will be extra ids when using mean of modes.
    novo <- lapply(ids, \(y) x[x@ids %in% y])

    return(novo)
}

#' Find de novo pairs
#'
#' After merging sets of interaction pairs, `deNovo`
#' identifies which pairs are specific to each input set.
#' If several pairs are clustered together, they are
#' considered de novo if they all belong to the same
#' source set.
#'
#' @inheritParams selectionMethod
#'
#' @return A `MergedGInteractions` object of de novo pairs.
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
#' deNovo(x)
#'
#' @rdname deNovo
#' @export
setMethod("deNovo", signature(x = "MergedGInteractions"),
          definition = .deNovo)


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
