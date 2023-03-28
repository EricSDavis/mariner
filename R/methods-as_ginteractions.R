#' Internal function for finding BEDPE format
#'
#' Can be either 6 or 10-column format.
#'
#' @inheritParams as_ginteractions
#'
#' @return A character string denoting BEDPE format.
#' @noRd
.findBedpeFormat <- function(df) {

    ## Use the first 10 columns of df
    df_cols <- colnames(df)[seq(1,10)]

    ## Define 10-column format
    cols <- c("seqnames1", "start1", "end1", "width1", "strand1",
              "seqnames2", "start2", "end2", "width2", "strand2")

    ## Check for 10 column format
    if (ncol(df) > 9) {

        if (identical(df_cols, cols) |
            identical(df_cols, gsub("seqnames", "chr", cols)) |
            identical(df_cols, gsub("seqnames", "chrom", cols))) {

            return("bedpe10")
        }
    }

    ## Otherwise, assume 6-column format
    return('bedpe6')
}

#' Internal function for splitting anchors
#'
#' Splits DataFrame-like object into anchors
#' and converts them to a GInteractions object.
#'
#' @param a1Coords vector of indices for the first anchor.
#' @param a2Coords vector of indices for the second anchor.
#' @inheritParams as_ginteractions
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom InteractionSet GInteractions
#' @importClassesFrom InteractionSet GInteractions
#'
#' @return A GInteractions object
#' @noRd
.splitAnchors <- function(df,
                          a1Coords=seq(1,3),
                          a2Coords=seq(4,6),
                          starts.in.df.are.0based) {

    ## Split into anchors
    a1 <- df[a1Coords] |> `colnames<-`(c('seqnames', 'start', 'end'))
    a2 <- df[a2Coords] |> `colnames<-`(c('seqnames', 'start', 'end'))

    ## Convert anchors to GRanges
    a1 <- makeGRangesFromDataFrame(df = a1,
                                   starts.in.df.are.0based =
                                       starts.in.df.are.0based)

    a2 <- makeGRangesFromDataFrame(df = a2,
                                   starts.in.df.are.0based =
                                       starts.in.df.are.0based)

    ## Return as GInteractions object
    return(GInteractions(a1, a2))

}

#' Internal function for `as_ginteractions` method
#' @inheritParams as_ginteractions
#'
#' @importFrom S4Vectors DataFrame
#' @importClassesFrom S4Vectors DataFrame
#' @importFrom glue glue
#' @importFrom rlang abort
#'
#' @noRd
.as_ginteractions <- function(df,
                              keep.extra.columns,
                              starts.in.df.are.0based) {

    ## Convert data.table/data.frame to DataFrame
    df <- DataFrame(df)

    ## Handle improper dimensions
    if (ncol(df) < 6) {
        msg <- c(glue("Improper dimensions in `df`."),
                 'i' = glue("There must be at least 6 columns."),
                 'x' = glue("You've supplied {ncol(df)} column(s)."))
        abort(msg)
    }

    ## Determine BEDPE 10 or 6-column format
    bedpeFormat <- .findBedpeFormat(df)

    if (identical(bedpeFormat, 'bedpe6')) {

        ## Split anchors
        gi <- .splitAnchors(df = df,
                            a1Coords = seq(1,3),
                            a2Coords = seq(4,6),
                            starts.in.df.are.0based =
                                starts.in.df.are.0based)

        ## Add in metadata columns
        if (keep.extra.columns & ncol(df) > 6) {
            mcols(gi) <- df[seq(7, ncol(df))]
        }
    }

    if (identical(bedpeFormat, 'bedpe10')) {

        ## Split anchors
        gi <- .splitAnchors(df = df,
                            a1Coords = seq(1,3),
                            a2Coords = seq(6,8),
                            starts.in.df.are.0based =
                                starts.in.df.are.0based)

        ## Add in metadata columns
        if (keep.extra.columns & ncol(df) > 10) {
            mcols(gi) <- df[seq(11, ncol(df))]
        }
    }

    ## Return GInteractions object
    return(gi)

}

#' @rdname as_ginteractions
#' @examples
#' ## data.frame
#' df <- data.frame(chr1 = "chr1", x1 = 10000, x2 = 20000,
#'                  chr2 = "chr1", y1 = 30000, y2 = 40000)
#' makeGInteractionsFromDataFrame(df)
#'
#' @export
setMethod("makeGInteractionsFromDataFrame",
          signature(df = 'DF_OR_df_OR_dt',
                    keep.extra.columns = 'logical_OR_missing',
                    starts.in.df.are.0based = 'logical_OR_missing'),
          definition = .as_ginteractions)

#' Convert DataFrames to GInteraction objects
#'
#' `as_ginteractions` takes
#' a paired-interaction (i.e. BEDPE) formatted
#' data-frame-like object and converts it to a
#' GInteractions object. For convenience,
#' `makeGInteractionsFromDataFrame` can be used as an alias.
#'
#' @param df A data.table, data.frame, or DataFrame object.
#'   Assumes that the first 6 colummns are in the format
#'   chr1, start1, end1 and chr2, start2, end2, representing
#'   each pair of interactions.
#' @param keep.extra.columns TRUE or FALSE (the default).
#'   If TRUE, the columns in df that are not used to form
#'   the genomic ranges of the returned GRanges object are
#'   then returned as metadata columns on the object.
#'   Otherwise, they are ignored. If df has a width column,
#'   then it's always ignored.
#' @param starts.in.df.are.0based TRUE or FALSE (the default).
#'   If TRUE, then the start positions of the genomic ranges
#'   in df are considered to be 0-based and are converted to
#'   1-based in the returned GRanges object. This feature is
#'   intended to make it more convenient to handle input that
#'   contains data obtained from resources using the
#'   "0-based start" convention. A notorious example of such
#'   resource is the UCSC Table Browser
#'   (http://genome.ucsc.edu/cgi-bin/hgTables).
#' @param ... Additional arguments.
#'
#' @return GInteraction object
#'
#' @rdname as_ginteractions
#' @aliases makeGInteractionsFromDataFrame
#'
#' @examples
#' ## data.frame
#' df <- data.frame(chr1 = "chr1", x1 = 10000, x2 = 20000,
#'                  chr2 = "chr1", y1 = 30000, y2 = 40000)
#' as_ginteractions(df)
#'
#' ## data.table
#' library(data.table)
#' df <- data.table::data.table(chr1 = "chr1", x1 = 10000, x2 = 20000,
#'                              chr2 = "chr1", y1 = 30000, y2 = 40000)
#' as_ginteractions(df)
#'
#' ## DataFrame
#' library(S4Vectors)
#' df <- DataFrame(chr1 = "chr1", x1 = 10000, x2 = 20000,
#'                 chr2 = "chr1", y1 = 30000, y2 = 40000)
#' as_ginteractions(df)
#'
#' ## Alias
#' df <- data.frame(chr1 = "chr1", x1 = 10000, x2 = 20000,
#'                  chr2 = "chr1", y1 = 30000, y2 = 40000,
#'                  pval = 0.05, dist = 10000)
#' makeGInteractionsFromDataFrame(df)
#'
#' ## Additional metadata
#' df <- data.frame(chr1 = "chr1", x1 = 10000, x2 = 20000,
#'                  chr2 = "chr1", y1 = 30000, y2 = 40000,
#'                  pval = 0.05, dist = 10000)
#' as_ginteractions(df)
#'
#' ## Remove additional metadata
#' as_ginteractions(df, keep.extra.columns = FALSE)
#'
#' ## Add 1 to starts (for 0-based programs)
#' as_ginteractions(df, starts.in.df.are.0based = TRUE)
#'
#' @importFrom S4Vectors DataFrame `mcols<-`
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom InteractionSet GInteractions
#' @importFrom data.table data.table
#' @importFrom rlang abort
#' @importFrom glue glue
#' @export
setMethod("as_ginteractions",
          signature(df = 'DF_OR_df_OR_dt',
                    keep.extra.columns = 'logical_OR_missing',
                    starts.in.df.are.0based = 'logical_OR_missing'),
          definition = .as_ginteractions)
