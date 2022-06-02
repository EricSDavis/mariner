#' @rdname as_ginteractions
#' @export
setGeneric("as_ginteractions",
           function(df,
                    keep.extra.columns = TRUE,
                    starts.in.df.are.0based = FALSE,
                    ...)
               standardGeneric("as_ginteractions"))

#' @rdname as_ginteractions
#' @export
setGeneric("makeGInteractionsFromDataFrame",
           function(df,
                    keep.extra.columns = TRUE,
                    starts.in.df.are.0based = FALSE,
                    ...)
               standardGeneric("makeGInteractionsFromDataFrame"))
