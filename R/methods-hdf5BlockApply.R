#' Internal HDF5-backed blockApply
#' @inheritParams hdf5BlockApply
#' @importFrom DelayedArray read_block write_block DelayedArray close
#' @returns DelayedArray object
#' @noRd
.hdf5BlockApply <- function(x, FUN, sink, grid, sink_grid, verbose) {
    nblocks <- length(grid)
    invisible({
        lapply(seq_along(grid), \(bid) {
            ## Read block
            if (verbose)
                message("/ Reading and realizing block ",
                        bid, "/", nblocks, " ... ",
                        appendLF=FALSE)
            x_viewport <- grid[[bid]]
            block <- read_block(x, x_viewport)
            if (verbose) message("OK")

            ## Apply function
            if (verbose)
                message("\\ Processing it ... ", appendLF=FALSE)
            block_ans <- FUN(block)

            ## Write to HDF5-backed sink
            sink_viewport <- sink_grid[[bid]]
            sink <- write_block(sink, sink_viewport, block_ans)
            if (verbose) message("OK")
        })
    })

    ## Close and coerce
    DelayedArray::close(sink)
    return(as(sink, "DelayedArray"))
}

#' HDF5-backed blockApply
#'
#' Read in array data in blocks, apply function,
#' and write back to an HDF5 file.
#'
#' @param x Delayed Array object.
#' @param FUN Function that takes one argument
#'  'block' and processes it.
#' @param sink HDF5RealizationSink object.
#' @param grid ArrayGrid over array `x`.
#' @param sink_grid ArrayGrid over `sink`.
#' @param verbose Logical - whether block processing
#'  progress should be displayed.
#' @returns DelayedArray object.
#' @rdname hdf5BlockApply
#' @export
setMethod("hdf5BlockApply",
          signature(x="DelayedArray"),
          definition = .hdf5BlockApply)
