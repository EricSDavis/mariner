
test_that("hdf5BlockApply works as expected", {

    ## Load utilities from DelayedArray
    library(DelayedArray)

    ## Create example array that is longer in the
    ## 3rd dimension (representing interactions)
    dims <- c(11L, 11L, 100L, 2L)
    a <- array(data=seq(1, prod(dims)), dim=dims)
    a <- DelayedArray(a)

    ## Define spacings, breaking up the longest dim
    ## Here we are processing in blocks of 10
    spacings <- dim(a)
    spacings[3] <- ceiling(spacings[3]/10)

    ## Define storage dimensions (all except those
    ## over which the function is being applied)
    storageDims <- dims[c(1,2,3)]

    ## Define chunk dimensions for writing to HDF5
    chunkDims <- storageDims
    chunkDims[3] <- spacings[3]

    ## Create grid for applying the data (grid)
    ## and grid for writing to the sink (sink_grid)
    grid <- RegularArrayGrid(dims, spacings)
    sink_grid <- RegularArrayGrid(storageDims, chunkDims)

    ## Create HDF5 file for writing
    h5 <- tempfile(fileext = ".h5")
    h5createFile(h5)

    ## Define compression for HDF5
    compressionLevel <- 0

    ## Create HDF5-backed realization sink
    sink <- HDF5RealizationSink(filepath=h5,
                                name="counts",
                                type="integer",
                                dim=storageDims,
                                chunkdim=chunkDims,
                                level=compressionLevel)

    ## Wrap function that operates on each block
    ## this can be anything, here it is sum
    FUN <- \(block) apply(block, c(1,2,3), sum)

    ## Read, apply, and write to HDF5
    ans <- hdf5BlockApply(x=a,
                          FUN=FUN,
                          sink=sink,
                          grid=grid,
                          sink_grid=sink_grid,
                          verbose=TRUE)

    expect_s4_class(ans, "HDF5Array")
    expect_identical(dim(ans), dims[c(1,2,3)])
    expect_identical(
        as.array(ans),
        apply(as.array(a), c(1,2,3), sum)
    )

})
