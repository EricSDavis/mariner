test_that("Develop InteractionJaggedArray", {


    ## Develop method for storing, visualizing,
    ## subsetting, and coercing jagged (a.k.a
    ## irregular) arrays.

    ## Can't store in traditional array
    ## because they can't handle different
    ## dimensions (without making them all
    ## the same but having NA values, which
    ## would be pretty messy).

    ## Alternative approach is to record the
    ## dimensions of each matrix and store
    ## the data in separate groups for each
    ## .hic file.

    ## This will require an entirely different
    ## class (not inherited from InteractionSet)
    ## or from Summarized Experiment since these
    ## rely on strictly rectangular data structures.

    ## Goal will be to make it look as similar as
    ## possible to existing data structures for
    ## consistency in the interface.

    ## Create hdf5 dataset
    h5 <- tempfile(fileext=".h5")
    h5createFile(h5)

    ## Create sample data
    set.seed(123)
    nint = 10000 # n interactions (i.e. loops)
    dim1 = sample(1:20, nint, TRUE) # row dimensions
    dim2 = sample(1:20, nint, TRUE) # col dimensions
    data = list( # data to fill for each hic file
        rep(seq_len(nint), dim1*dim2),
        rep(seq_len(nint), dim1*dim2),
        rep(seq_len(nint), dim1*dim2),
        rep(seq_len(nint), dim1*dim2),
        rep(seq_len(nint), dim1*dim2),
        rep(seq_len(nint), dim1*dim2),
        rep(seq_len(nint), dim1*dim2),
        rep(seq_len(nint), dim1*dim2),
        rep(seq_len(nint), dim1*dim2),
        rep(seq_len(nint), dim1*dim2),
        rep(seq_len(nint), dim1*dim2)
    )
    names(data) <- as.character(seq_along(data)) # character for hdf5 storage
    offset2 = cumsum(dim1*dim2) # offset in data for each matrix
    offset1 = c(0L, offset2[-length(offset2)])+1
    slices <- data.table(dim1, dim2, offset1, offset2) # use this as DelayedArray
    slices <- as.matrix(slices)

    ## Write data to hdf5
    h5write(slices, h5, "slices")
    h5write(data, h5, "data")


    ## Example for accessing a subset of these data
    i = 70:73
    files = c(1, 10)
    if (is(files, "numeric")) { #fix this logic to be flexible
        files <- paste0("data/", files)
    }
    slices <- h5read(h5, 'slices', index=list(i, seq_len(4)))
    res <-
        lapply(files, \(k) {
            lapply(seq_len(nrow(slices)), \(j){
                slice = seq(from=slices[j,3], to=slices[j,4])
                data = h5read(h5, k, index=list(slice))
                matrix(data=data, nrow=slices[j,1], ncol=slices[j,2])
            })
        })

    if (length(res) == 1) {
        res <- res[[1]]
    }
    str(res)

    res
    slices


})
