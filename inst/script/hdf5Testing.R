library(rhdf5)

## Create a temporary file
fp <- tempfile(pattern = "mariner",
               fileext = '.h5')

## Create hdf5 file
h5createFile(fp)


## Generate sample data
set.seed(123)
s <- sample(1:100, 100, TRUE)
m <- matrix(data = s, nrow = 10, ncol = 10)

## Create group (for more structure)
# h5createGroup(file = fp, group = "apa")

## Adding data
h5write(obj = m,
        file = fp,
        name = "apa")

## Open file (needed for adding attributes)
# fid <- H5Fopen(name = fp)
# did <- H5Dopen(fid, "apa")
# H5Dclose(did)
# H5Fclose(fid)

## Read file
h5read(file = fp,
       name = 'apa',
       index = list(1:3, 1:3))


## View structure of the file
h5ls(fp)
h5delete(fp, 'a')
h5read(fp, "a")

## Create multidimensional array
h5createDataset(fp, 'a', c(10, 10, 3))

## Add data
h5write(m, file = fp, name = 'a', index = list(NULL, NULL, 3))

## Aggregate
apply(h5read(fp, 'a'), c(1,2), sum)

## Aggregate subset
apply(h5read(fp, 'a', index = list(1:3, 1:3, NULL)), c(1,2), sum)


## Generate large sample data ------------------------------------------------------------
library(rhdf5)

## Create a temporary file
fp <- tempfile(pattern = "mariner",
               fileext = '.h5')

## Create hdf5 file
h5createFile(fp)

## Create multidimensional array
n <- 2e6
h5createDataset(fp, 'a', c(10, 10, n))

## Generate sample array
set.seed(123)
s <- sample(1:100, 100*n, TRUE)
a <- array(data = s, dim = c(10, 10, n))
format(object.size(a), units = "auto")

## Write to hdf5 file
h5write(obj = a,
        file = fp,
        name = 'a',
        index = list(NULL, NULL, 1:n))

## Aggregate (from hdf5 file)
h5read(file = fp,
       name = 'a',
       index = list(NULL, NULL, 1:n)) |>
  apply(c(1,2), sum) |>
  system.time()

## Aggregate (on array directly)
apply(a[,,1:n], c(1,2), sum) |>
  system.time()

## Random access testing
set.seed(123)
r <- sample(x = 1:n, 1e6, FALSE)
apply(a[,, r], c(1,2), sum) |>
  system.time()
