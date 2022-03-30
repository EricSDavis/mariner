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

## Create group
h5createGroup(file = fp, group = "group1")

## Adding data
h5write(obj = m,
        file = fp,
        name = "group1/matrix")

## Open file
fid <- H5Fopen(name = fp)
did <- H5Dopen(fid, "group1/matrix")

H5Fclose(fid)

## View structure of the files
h5ls(fp)
