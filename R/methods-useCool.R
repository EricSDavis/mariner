#' Function for checking if a file is a `.cool` or `.mcool` file, giving an
#' error 
#' 
#' @param fname path to file
#' 
#' @importFrom glue glue glue_collapse
#' @importFrom rlang abort
#' @importFrom rhdf5 H5Fis_hdf5 h5ls
#' 
#' @returns character string of file ending if file is `.cool` or `.mcool`, 
#' aborts with error message if not
#' 
.checkIfCool <- function(fname){
  
  ## Check if file is hdf5
  isHDF5 <- H5Fis_hdf5(fname)
  
  if(!isHDF5){
    abort(c("File must be a `.cool` or `.mcool` file",
            "x" = glue("{fname} is not in HDF5 format")))
  }
  
  ## Check if it is mcool or cool
  lvl1contents <- h5ls(fname, recursive = 1)
  
  expectedContents <- c("bins","chroms","indexes","pixels")
  
  if(nrow(lvl1contents) == 1){
    ## get deeper contents
    lvl2contents <- h5ls(fname, recursive = 3) |>
      pull(name) |>
      unique()
    
    ## check that the top level is resolution 
    ## and file contains expected datasets
    if(all(expectedContents %in% lvl2contents) & 
       lvl1contents$name == "resolutions"){
      return(".mcool")
    } 
  } else {
    ## check that file contains expected datasets
    if(all(expectedContents %in% lvl1contents$name)){
      return(".cool")
    } 
  }
  
  ## Abort with message if not a cool file
  expectedNames <- glue_collapse(glue("`{expectedContents}`"), sep = ", ")
  abort(c("File must be a `.cool` or `.mcool` file",
        "x" = glue("{fname} does not contain expected datasets"),
        "i" = glue("Expected datasets {expectedNames}")))
          
}

#' Function for reading available normalizations from .cool or
#' .mcool file 
#' 
#' @details
#' The "BALANCE" normalization refers to applying the pre-calculated matrix 
#' balancing weights in the `weight` dataset of `fname`, typically present
#' in files created using cooler. VC is vanilla coverage, 
#' VC_SQRT is square root of vanilla coverage, 
#' and KR is Knight-Ruiz normalization.
#' 
#' @param fname path to .cool or .mcool file
#' @param resolutions optional, specify which resolution(s) to read 
#' normalization types from. Default is all resolutions in `fname`
#' 
#' @importFrom rlang arg_match
#' @importFrom rhdf5 h5ls
#'
#' @returns Vector of available normalizations
#' 
#' @export
readCoolNormTypes <- function(fname, resolution){
  fileType <- .checkIfCool(fname)
  
  if (!resolution %in% readCoolBpResolutions(fname)) {
    abort(c(
      glue("resolution={resolution} is not valid."),
      "i"= glue("Use `readCoolBpResolutions()` to \\
                          see allowed values.")
    ))
  }
  
  datasetPath <- ""
  ## if `.mcool` file, add the proper resolution to path
  if(fileType == ".mcool"){
    datasetPath <- paste0("/resolutions/", as.integer(resolution))
  }
  
  ## get normalizations from the names of datasets in bins
  fileNorms <- h5read(fname, name = paste0(datasetPath,"/bins")) |> 
    names() 
  
  fileNorms <- fileNorms[!(fileNorms %in% c("chrom","start","end"))] |>
    gsub("weight","BALANCE", x = _) |> # replace "weight" with "BALANCE"
    c("NONE") # add "NONE" for raw counts
  
  return(fileNorms)
}

#' Function for reading chromosomes from a `.cool` or `.mcool` file
#'
#' @param fname path to .cool or .mcool file
#' @param resolution for .mcool files, 
#' 
#' @importFrom rhdf5 h5read
#'
#' @returns Data frame of chromosome names and lengths
#'
#' @export
readCoolChroms <- function(fname, resolution){
  ## check if input is a cool file
  fileType <- .checkIfCool(fname)
  
  ## setup correct dataset path
  datasetPath <- ""
  ## if `.mcool` file, add the proper resolution to path
  if(fileType == ".mcool"){
    validResolutions <- readCoolBpResolutions(fname)
    if(missing(resolution)){
      resolution = validResolutions[1]
    } else if (!resolution %in% validResolutions) {
      abort(c(
        glue("resolution={resolution} is not valid."),
        "i"= glue("Use `readCoolBpResolutions()` to \\
                          see allowed values.")
      ))
    }
    
    datasetPath <- paste0("/resolutions/",as.integer(resolution))
  }
  
  chromInfo <- h5read(fname, name = paste0(datasetPath,"/chroms")) |>
    as.data.frame()
  
  ## add index column of the order of chromosomes
  chromInfo$index <- 1:nrow(chromInfo)
  
  return(chromInfo)
}

#' Function for reading basepair resolutions from .cool or .mcool file
#'
#' @param fname path to .cool or .mcool file
#' @importFrom rhdf5 h5ls h5read
#'
#' @returns Vector of basepair resolutions
#' 
#' @export
readCoolBpResolutions <- function(fname){
  ## Check if file is `.mcool` or `.cool`
  fileType <- .checkIfCool(fname)
  
  ## find all valid resolutions
  if(fileType == ".mcool"){
    ## get contents from cool file
    contents <- h5ls(fname, recursive = 2)
    
    res <- contents[contents$group == "/resolutions","name"] |>
      as.integer() |>
      sort()
  } else {
    ## .cool files only have one resolution, calculate from bins
    start <- h5read(fname, name = "/bins/start", index = list(1))
    end <- h5read(fname, name = "/bins/end", index = list(1))
    res <- as.integer(end-start)
  }
  
  return(res)
}


#' Check cool arguments are valid against files
#'
#' @inheritParams pullHicMatrices
#' @importFrom rhdf5 h5read
#' @importFrom rlang abort arg_match
#' @importFrom glue glue
#' @noRd
.checkCoolArgs <- function(files, half, norm, binSize, matrix) {
  
  ## Check norm and binSize against each file
  for(f in files) {
    if (!binSize %in% readCoolBpResolutions(f)) {
      abort(c(
        glue("binSize={binSize} is not valid."),
        "i"= glue("Use `readCoolBpResolutions()` to \\
                          see allowed values.")
      ))
    }
    
    arg_match(norm, readCoolNormTypes(f, binSize))
  }
  
  ## Check half and matrix arguments
  arg_match(half, c("both", "upper", "lower"))
  
  if(matrix != "observed"){
    abort(c(glue("matrix=\"{matrix}\" is not valid."),
            "i"=glue("Only the \"observed\" matrix can be extracted from ", 
                     ".cool or .mcool files")))
  }
  return(NULL)
}

#' Check chromosomes in files
#'
#' Ensure they are identical for all files
#' and all chromosomes in `x` are contained
#' in these files.
#'
#' @inheritParams pullHicMatrices
#' @importFrom GenomeInfoDb seqinfo seqnames
#' @importFrom rlang abort
#' @importFrom glue glue
#' @returns Error if there is a chromosome issue
#' @noRd
.checkCoolChroms <- function(x, files, binSize) {
  
  ## Ensure all chromosome names and lengths
  chrInfo <- lapply(files, \(f) readCoolChroms(f, binSize))
  if (!all(vapply(chrInfo, identical, chrInfo[[1]],
                  FUN.VALUE = logical(1L)))) {
    abort(c(
      "Chromosomes in `files` are not identical",
      "i"="Check this with `readCoolChroms(files[1])`",
      "*"=glue("This is essential to ensure interactions ",
               "are ordered correctly."),
      "*"=glue("Call this function multiple times ",
               "or reprocess the Hi-C maps in the ",
               "same way to proceed.")
    ))
  }
  
  ## Extract chromosomes from x and files
  chromsInX <- seqnames(seqinfo(x))
  chromsInFile <- chrInfo[[1]]$name
  
  ## Ensure all chromosomes in x are in file
  if (!all(chromsInX %in% chromsInFile)) {
    abort(c(
      "There's some chr-chr-craziness going on.",
      'x' = glue("seqnames in `x` are not correctly ",
                 "formatted or do not exist in `files`."),
      "Try the following steps:",
      '>' = "Check `x` with `GenomeInfoDb::seqinfo(x)`.",
      '>' = "Check `file` with `readCoolChroms(files)`.",
      '>' = "Edit seqnames in `x` to match chromosomes in `files`.",
      "Hint:",
      '>' = glue("`GenomeInfoDb::seqlevelsStyle(x)",
                 " <- 'UCSC'` for 'chr' prefix."),
      '>' = glue("`GenomeInfoDb::seqlevelsStyle(x)",
                 " <- 'ENSEMBL'` without 'chr' prefix.")
    ))
  }
}


#' Equivalent to `strawr::straw` for `.cool` and `.mcool` files
#' 
#' @description
#' Reads the .hic file, finds the appropriate matrix and slice of data, 
#' and outputs as data.frame in sparse upper triangular format. 
#' Currently only supporting "observed" matrixes.
#' 
#' @param norm Normalization to apply. Must be one of the normalizations in the
#' given file. Use `readCoolNormTypes(fname)` for accepted normalization types
#' @param fname path to .cool or .mcool file
#' @param chr1loc first chromosome location, in the format "chr1:start1:end1"
#' @param chr2loc second chromosome location, in the format "chr2:start2:end2"
#' @param binsize The bin size in basepairs. Must be one of the resolutions in
#' the given file. Use `readCoolBpResolutions(fname)` for accepted binsizes
#'
#' @importFrom rlang abort
#' @importFrom glue glue
#' @importFrom rhdf5 h5ls h5read
#' @importFrom stringr str_extract
#' 
#' @returns data.frame of a sparse matrix of data from cool file. x,y,counts
#' 
#' @export
coolStraw <- function(norm, fname, chr1loc, chr2loc, binsize){
  ## parameter checking---------------------------------------------------------
  
  ## check that norm and binsize are in fname
  resolutions <- readCoolBpResolutions(fname)
  if(!binsize %in% resolutions){
    abort(glue("binsize={binsize} is not valid."),
          "i"= glue("Use `readCoolBpResolutions()` to \\
                          see allowed values."))
  }
  
  normalizations <- readCoolNormTypes(fname, binsize)
  if(!norm %in% normalizations){
    abort(glue("norm={norm} is not valid."),
          "i"= glue("Use `readCoolNormTypes()` to \\
                          see allowed values."))
  }
  
  ## check that chr1loc and chr2loc are in the correct format
  if(!(grepl("^\\w+\\:[0-9e]+\\:[0-9e]+$",chr1loc))) {
    abort(glue("chr1loc={chr1loc} is not valid."),
          i = "Chromosome locations must be in the format chr:start:end")
  }
  
  if(grepl(":",chr2loc)){
    if(!(grepl("^\\w+\\:[0-9e]+\\:[0-9e]+$",chr2loc))) {
      abort(glue("chr1loc={chr1loc} is not valid."),
            i = "Chromosome locations must be in the format chr:start:end")
    }
  }
  
  ## parse and adjust chromosome locations--------------------------------------
  
    ## parse chromosome locations
    chr1loc <- stringr::str_split(chr1loc,":")[[1]]
    chr1 <- chr1loc[1]
    start1 <- as.numeric(chr1loc[2])
    end1 <- as.numeric(chr1loc[3])
    
    chr2loc <- stringr::str_split(chr2loc,":")[[1]]
    chr2 <- chr2loc[1]
    start2 <- as.numeric(chr2loc[2])
    end2 <- as.numeric(chr2loc[3])
    
    ## Check that chr1 and chr2 are in the cool file
    chroms <- readCoolChroms(fname, binsize)
    if(!chr1 %in% chroms$name){
      abort(c(glue("Chromosome \"{chr1}\" not found in the file."),
            i = "Read chroms in file by using `readCoolChroms()`"))
    }
    if(!chr2 %in% chroms$name){
      abort(c(glue("Chromsome \"{chr2}\" not found in the file."),
            i = "Read chroms in file by using `readCoolChroms()`"))
    }
    
    ## get file type
    fileEnding <- .checkIfCool(fname)
    datasetPath <- ""
    ## if `.mcool` file, add the proper resolution to path
    if(fileEnding == ".mcool"){
      datasetPath <- paste0("/resolutions/", as.integer(binsize))
    }
    
    ## check that chr locations are in-bounds
    chrInfo <- readCoolChroms(fname, binsize)
    chr1length <- chrInfo[chrInfo$name==chr1,"length"]
    chr2length <- chrInfo[chrInfo$name==chr2,"length"]
    
    # set ends to chrom lengths if longer than chroms
    end1 <- ifelse(end1 <= chr1length, end1, chr1length)
    end2 <- ifelse(end2 <= chr2length, end2, chr2length)
    
    if(start1 < 0 | start2 < 0 | start1 > end1 | start2 > end2){
      rlang::abort(glue::glue("Chrom starts must be between 0 and ",
                              "chrom ends."))
    }
    
    ## Find bin ids for both chrom locations------------------------------------
    
    ## adjust start to start of bin it's in and end to start of bin it's in
    start1 <- (start1 %/% binsize) * binsize
    start2 <- (start2 %/% binsize) * binsize
    
    end1 <- ((end1 %/% binsize) * binsize)
    end2 <- ((end2 %/% binsize) * binsize)
    
    ## Get vector of chromosome offsets
    chrom_offsets <- h5read(fname, 
                            name=paste0(datasetPath,"/indexes/chrom_offset"))
    
    ## Find bin id for chr1 start and end
    chr1idx <- chrInfo[chrInfo$name==chr1,"index"]
    # get start1 bin
    chr1_starts <- h5read(fname, name=paste0(datasetPath,"/bins/start"),
                          index = list((chrom_offsets[chr1idx]+1):
                                         chrom_offsets[chr1idx+1]))
    start1bin <- which(chr1_starts == start1) + chrom_offsets[chr1idx] - 1
    # get end1 bin (i.e. bin where end1 starts)
    end1bin <- which(chr1_starts == end1) + chrom_offsets[chr1idx] - 1
    
    
    ## Find bin id for chr2 start and end
    chr2idx <- chrInfo[chrInfo$name==chr2,"index"]
    # get start2 bin
    chr2_starts <- h5read(fname, name=paste0(datasetPath,"/bins/start"),
                          index = list((chrom_offsets[chr2idx]+1):
                                         chrom_offsets[chr2idx+1]))
    start2bin <- which(chr2_starts == start2) + chrom_offsets[chr2idx] - 1
    # get end2 bin
    end2bin <- which(chr2_starts == end2) + chrom_offsets[chr2idx] - 1
    
    ## Extract counts and match locations---------------------------------------
    
    ## Get bin offsets for counts
    bin_offsets <- h5read(fname, name=paste0(datasetPath,"/indexes/bin1_offset"))
    
    ## Pull all bin2 ids for interactions with all bin1s between start1 and end1
    ## Read 5 million indices at a time for more efficient calling
    allBin1s <- (bin_offsets[start1bin+1]+1):(bin_offsets[end1bin+1]+1)
    
    binChunkSize <- 5e6
    
    if(length(allBin1s) > binChunkSize){
      binChunks <- seq(bin_offsets[start1bin+1]+1,
                       bin_offsets[end1bin+1]+1,
                       binChunkSize)
      
      count_idx <- numeric(0)
      for(binChunk in binChunks){
        binChunkEnd <- binChunk + binChunkSize - 1
        if(binChunkEnd > bin_offsets[end1bin+1]+1){
          binChunkEnd <- bin_offsets[end1bin+1]+1
        }
        if(binChunkEnd > binChunk){
        bin2s <- h5read(fname, name=paste0(datasetPath,"/pixels/bin2_id"),
                        index = list(binChunk:binChunkEnd)) 
        
        ## Find indexes for interactions with all bin2s between start2 and end2
        count_idx <- c(count_idx,
                       which(bin2s %in% start2bin:end2bin) + 
                         binChunk - 1)
        }
      } 
    }
    else {
      bin2s <- h5read(fname, name=paste0(datasetPath,"/pixels/bin2_id"),
                      index = list(allBin1s)) 
      
      ## Find indexes for interactions with all bin2s between start2 and end2
      count_idx <- which(bin2s %in% start2bin:end2bin) + 
        bin_offsets[start1bin+1]
    }
    
    ## Pull out bin1 ids, bin2 ids, and counts for all interactions in slice
    bin1ids <- h5read(fname, name=paste0(datasetPath,"/pixels/bin1_id"), 
                      index = list(count_idx))
    
    bin2ids <- h5read(fname, name=paste0(datasetPath,"/pixels/bin2_id"), 
                      index = list(count_idx))
    
    counts <- h5read(fname, name=paste0(datasetPath,"/pixels/count"), 
                     index = list(count_idx)) |>
      #defaults to integer, needs to be numeric to convert to NA correctly
      as.numeric() 
    
    ## Multiply by normalization factors, if applicable
    if(norm == "BALANCE"){
      bin1norm <- h5read(fname, name=paste0(datasetPath,"/bins/weight"),
                         list(as.numeric(bin1ids+1)))
      bin2norm <- h5read(fname, name=paste0(datasetPath,"/bins/weight"),
                         list(as.numeric(bin2ids+1)))
      counts <- counts / (bin1norm * bin2norm)
    } else if(norm != "NONE"){
      bin1norm <- h5read(fname, name=paste0(datasetPath,"/bins/",norm),
                         list(as.numeric(bin1ids+1)))
      bin2norm <- h5read(fname, name=paste0(datasetPath,"/bins/",norm),
                         list(as.numeric(bin2ids+1)))
      counts <- counts / (bin1norm * bin2norm)
    }
    
    ## Get genomic locations for each bin id 
    xs <- h5read(fname, name=paste0(datasetPath,"/bins/start"),
                index = list(as.numeric(bin1ids+1)))
    ys <- h5read(fname, name=paste0(datasetPath,"/bins/start"),
                 index = list(as.numeric(bin2ids+1)))
    
    ## create sparse matrix
    sparseMatrix <- data.frame(x = xs, y = ys, counts = counts)
    
    return(sparseMatrix)
}
