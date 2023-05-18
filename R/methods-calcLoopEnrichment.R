#' Internal for modifyEnrichFun
#' @param FUN custom function to modify
#' @param eframe environment with local variables to use in the function
#' @importFrom rlang parse_expr abort inform
#' @importFrom glue glue
#' @return function with parameters x, fg, and bg. All other parameters
#'  are converted to either default values or local variable values
#' @noRd
.modifyEnrichFun <- function(FUN, eframe=parent.frame(2)) {
  ## check that FUN is a function
  if(!is.function(FUN)){
    abort(c("`FUN` must be a function.",
            "x"=glue("`class(FUN)` is {class(FUN)}.")))
  }
  
  ## get argument names of function
  argNames <- formalArgs(FUN) 
  
  ## check if any arguments are function names
  areFunNames <- vapply(argNames, function(x) {
    exists(x) && is.function(get(x))
  }, logical(1))
  
  if(sum(areFunNames)>0){
    funNames <- glue_collapse(glue("`{argNames[areFunNames]}`"), sep = ", ")
    abort(c("`FUN` arguments cannot be the same as any function name.",
            "x"=glue("{funNames} share name(s) with a function.")))
  }
  
  ## rename arguments to avoid conflict
  x <- reconcileArgs("x", argNames)
  
  ## extract arguments and argument names from given function
  args <- formals(FUN)
  fg <- argNames[1]
  bg <- argNames[2]
  extraArgs <- args[-seq_len(2L)]
  
  ## Define functions to apply
  a <- parse(text=deparse1(FUN))
  subFg <- parse_expr(paste0(x,"[",fg,"]"))
  subBg <- parse_expr(paste0(x,"[",bg,"]"))
  subList <- list(subFg,subBg)
  names(subList) <- c(fg,bg)
  newFUN <- do.call("substitute", list(a[[1]], subList))
  
  ## Get missing args
  margs <- names(extraArgs) %in% ls(envir = eframe)
  extraArgs[margs] <- 
    lapply(names(extraArgs[margs]), get, envir=eframe)
  
  missingArgs <- vapply(extraArgs, is, class="name", logical(1L))
  if(sum(missingArgs) > 0){
    missingArgNames <- glue_collapse(glue("`{names(missingArgs)}`"), sep = ", ")
    abort(glue("object(s) {missingArgNames} not found"))
  }
  
  ## substitute user-supplied variables
  b <- parse(text=deparse1(newFUN))
  newFUN <- do.call("substitute", list(b[[1]], extraArgs))
  newFUN <- as.function(eval(newFUN))
  extraArgNames <- paste(names(extraArgs), collapse=", ")
  header <- ifelse(length(extraArgs) > 0,
    sprintf("function(%s, %s, %s, %s)", x, fg, bg, extraArgNames),
    sprintf("function(%s, %s, %s)", x, fg, bg))
  body <- deparse1(body(newFUN))
  
  newFUN <- eval(parse(text=paste(header,body)))
  return(newFUN)
}

#' Internal for applyEnrichment
#' @inheritParams calcLoopEnrichment
#' @importFrom DelayedArray RegularArrayGrid DelayedArray blockApply
#' @importFrom utils capture.output
#' @importFrom stats median
#' @importFrom rlang inform
#' @noRd
.applyEnrichment <- function(x, fg, bg, FUN, nBlocks, verbose,
                                BPPARAM, ...) {
  cnts <- counts(x)
  
  ## Build array grid
  spacings <- dims <- dim(cnts)
  spacings[3] <- ceiling(spacings[3]/nBlocks)
  grid <- RegularArrayGrid(dims, spacings)
  
  ## Define functions to apply
  blockApplyFUN <- \(x) apply(x, c(3,4), \(x) FUN(x, fg, bg))
  combineFUN <- \(x) do.call("rbind", args=x)
  
  ## Apply in blocks
  blocks <- blockApply(x=cnts,
                       FUN=blockApplyFUN,
                       grid=grid,
                       verbose=verbose,
                       BPPARAM=BPPARAM)
  
  ans <- DelayedArray(combineFUN(blocks))
  
  if(!identical(dim(x),dim(ans))){
    dimExp <- glue_collapse(glue("{dim(x)}"), sep = " x ")
    dimAns <- glue_collapse(glue("{dim(ans)}"), sep = " x ")
    warn(c("Dimensions of output are not as expected.",
           "i"=glue("Expected dim: {dimExp}. Actual dim: {dimAns}."),
           "*"="Does `FUN` return more than one value?"))
  }
  return(ans)
}

#' Internal for calcLoopEnrichmentFromFiles
#' @inheritParams calcLoopEnrichment
#' @importFrom rlang abort inform
#' @importFrom glue glue
#' @importFrom DelayedArray RegularArrayGrid DelayedArray blockApply
#' @importFrom utils capture.output
#' @importFrom stats median
#' @noRd
.calcLoopEnrichmentFromFiles <- function(x, files, fg, bg, FUN, nBlocks, verbose,
                                         BPPARAM, ...) {
  
  ## Parameter parsing
  if (nBlocks <= 0) abort("`nBlocks` must be > 0.")
  
  if(class(fg)[1] != "MatrixSelection") {
    abort(c("`fg` must be a MatrixSelection object.",
            "x"=glue("`class(fg)` is {class(fg)}.")))
  }
  if(class(bg)[1] != "MatrixSelection") {
    abort(c("`bg` must be a MatrixSelection object.",
            "x"=glue("`class(bg)` is {class(bg)}")))
  }
  
  ## Determine resolution from x
  ## and ensure all pixels are same res.
  binSize <- .getBinSize(x)
  if (length(binSize) != 1L) {
    abort(c(glue("All interactions in `x` must be \\
                     the same width."),
            "i"="Check this with `width(x)`.",
            "i"="Set binSize with `binPairs(x, binSize)`."))
  }
  
  ## Check & show selection
  .checkBuffer(fg$buffer, bg$buffer)
  buffer <- fg$buffer
  fg <- fg$x
  bg <- bg$x
  if(verbose){ 
    .showMultiSelection(fg=fg, bg=bg, buffer=buffer)
  }
  
  ## Pull Hi-C matrices
  mr <- pixelsToMatrices(x, buffer)
  
  ## Pull count matrices
  ## Use nBlocks to set blockSize if not provided?
  iarr <- pullHicMatrices(mr, files, binSize, ...)
  
  ## Call apply enrichment to calculate scores on iarr
  newFun <- .modifyEnrichFun(FUN, eframe = parent.frame())
  ans <- .applyEnrichment(iarr, fg, bg, newFUN, nBlocks, verbose, BPPARAM)
  return(ans)
}

#' Internal for calcLoopEnrichmentFromIA
#' @inheritParams calcLoopEnrichment
#' @importFrom rlang abort inform
#' @importFrom DelayedArray RegularArrayGrid DelayedArray blockApply
#' @importFrom utils capture.output
#' @importFrom stats median
#' @noRd
.calcLoopEnrichmentFromIA <- function(x, fg, bg, FUN, nBlocks, verbose,
                                      BPPARAM, ...) {
  ## Parameter parsing
  if (nBlocks <= 0) abort("`nBlocks` must be > 0.")
  
  if(class(fg)[1] != "MatrixSelection") {
    abort(c("`fg` must be a MatrixSelection object.",
            "x"=glue("`class(fg)` is {class(fg)}.")))
  }
  if(class(bg)[1] != "MatrixSelection") {
    abort(c("`bg` must be a MatrixSelection object.",
          "x"=glue("`class(bg)` is {class(bg)}")))
  }
  
  ## If no foreground or background supplied,
  ## set to match buffer of count matrices
  buffer <- defaultBuffer(x)
  
  if(missing(fg)){
    fg <- selectCenterPixel(mhDist=1, buffer=defaultBuffer(x))  
  } 
  if(missing(bg)){
    bg <- selectTopLeft(n=4, buffer=defaultBuffer(x)) +
      selectBottomRight(n=4, buffer=defaultBuffer(x))
  }
  
  ## Check buffer & show selection
  .checkBuffer(buffer, fg$buffer)
  .checkBuffer(fg$buffer, bg$buffer)
  fg <- fg$x
  bg <- bg$x
  if(verbose){
    .showMultiSelection(fg=fg, bg=bg, buffer=buffer)
  }
  
  newFun <- .modifyEnrichFun(FUN)
  ans <- .applyEnrichment(x, fg, bg, FUN=newFun, nBlocks, verbose, BPPARAM)
  return(ans)
}

#' Calculate loop enrichment over background.
#'
#' Pulls Hi-C pixels and calculates the enrichment of
#' the selected foreground (`fg`) over the selected
#' background (`bg`).
#'
#' @param x GInteractions object or an InteractionArray object.
#' @param files Character file paths to `.hic` files. Required only if
#'  GInteractions object is supplied for x.
#' @param fg MatrixSelection object of matrix indices for the foreground.
#' @param bg MatrixSelection object of matrix indices for the background.
#' @param FUN Function with at least two parameters (i.e., `fg`, `bg`)
#'  defining how enrichment should be calculated. 
#'  Must produce a single value (numeric of length one).
#'  The first and second parameters must represent 
#'  fg and bg, respectively.
#' @param nBlocks Number of blocks for block-processing
#'  arrays. Default is 5. Increase this for large
#'  datasets. To read and process all data at once, set
#'  this value to 1.
#' @param verbose Boolean (TRUE or FALSE) describing
#'  whether to report block-processing progress.
#' @param BPPARAM Parallelization params (passed to
#'  `BiocParallel::bplapply()`). Default is the result
#'  of `BiocParallel::bpparams()`. Parallel processing
#'  is not available when `by=interactions`.
#' @param ... Additional arguments passed to
#'  `pullHicMatrices`. See ?[`pullHicMatrices`].
#' @returns A DelayedMatrix of enrichment scores
#'  where rows are interactions (i.e. loops) and
#'  columns are Hi-C files.
#'
#' @examples
#' ## Load marinerData
#' if (!require("marinerData", quietly = TRUE))
#'     BiocManager::install("marinerData")
#'
#' ## Read .hic file paths
#' hicFiles <- c(
#'     marinerData::LEUK_HEK_PJA27_inter_30.hic(),
#'     marinerData::LEUK_HEK_PJA30_inter_30.hic()
#' )
#' names(hicFiles) <- c("FS", "WT")
#'
#' ## Read in loops as GInteractions object
#' loops <-
#'     WT_5kbLoops.txt() |>
#'     setNames("WT") |>
#'     read.table(header=TRUE) |>
#'     as_ginteractions(keep.extra.columns=FALSE)
#'
#' ## Removes the "chr" prefix for compatibility
#' ## with the preprocessed hic files
#' GenomeInfoDb::seqlevelsStyle(loops) <- 'ENSEMBL'
#'
#' ## Expand binSize of loops
#' loops <- binPairs(x=loops, binSize=100e3)
#'
#' ## Calculate loop enrichment
#' calcLoopEnrichment(x=loops[1:10],
#'                    files=hicFiles)
#'
#' ## Customize different foreground/background
#' ## with selection functions
#' buffer <- 10 # choose pixel radius around center
#' fg <- selectCenterPixel(mhDist=seq(0,4), buffer=buffer)
#' bg <- selectCorners(n=6, buffer=buffer) +
#'     selectOuter(n=2, buffer=buffer)
#' 
#' ## Calculate loop enrichment
#' calcLoopEnrichment(x=loops[1:10],
#'                    files=hicFiles,
#'                    fg=fg,
#'                    bg=bg)
#'
#' @rdname calcLoopEnrichment
#' @export
setMethod("calcLoopEnrichment",
          signature(x="GInteractions",
                    files="character"),
          definition=.calcLoopEnrichmentFromFiles)

#' Calculate loop enrichment over background.
#'
#' @examples
#' ## Extract count matrices first
#' mats <- binPairs(loops[1:10],100e3) |>
#'   pixelsToMatrices(buffer=10) |>
#'     pullHicMatrices(
#'     files=hicFiles,
#'     binSize=100e3)
#' 
#' ## Calculate loop enrichment from count matrices 
#' calcLoopEnrichment(x = mats)
#'
#' @rdname calcLoopEnrichment
#' @export
setMethod("calcLoopEnrichment",
          signature(x="InteractionArray",
                    files="missing"),
          definition=.calcLoopEnrichmentFromIA)



