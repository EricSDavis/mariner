#' Parameter parsing function
#'
#' Place this function inside the parent function
#'
#' @param params `pgParams` object to override
#'  default arguments of parent function
#' @param defaultArgs List of defaults for each
#'  argument of parent function.
#' @param declaredArgs List of arguments to
#'  override all others
#'
#' @importFrom rlang abort
#' @importFrom stats na.omit
#'
#' @returns Final set of overridden arguments
#' @noRd
.parseParams <- function(params = params,
                         defaultArgs = formals(eval(match.call()[[1]])),
                         declaredArgs = lapply(match.call()[-1],
                                               eval.parent,
                                               n=2)) {

    ## Remove 'params' and '...' from defaultArgs and declaredArgs
    defaultArgs[["params"]] <- NULL
    declaredArgs[["params"]] <- NULL
    if ("..." %in% names(defaultArgs)) defaultArgs[["..."]] <- NULL
    if ("..." %in% names(declaredArgs)) declaredArgs[["..."]] <- NULL

    ## If pgParams are supplied override matching defaultArguments
    if (!is.null(params)) {

        if (!is(params, "pgParams")) {
            abort("`params` must be a `pgParams` object.")
        }

        ## Replace matching defaultArgs with params
        idx <- na.omit(sort(match(names(defaultArgs), names(params))))
        matchedParams <- params[idx]

        idx <- na.omit(match(names(params), names(defaultArgs)))
        defaultArgs[idx] <- matchedParams
    }

    ## Replace default args with declared args
    if (length(declaredArgs) != 0) {
        defaultArgs[names(defaultArgs) %in%
                        names(declaredArgs)] <-  declaredArgs
    }

    ## Set arguments without default to NULL
    unset <- unlist(lapply(defaultArgs, is.name))
    defaultArgs[unset] <-
        lapply(defaultArgs[unset], deparse) |>
        lapply(as.null)

    ## Return final set of overriden arguments
    return(defaultArgs)

}

#' Function to check zrange
#' @inheritParams plotMatrix
#' @importFrom rlang abort
#' @noRd
.checkZrange <- function(zrange) {

    if (is.null(zrange)) return(invisible(NULL))

    if (!is.vector(zrange)){
        abort("`zrange` must be a length 2 vector.")
    }

    if (length(zrange) != 2){
        abort("`zrange` must be a length 2 vector.")
    }

    if (!is.numeric(zrange)){
        abort("`zrange` must be a vector of two numbers.")
    }

    if (zrange[1] >= zrange[2]){
        abort("`zrange[2]` must be > `zrange[1]`.")

    }
}

#' Function that sets the zrange
#'
#' Sets the zrange based on the
#' matrix data if it is NULL.
#'
#' @param x `MatrixPlot` object with a slot for `zrange`
#'  and `data`.
#' @returns A `MatrixPlot` object with a set
#'  `zrange`.
#' @noRd
.setZrange <- function(x){

    ## Set zrange if it is null
    if (is.null(x$zrange)) {

        ## Set matrix counts variable for convenience
        cnts <- na.omit(as.vector(x$data))

        ## Only one value
        if (length(unique(cnts)) == 1) {
            x$zrange <- c(unique(cnts), unique(cnts))
        }

        ## Multiple values
        if (length(unique(cnts)) > 1) {

            ## Divergent data (center on 0 by default)
            if (all(c(-1, 1) %in% unique(sign(cnts))))
            {
                x$zrange <- c(floor(-max(abs(range(cnts)))),
                                     ceiling(max(abs(range(cnts)))))
            }

            ## Positive sequential data
            else if (!-1 %in% unique(sign(cnts)))
            {
                x$zrange <- c(0, ceiling(max(cnts)))
            }

            ## Negative sequential data
            else
            {
                x$zrange <- c(floor(min(cnts)), 0)
            }
        }
    }

    return(x)

}

#' Internal plotMatrix function
#' @inheritParams plotMatrix
#' @importFrom RColorBrewer brewer.pal
#' @importFrom plotgardener mapColors
#' @importFrom grid viewport unit grid.newpage gTree grid.ls rasterGrob
#'  addGrob grid.draw
#' @importFrom rlang inform
#' @importFrom grDevices colorRampPalette
#' @noRd
.plotMatrix <- function(data, params, x, y, width, height,
                        just, default.units, draw, palette, zrange,
                        na.color) {

    ## Parse parameters & Create Object ------

    ## Get parsed arguments
    parsedArgs <- .parseParams(
        params = params,
        defaultArgs = formals(eval(match.call()[[1]])),
        declaredArgs = lapply(match.call()[-1], eval.parent, n=2)
    )

    ## Evaluate parsed arguments
    parsedArgs <- lapply(parsedArgs, eval)

    ## Initialize object
    matrixPlot <- structure(
        .Data = c(parsedArgs,
                  list(
                      color_palette=parsedArgs$palette,
                      grobs=NULL
                  )
        ),
        class="MatrixPlot"
    )

    ## Set attributes for object
    attr(x=matrixPlot, which="plotted") <- parsedArgs$draw

    ## Parse units
    matrixPlot <-
        plotgardener:::defaultUnits(object=matrixPlot,
                                    default.units=parsedArgs$default.units)


    ## Parse matrix format ----------

    ## Realize DelayedMatrix as matrix
    if (is(matrixPlot$data, "DelayedMatrix")) {
        matrixPlot$data <- as.matrix(matrixPlot$data)
    }

    ## Map data ranges & colors ------

    ## Check for zrange errors
    .checkZrange(matrixPlot$zrange)

    ## Set zrange using matrix data if its null
    matrixPlot <- .setZrange(matrixPlot)

    ## Map values to colors
    colv <- mapColors(vector=as.vector(matrixPlot$data),
                      palette=matrixPlot$color_palette,
                      range=matrixPlot$zrange)

    ## Replace NA with na.color
    colv[is.na(colv)] <- na.color

    ## Format color vector back to matrix
    m <- matrix(data=colv,
                nrow=nrow(matrixPlot$data),
                ncol=ncol(matrixPlot$data))

    ## Viewports ---------------------------

    ## If placing information is provided but plot == TRUE,
    ## set up it's own viewport separate from pageCreate

    ## Not translating into page_coordinates
    if (is.null(matrixPlot$x) & is.null(matrixPlot$y)){

        vp <- viewport(x=unit(0.5, "npc"),
                       y=unit(0.5, "npc"),
                       height=unit(1, "snpc"),
                       width=unit(1, "snpc"),
                       clip="on",
                       just="center",
                       name="MatrixPlot1")

        if (matrixPlot$draw){
            grid.newpage()
        }

    } else {
        
        ## Get viewport name
        currentViewports <- plotgardener:::current_viewports()
        nVp <- length(grep("MatrixPlot", currentViewports))
        vp_name <- paste0("MatrixPlot", nVp + 1)

        ## Check that plotgardener page exists
        plotgardener:::check_page(
            "Use pageCreate() to make a plotgardener page to place a plot."
            )

        ## Convert coordinates into same units as page
        page_coords <- plotgardener:::convert_page(matrixPlot)


        ## Make viewport
        vp <- viewport(x=page_coords$x,
                       y=page_coords$y,
                       height=page_coords$height,
                       width=page_coords$width,
                       clip="on",
                       just=matrixPlot$just,
                       name=vp_name)
    }

    ## Handle graphical objects ----------------

    ## Initialize gTree for grobs
    assign("MatrixPlotGrobs", gTree(vp=vp), envir=plotgardener:::pgEnv)

    ## Get the number of MatrixPlot grobs
    nObjs <-
        grep(pattern="MatrixPlotGrobs",
             x=grid.ls(print=FALSE, recursive=FALSE)) |>
        length()

    ## Assign name to grob
    name <- paste0("MatrixPlot", nObjs + 1)

    ## Make grobs
    mpRaster <- rasterGrob(image=m, interpolate=FALSE, name=name)

    ## Assign grobs to gTree
    assign(x="MatrixPlotGrobs",
           value=addGrob(gTree=get(x="MatrixPlotGrobs",
                                   envir=plotgardener:::pgEnv),
                         child=mpRaster),
           envir=plotgardener:::pgEnv)

    ## Add grobs to object
    matrixPlot$grobs <- get("MatrixPlotGrobs", envir=plotgardener:::pgEnv)

    ## Plot grobs
    if (matrixPlot$draw) {
        grid.draw(matrixPlot$grobs)
    }

    ## Return object --------------------------
    inform(paste0("MatrixPlot[", mpRaster$name, "]"))
    invisible(matrixPlot)
}

#' Plot matrix
#'
#' Used to plot single or aggregate matrix
#' such as aggregate peak analysis.
#'
#' @param params Optional `pgParams` object
#'  containing relevant function parameters.
#' @param data `DelayedMatrix`, `matrix`, list of matrices,
#'  or 3 column `data.frame` of APA results.
#' @param x Numeric or unit object specifying
#'  the x-location of plot.
#' @param y Numeric or unit object specifying
#'  the y-location of plot.
#' @param width Numeric or unit object specifying
#'  the width of plot.
#' @param height Numeric or unit object specifying
#'  the height of plot.
#' @param just String or numeric vector specifying
#'  the justification of the viewport relative to
#'  its (x, y) location.
#' @param default.units String indicating the
#'  default units to use if `x`, `y`, `width`, or
#'  `height` are only given as numeric vectors.
#' @param draw Logical value indicating whether
#'  graphics output should be produced.
#' @param palette `colorRampPalette` function
#'  to use for mapping values to colors.
#' @param zrange Vector of length 2;
#'  max and min values to set color scale
#' @param na.color String indicating the color
#'  to use for mapping NA values.
#'
#' @return Function will draw a color-mapped matrix
#'  and return an S3 object of class `MatrixPlot`.
#'
#' @examples
#'
#' library(plotgardener)
#' library(RColorBrewer)
#'
#' ## Create divergent matrix ####
#' m <- matrix(data=rnorm(n=21*21, mean=0, sd=2), nrow=21, ncol=21)
#'
#' ## Define parameters
#' p <- pgParams(width=3, height=3, default.units="inches")
#'
#' ## Create page
#' pageCreate(params=p)
#'
#' ## Plot apa
#' plot <- plotMatrix(data=m,
#'                    x=p$width/2,
#'                    y=p$height/2,
#'                    width=p$width*0.5, height = p$width*0.5,
#'                    just=c("center", "center"),
#'                    palette=colorRampPalette(c("blue", "white", "red")),
#'                    zrange=NULL)
#'
#' ## Annotate legend
#' annoHeatmapLegend(plot=plot,
#'                   x=2.3,
#'                   y=0.75,
#'                   width=0.1,
#'                   height=0.75)
#'
#'
#' ## Create sequential matrix
#' m <- matrix(data=sample(0:100, 21*21, replace=TRUE), nrow=21, ncol=21)
#'
#' ## Define parameters
#' p <- pgParams(width=3, height=3, default.units="inches")
#'
#' ## Create page
#' pageCreate(params=p)
#'
#' ## Plot apa
#' plot <- plotMatrix(data=m,
#'                    x=p$width/2,
#'                    y=p$height/2,
#'                    width=p$width*0.5,
#'                    height=p$width*0.5,
#'                    just=c("center", "center"),
#'                    palette=colorRampPalette(c("white", "dark red")),
#'                    zrange = NULL)
#'
#' ## Annotate legend
#' annoHeatmapLegend(plot=plot,
#'                   x=2.3,
#'                   y=0.75,
#'                   width=0.1,
#'                   height=0.75)
#'
#' @rdname plotMatrix
#' @export
setMethod("plotMatrix",
          signature(data="DelayedMatrix_OR_matrix"),
          definition=.plotMatrix)
