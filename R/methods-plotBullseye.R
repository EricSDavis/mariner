#' Define a function to draw a bullseye plot
#' @param x matrix
#' @param title Title for plot
#' @param nl Number of loops to divide matrix by
#' @param cols passed to `colorRampPalette` function
#'  to use for mapping values to colors.
#'
#' @importFrom colourvalues color_values
#' @importFrom grDevices colorRamp
#' @importFrom graphics polygon
#' @importFrom graphics par
#' @importFrom graphics rect
#' @importFrom graphics segments
#' @importFrom graphics text
#'
#' @examples
#' plotBullseye(x = matrix(1:(21*21), 21, 21), title="test")
#' @noRd
plotBullseye <- function(x,
                         title = "My plot",
                         nl = 1,
                         zrange=NULL,
                         cols = c("white", "red")){

    ## Collect points by manhattan distance ####
    ## Read in and clean data
    dat <- as.matrix(x)

    ## Set scale to 0 and
    ## divide aggregate values by number of loops
    ## to get average loop strength
    dat <- dat - min(dat)
    dat <- dat/nl

    ## Define color function
    col_fun <- colorRampPalette(cols)

    ## Adjust scale to zrange ####
    if(is.null(zrange)){
        ## Convert to colorscale
        grps <- cut(dat, seq(min(dat), max(dat), length.out = 1000),
                    include.lowest = T)
        cdat <- matrix(col_fun(1000)[grps], nrow = 21, ncol = 21)

        ## Create color values for legend scale
        scale_cols <- col_fun(1000)

    } else {
        ## Rescale data to zrange
        dat[dat >= zrange[2]] <- zrange[2]
        dat[dat <= zrange[1]] <- zrange[1]

        ## Convert to colorscale
        grps <- cut(dat, seq(zrange[1], zrange[2],
                             length.out = 1000), include.lowest = T)
        cdat <- matrix(col_fun(1000)[grps], nrow = 21, ncol = 21)

        ## Create color values for legend scale
        scale_cols <- colourvalues::color_values(
            seq(zrange[1],zrange[2], zrange[2]*0.001),
            palette = colorRamp(cols)((seq(1,10))/10)
        )

        scale_cols <- col_fun(1000)

    }

    ## Define lists to collect information
    pos1 <- c(11)
    pos2 <- c(11)
    shell <- c(0)
    value <- c(dat[11,11])
    color <- c(cdat[11,11])

    ## Define center point and shell layers
    n1 = 11
    n2 = 11
    j = 10 # total number of shells

    ## Walk along each diagonal (counterclock-wise) to get points
    for (n in seq(1,j)){
        for(i in 0:(n-1)){
            ## NW diagonal
            pos1 <- c(pos1, n1-n+i)
            pos2 <- c(pos2, n2-i)
            shell <- c(shell, n)
            value <- c(value, dat[n1-n+i, n2-i])
            color <- c(color, cdat[n1-n+i, n2-i])
        }
        for(i in 0:(n-1)){
            ## SW diagonal
            pos1 <- c(pos1, n1+i)
            pos2 <- c(pos2, n2-n+i)
            shell <- c(shell, n)
            value <- c(value, dat[n1+i, n2-n+i])
            color <- c(color, cdat[n1+i, n2-n+i])
        }
        for(i in 0:(n-1)){
            ## SE diagonal
            pos1 <- c(pos1, n1+n-i)
            pos2 <- c(pos2, n2+i)
            shell <- c(shell, n)
            value <- c(value, dat[n1+n-i, n2+i])
            color <- c(color, cdat[n1+n-i, n2+i])
        }
        for(i in 0:(n-1)){
            ## NE diagonal
            pos1 <- c(pos1, n1-i)
            pos2 <- c(pos2, n2+n-i)
            shell <- c(shell, n)
            value <- c(value, dat[n1-i, n2+n-i])
            color <- c(color, cdat[n1-i, n2+n-i])
        }
    }

    ## Bind into data frame
    df <- data.frame(cbind(pos1, pos2, shell, value, color),
                     stringsAsFactors = F, row.names = NULL)

    ## Change columns to numeric
    df[,seq(1,3)] <- apply(df[,seq(1,3)], 2, as.numeric)

    ## Define function to plot shell layers ###
    plotShell <- function(cx=0.5,
                          cy=0.5,
                          r1=1*0.025,
                          r2=2*0.025,
                          slices=4,
                          cols=seq(1,length(breaks)),
                          shell=1){
        degrees = 360*slices
        res = 2*pi/slices

        xarc1 <- cos(seq((shell*2-1)*pi/slices,
                         (shell*2-1)*pi/slices+2*pi,
                         length.out = degrees))*r1+cx
        yarc1 <- sin(seq((shell*2-1)*pi/slices,
                         (shell*2-1)*pi/slices+2*pi,
                         length.out = degrees))*r1+cy

        xarc2 <- cos(seq((shell*2-1)*pi/slices,
                         (shell*2-1)*pi/slices+2*pi,
                         length.out = degrees))*r2+cx
        yarc2 <- sin(seq((shell*2-1)*pi/slices,
                         (shell*2-1)*pi/slices+2*pi,
                         length.out = degrees))*r2+cy

        breaks = seq(1, degrees, (degrees/slices))

        for(i in seq(1,(length(breaks)-1))){
            polygon(x = c(xarc1[seq(breaks[i],breaks[i+1])],
                          rev(xarc2[seq(breaks[i],breaks[i+1])])),
                    y = c(yarc1[seq(breaks[i],breaks[i+1])],
                          rev(yarc2[seq(breaks[i],breaks[i+1])])),
                    col = cols[i], border = NA)
        }
        polygon(x = c(xarc1[seq(breaks[i+1],length(xarc1))],
                      rev(xarc2[seq(breaks[i+1],length(xarc2))])),
                y = c(yarc1[seq(breaks[i+1],length(yarc1))],
                      rev(yarc2[seq(breaks[i+1],length(yarc2))])),
                col = cols[i+1], border = NA)
    }

    ## Create bulls-eye plot ####
    ## Define variables for plotting
    len = 0.045 # radius length
    center = 0.5 # center of the plot (could define cx and cy separately if desired)

    par(pty="s") ## keep square asp ratio when not plotting axes
    plot(center, center, 'n',
         xlim=c(0,1.2), ylim=c(0,1.2),
         axes=F, main = title,
         asp=1, xlab="", ylab="")

    ## Plot center pixel
    polygon(x = cos(seq(0,2*pi,length.out = 360))*len+center,
            y = sin(seq(0,2*pi,length.out = 360))*len+center,
            col = df$color[1], border = NA)

    ## Plot each shell (increasing by 4 each time)
    shells = 10
    seg = 4
    for(n in 0:shells){
        plotShell(cx = center, cy = center,
                  r1 = (n+1)*len, r2 = (n+2)*len,
                  slices = seg, cols = df$color[df$shell == n+1],
                  shell=n+1)
        seg <- seg + 4
    }

    ## Add outer circle border
    polygon(x = cos(seq(0,2*pi,length.out = 360))*(shells+1)*len+center,
            y = sin(seq(0,2*pi,length.out = 360))*(shells+1)*len+center,
            border = "black")

    ## Add scale bar ####
    ## Create breaks to plot the colors
    breaks <- seq(center-0.05, 0.95, length.out = length(scale_cols))

    ## Draw rectangles for each color in the scale
    for(i in seq(1,(length(scale_cols)-1))){
        rect(xleft = 1.025, xright = 1.1,
             ybottom = breaks[i], ytop = breaks[i+1],
             col=scale_cols[i], border = NA)
    }

    ## Add text for color scale
    if(is.null(zrange)){
        labels <- round(seq(min(dat), max(dat), length.out = 5))
    } else {
        labels <- round(seq(zrange[1], zrange[2], length.out = 5))
    }

    positions <- seq(center-0.05, 0.95, length.out = 5)

    for(i in seq(1,length(positions))){
        segments(x0=1.025, x1=1.045, y0=positions[i], y1=positions[i])
        segments(x0=1.08, x1=1.1, y0=positions[i], y1=positions[i])
        text(x = 1.11, y=positions[i], labels = labels[i],
             cex = 0.6, adj = c(0, 0.5))
    }

    ## Plot outer border
    rect(xleft = 1.025, xright = 1.1, ybottom = center-0.05, ytop = 0.95)

    ## Return dataframes
    return(
        list(
            mat = dat,
            df = df
        )
    )
}

