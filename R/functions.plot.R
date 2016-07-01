#!/usr/bin/Rscript
#=======================================================================
#
#   Function for related to plotting
#
#=======================================================================

# load some useful libraries
require(stringr)
require(RColorBrewer)
require(gridExtra)      # for dotplot with denisty at axis
require(gplots)         # for heatmap.2
require(ggplot2)        # for nice plots



#-----------------------------------------------------------------------
# get group labels from counts by grouping all elements >= maxCounts together.
#-----------------------------------------------------------------------
getGroupsByCounts <- function(counts, maxCount=10){
    groups = as.character(counts)
    groups[counts >= maxCount] = paste0(maxCount, '+')
    return(groups)
}

#-----------------------------------------------------------------------
# function to plot error bars (from: http://monkeysuncle.stanford.edu/?p=485 )
#-----------------------------------------------------------------------
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
    arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

#-----------------------------------------------------------------------
# Modified boxplot function with rotated x-labels 
#-----------------------------------------------------------------------
my.boxplot <- function(valueList, offset=0.075, names=names(valueList), log = "", returnBp=FALSE, ...){

    op = par(mar = c(7, 4, 4, 2) + 0.1)
    boxP = boxplot(valueList, names=NA, log=log, ...)
    
    if (log==""){
        yValueRange = max(unlist(valueList), na.rm=TRUE) - min(unlist(valueList), na.rm=TRUE)
        yTextPos = par("usr")[3] - offset*yValueRange

    }else{
        yValueRange = log(max(unlist(valueList))) - log(min(unlist(valueList)))
        yTextPos = par("usr")[3] + offset * yValueRange
    }
    
    
    text(seq(length(valueList)), yTextPos, srt = 45, adj = 1, labels = names, xpd = TRUE)

    par(op)
    if(returnBp){
        return(boxP)
    }
}

#-----------------------------------------------------------------------
# Modified barplot function with rotated x-labels 
#-----------------------------------------------------------------------
my.barplot <- function(values, offset=0.1, beside=TRUE, addValues=FALSE, customValues=values, digits=3, yMax = ifelse(addValues, 1.2*max(values), max(values)), names=NULL, axis.lwd=0, mar=c(7, 4, 4, 2) + 0.1, srt=45, adj=1, ...){
    
    # if input is single vector, convert it to matrix:
    if( class(values) == "numeric"){
        valuesMat=t(matrix(values))
    }else{
        valuesMat = values
    }

    op = par(mar = mar)

    bp = barplot(values, beside=beside, names.arg=rep(NA, ncol(valuesMat)), ylim=c(0,yMax), ...)

    # add x-axis labels
    yValueRange = max(values)
    yTextPos = par("usr")[3] - offset*yValueRange

    if( class(values) == "numeric"){
        xpos=bp
    }else{
        xpos = colMeans(bp)
    }


    axis(1, at=xpos , labels = FALSE, lwd=axis.lwd, lwd.ticks=axis.lwd)
    text(xpos , yTextPos, srt = srt, adj=adj, labels = names, xpd = TRUE)
    
    if (is.numeric(customValues)){
        if (digits){
            valueLabels = signif(customValues, digits)
        }else{
            valueLabels = customValues
        }
    }else{
        valueLabels = customValues
    }

    if(addValues){
        text(bp, values, valueLabels , pos=3)
    }

    par(op)
    return(bp)
}

#-----------------------------------------------------------------------
# function to add p-values to a barplot 
#-----------------------------------------------------------------------
add_pval_two_pairs <- function(bp, values, pvalues, useStars=FALSE, offset=.1*max(values), min_pval=10^-100, digits=2, ...){


    stopifnot(nrow(values) == 2)

    # if option given replace actual p-values by star symbols
    if (useStars){
        stars = rep(NA, length(pvalues))
        stars[pvalues >=0.05] = "n.s."
        stars[pvalues < 0.05] = "*"
        stars[pvalues < 0.01] = "**"
        stars[pvalues < 0.001] = "***"
        stars[pvalues < 2.2*10^-16] = "****"
    }
    
    vertLen = .2 * offset
    
    ypos = apply(values, 2, max) + offset
        
    for (j in 1:ncol(values)){
        
        # add horizontal line:
        lines(c(bp[1,j], bp[2,j]), c(ypos[j], ypos[j]), ...)

        # add vertical lines
        lines(c(bp[1,j],bp[1,j]), c(ypos[j], ypos[j]-vertLen), ...)
        lines(c(bp[2,j],bp[2,j]), c(ypos[j], ypos[j]-vertLen), ...)
        
        # add p-value:
        if (useStars){
            p_str=stars[j]
        }else{
            if (pvalues[j] >= min_pval | is.na(pvalues[j]) ) {
                p_str = paste("p=", signif(pvalues[j], digits), sep="")
            }else{
                p_str = paste("p<", signif(min_pval,2), sep="" )
            }
        }
        text(mean(bp[,j]), ypos[j], p_str, pos=3, ...)
     
    }
}

#-----------------------------------------------------------------------
# function to add p-values to a barplot 
#-----------------------------------------------------------------------
add_pval_all_allPairs <- function(xpos, ymax, pval_mat, offset_fac=.125, min_pval=10^-100, digits=3, ...){
    ####################################################################
    # Annotates a boxplot or barplot with pvalues:
    # xpos := the barplot object, or a list of x-axis centers
    # ymax := the maximal y-value in the plot
    # pval_mat := a (upper triangular) matrix with pairwise p-values
    # offset_frac := fraction of maximal height value used space to the lines
    # min_pval := the minimal p-value to be shown 
    #... := other graphical parameter like cex=1.3
    ####################################################################

    # round p-values to 'digits' significant digits
    pval_mat = signif(pval_mat, digits)

    # take fraction of cex as vertical offset for annotation line
    vert_offset = offset_fac * ymax

    # the vertical height of 'spanner' are 25% of the offset
    vert = .25 * vert_offset

    n = length(xpos)

    lineCnt = 0
    

    for (i in seq(1, n-1)){
                
        # iterate over all other bars:
        for (j in seq(i+1, n)){
        
            # get ploting height
            y = ymax + lineCnt * vert_offset + vert
            
            # add horizontal line:
            lines(c(xpos[i], xpos[j]), c(y, y) + vert, ...)
            
            # add vertical lines
            lines(c(xpos[i],xpos[i]), c(y, y+vert), ...)
            lines(c(xpos[j],xpos[j]), c(y, y+vert), ...)
            
            # add p-value:
            if (pval_mat[i,j] >= min_pval | is.na(pval_mat[i,j]) ) {
                p_str = paste("p=", pval_mat[i,j], sep="")
            }else{
                p_str = paste("p<", signif(min_pval,2), sep="" )
            }
            text(mean(c(xpos[i], xpos[j])), y+vert, p_str, pos=3, ...)
            
            lineCnt = lineCnt+1
        }
    }
}
# add_pval_all(seq(length(bp))[5:8], sapply(logFC_list, max)[[5:8]], f_mat)

#-----------------------------------------------------------------------
# plot strand-specific counts around boundary
#-----------------------------------------------------------------------
plotStrandedCountsAroundBoundaries <- function(strandList, outFile, windowSize=2*10^5, nbins=20, strandCols="gray", main="", ylab="Average counts", cex=1.5, lwd=3, ...){
    
    if (length(strandCols) == 1){
        strandCols = rep(strandCols, 2)
    }

    ymax = max(unlist(strandList))
    ymin = min(unlist(strandList))

    cairo_pdf(outFile)

        par(cex=cex, lwd=lwd)
        plot(strandList[["+"]], type="b", pch=-8853, ylim=c(ymin, ymax), col=strandCols[1], xaxt = "n", ylab=ylab, xlab="", main=main, ...)
        points(strandList[["-"]], type="b", pch=-8854, col=strandCols[2], ...)
        
        xlab = c(-.5, -.25, 0, .25, .5) * windowSize / 10^3
        xlabAT = seq(0,nbins, nbins/4)+.5
        axis(1, at=xlabAT, labels=xlab, line=1, cex=cex, lwd=lwd)
        mtext("Distance from TAD boundary (kb)", side=1, line=4, cex=cex)
        
        par(xpd=TRUE)
        valueRange = ymax - ymin
        yMin = ymin - .06 * valueRange
        yLow = ymin - .1 * valueRange
        polygon(xlabAT[c(1, 1,3, 3)], c(yLow, yMin, yMin, yLow), col="black")
        par(xpd=FALSE)
        
        abline(v=nbins/2+.5, lwd=lwd, col="gray", lty=3)
    
    dev.off() 
}

#-----------------------------------------------------------------------
# apply function to subset of data frame according to a specific group column
# Usefull for ggplot2 style of data frames to generate some annotations
#-----------------------------------------------------------------------
applyToSubset <- function(df, fun, valueCol, groupCol, ...){

    sapply(unique(df[,groupCol]), function(grp){
        fun(df[df[,groupCol] == grp, valueCol], ...)
    })
}


#-----------------------------------------------------------------------
# Dotplot with density at axis
# according to: http://www.r-bloggers.com/ggplot2-cheatsheet-for-visualizing-distributions/
# use grid::grid.draw() to plot the return value
#-----------------------------------------------------------------------
dotplotWithDensityLogXY <- function(plotDF, xvar, yvar, zvar, COL=c("orange", "purple"), ALPHA=.5, fit=FALSE, xlog=TRUE, ylog=TRUE, ylab=yvar, xlab=xvar){

    #placeholder plot - prints nothing at all
    empty <- ggplot()+geom_point(aes(1,1), colour="white") +
         theme(                              
           plot.background = element_blank(), 
           panel.grid.major = element_blank(), 
           panel.grid.minor = element_blank(), 
           panel.border = element_blank(), 
           panel.background = element_blank(),
           axis.title.x = element_blank(),
           axis.title.y = element_blank(),
           axis.text.x = element_blank(),
           axis.text.y = element_blank(),
           axis.ticks = element_blank()
         )
    
    #scatterplot of x and y variables
    if (fit){
        scatter <- ggplot(plotDF, aes_string(x=xvar, y=yvar, color=zvar, fill=zvar)) + 
          geom_point(alpha=ALPHA, size=.5, shape=20)  +
          geom_smooth(alpha=ALPHA) + 
          scale_color_manual(values = COL) + scale_fill_manual(values = COL) + theme_bw() + 
          theme(legend.position=c(1,1),legend.justification=c(-.01,-.01), plot.margin =rep(grid::unit(c(0,0,1,1), "cm"),4)) + labs(x=xlab, y=ylab)
    }else{
        scatter <- ggplot(plotDF, aes_string(x=xvar, y=yvar)) + 
          geom_point(aes_string(color=zvar), alpha=ALPHA) + 
          scale_color_manual(values = COL) + scale_fill_manual(values = COL) + theme_bw() + 
          theme(legend.position=c(1,1),legend.justification=c(-.01,-.01), plot.margin =rep(grid::unit(c(0,0,1,1), "cm"),4)) + labs(x=xlab, y=ylab) 
    }
    #marginal density of x - plot on top
    plot_top <- ggplot(plotDF, aes_string(xvar, fill=zvar)) + 
      geom_density(alpha=.5) + theme_bw() + #+ scale_x_log10() 
      scale_fill_manual(values = COL) + 
      theme(legend.position = "none", axis.title.x = element_blank(), plot.margin = rep(grid::unit(c(1,0,0,1), "cm"),4))
      
    #marginal density of y - plot on the right
    plot_right <- ggplot(plotDF, aes_string(yvar, fill=zvar)) + 
      geom_density(alpha=.5) + coord_flip() + theme_bw() +  # + scale_x_log10()
      scale_fill_manual(values = COL) + 
      theme(legend.position = "none", axis.title.y = element_blank(), plot.margin = rep(grid::unit(c(0,1,1,0), "cm"),4))
    
    if(xlog){
        scatter <- scatter + scale_x_log10()
        plot_top <- plot_top + scale_x_log10()
    }
    if(ylog){
        scatter <- scatter + scale_y_log10()
        plot_right <- plot_right + scale_x_log10()
    }
    
    #arrange the plots together, with appropriate height and width for each row and column
    arrangeGrob(plot_top, empty, scatter, plot_right, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
}
#grid::grid.draw(dotplotWithDensityLogXY(plotDFcloseHiC, "dist", "HiCraw", "group", COL))

#-----------------------------------------------------------------------
# reverse a data.frame row-wise. 
#-----------------------------------------------------------------------
revDF <- function(df){
    return(df[nrow(df):1,])
}


#-----------------------------------------------------------------------
# Add pictures as faced grid panel labels (plot with  grid.draw(g))
#-----------------------------------------------------------------------
addPictureLabels <- function(p, figPath){
    
    # load required libraries
    require(png)
    require(grid)
    require(gtable)

    p <- p + theme(
        strip.background = element_rect(colour = "black", fill = "white", size=.5)
        )
    
    # parse PNG figures
    figs <- sapply(figPath, readPNG)

    g <- ggplot_gtable(ggplot_build(p))

    strips <- grep("strip", g$layout$name)
    
    new_grobs <- lapply(figs, rasterGrob, width=.95, height=.45)

    
    g <- with(g$layout[strips,],
              gtable_add_grob(g, new_grobs,
                              t=t, l=l, b=b, r=r, name="strip_figs") )
    return(g)
}

