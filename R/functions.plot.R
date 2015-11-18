#!/usr/bin/Rscript
#=======================================================================
#
#   Function for related to plotting
#
#=======================================================================

# load some useful libraries
require(stringr)
require(RColorBrewer)
require(gplots) 
require(ggplot2)        # for nice plots
require(gridExtra)      # for dotplot with denisty at axis

#-----------------------------------------------------------------------
# Hi-C plot
# to rotate plot see: http://stackoverflow.com/questions/3792803/is-it-possible-to-rotate-a-plot-in-r-base-graphics
#-----------------------------------------------------------------------
myHiCplot <- function(intMat, outFile="test_plot.pdf"){
    
    require(grid)   # for customized Hi-C plot

    op <- par(mar = rep(0, 4))
    graphics::image(x=1:10, y=1:10, z=as.matrix(intMat), xlab="", ylab="")
    par(op)

    cap <- grid.cap()
    dev.off()
    
    pdf(outFile)
        
        grid.newpage()
        
        # add the heatmap rotated by 45 degrees
        grid.raster(cap, width=1/sqrt(2), vp=viewport(angle=45))
    
        # plot with rectangular to overwrite lower part of the heatmap
        grid.raster(1, x=c(0, 0.5), y=c(0,0), width=1, height=.5, just=c("left", "bottom") )    
    
    dev.off()
    
    #TODO: check #?viewport
    
}
#~ intMat = intdata(subMapNorm.binned)[1:10,1:10]

#-----------------------------------------------------------------------
# calculate the standard error of the mean for an input vector
#-----------------------------------------------------------------------
std.err <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))

#-----------------------------------------------------------------------
# calculate for each column in the matrix the standard error of the mean
#-----------------------------------------------------------------------
standardError <- function(mat, ...){
    apply( mat, 2, function(x){sd(x, ...) / sqrt(length(x))} )
}

#-----------------------------------------------------------------------
# get the density values from a histogram plot
#-----------------------------------------------------------------------
getDiscreteDensity <- function(counts, breaks=20){
    unlist(hist(x=counts, breaks=breaks, plot=FALSE)["density"])
}

#-----------------------------------------------------------------------
# get group labels from counts by grouping all elements >= maxCounts together.
#-----------------------------------------------------------------------
getGroupsByCounts <- function(counts, maxCount=10){
    groups = as.character(counts)
    groups[counts >= maxCount] = paste0(maxCount, '+')
    return(groups)
}

#-----------------------------------------------------------------------
# add wiskers to an barplot
#-----------------------------------------------------------------------
addWiskers <- function(xes, height, wiskers, length=.1, ...){
    arrows(xes, height, xes, height+wiskers, angle=90, length=length, ...)
    arrows(xes, height, xes, height-wiskers, angle=90, length=length, ...)
}

#-----------------------------------------------------------------------
# function to plot error bars (from: http://monkeysuncle.stanford.edu/?p=485 )
#-----------------------------------------------------------------------
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
    arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

#-----------------------------------------------------------------------
# add polygon to indicate deviation or error ranges
#-----------------------------------------------------------------------
addPolygon <- function(xes, height, wiskers, col, alpha=.5, ...){
    polygon(c(xes, rev(xes)), c(height+wiskers, rev(height-wiskers)), 
        col=add.alpha(col, alpha), border = NA, ...)
}


#-----------------------------------------------------------------------
# Add an alpha value to a colour 
# SOURCE: http://www.magesblog.com/2013/04/how-to-change-alpha-value-of-colours-in.html)
#-----------------------------------------------------------------------
add.alpha <- function(col, alpha=1){
    if(missing(col))
    stop("Please provide a vector of colours.")
    apply(sapply(col, col2rgb)/255, 2,
    function(x)
    rgb(x[1], x[2], x[3], alpha=alpha))
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
# function to convert P-values into star notation 
#-----------------------------------------------------------------------
pValToStars <- function(pvalues, nsLabel="n.s."){
    stars = rep(NA, length(pvalues))
    stars[pvalues >=0.05] = nsLabel
    stars[pvalues < 0.05] = "*"
    stars[pvalues < 0.01] = "**"
    stars[pvalues < 0.001] = "***"
    stars[pvalues < 2.2*10^-16] = "****"
    return(stars)
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
#~         for (j in seq(n, i+1)){
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
# plot data around boundary
#-----------------------------------------------------------------------
plotCountsAroundBoundaries <- function(meanCount, outFile, windowSize=2*10^5, nbins=20, col="gray", main="", ylab="Average counts", cex=1.5, lwd=3, type="o", pch=20, ...){
    cairo_pdf(outFile)

        par(cex=cex, lwd=lwd)
        plot(meanCount, type=type, pch=pch, col=col, xaxt = "n", 
           ylab=ylab, xlab="", main=main, ...)
        
        xlab = c(-.5, -.25, 0, .25, .5) * windowSize / 10^3
        xlabAT = seq(0,nbins, nbins/4)+.5
        axis(1, at=xlabAT, labels=xlab, line=1, cex=cex, lwd=lwd)
        mtext("Distance from TAD boundary (kb)", side=1, line=4, cex=cex)
        
        par(xpd=TRUE)
        valueRange = max(meanCount) - min(meanCount)
        yMin = min(meanCount) - .06 * valueRange
        yLow = min(meanCount) - .1 * valueRange
        polygon(xlabAT[c(1, 1,3, 3)], c(yLow, yMin, yMin, yLow), col="black")
        par(xpd=FALSE)
        
        abline(v=nbins/2+.5, lwd=lwd, col="gray", lty=3)
    
    dev.off() 
}

#-----------------------------------------------------------------------
# plot multiple count data around boundary
#-----------------------------------------------------------------------
plotMultipleCountsAroundBoundaries <- function(countList, outFile, windowSize=2*10^5, nbins=20, cols=gray.colors, main="", ylab="Average counts", cex=1.5, lwd=3, pch=20, type="o", ...){
    cairo_pdf(outFile)
        
        ymax = max(sapply(countList, max))
        ymin = min(sapply(countList, min))
        regCols = cols(length(countList))
        
        par(cex=cex, lwd=lwd)
        plot(countList[[1]], type="n", ylim=c(ymin, ymax), col=col, xaxt = "n", ylab=ylab, xlab="", main=main, ...)

        # add data to empty plot
        for (i in 1:length(countList)){
            
            lines(countList[[i]], type=type, col=regCols[i], ...)
        
        }

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
#plotMultipleCountsAroundBoundaries(countList, "test.out.pdf")

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
# plot counts in Bins of relative size along window around TAD
#-----------------------------------------------------------------------
plotCountsAroundTADwindow <- function(meanCount, outFile, errorCount=NULL, nbins=20, col="darkgray", main="", ylab="Average counts", cex=1.5, lwd=3, pch=20, type="o", ...){
    cairo_pdf(outFile)

        par(cex=cex, lwd=lwd)
        if (!is.null(errorCount)){
            valueRange = range(c(meanCount+errorCount, meanCount-errorCount))
        }else{
            valueRange = range(meanCount)
        }
        
        plot(meanCount, type="n", xaxt = "n", ylim=valueRange, pch=pch, xaxt = "n", ylab=ylab, xlab="", main=main, ...)
        
        xlab = c("-50%", "start", "TAD", "end", "+50%")
        xlabAT = seq(0,nbins, nbins/4)+.5
        axis(1, at=xlabAT, labels=xlab, line=1, cex=cex, lwd=lwd)
        
        if (!is.null(errorCount)){
            addWiskers(1:nbins, meanCount, errorCount, length=.05, col="gray")
            addPolygon(1:nbins, meanCount, errorCount, col="gray")
        }
        
        lines(meanCount, type=type, pch=pch, col=col)
        
        par(xpd=TRUE)
        yMin = valueRange[1] - .06 * (valueRange[2] - valueRange[1])
        yLow = valueRange[1] - .1 * (valueRange[2] - valueRange[1])
        polygon(xlabAT[c(2, 2,4, 4)], c(yLow, yMin, yMin, yLow), col="black")
        par(xpd=FALSE)
        
        abline(v=xlabAT[c(2,4)], lwd=lwd, col="gray", lty=3)
    
    dev.off() 
}

#-----------------------------------------------------------------------
# Plot stacked area plot around TAD window bins
# TODO: FIXME: arguments seems to be not passed, can not find myGroups
#-----------------------------------------------------------------------
plotStackedAreaAroundTAD <- function(bins, values, myGroups, colValues, outFile, nbins=20, ylab="Average counts"){

    # create proper x-axis tick labels
    xlabs = rep("", nbins+1)
    names(xlabs) = as.character(1:(nbins+1))
    xlabs[seq(0,nbins, nbins/4)+1] = c("-50%", "start", "TAD", "end", "+50%")
    
    # mark TAD body with rectengular
    TADrect <- data.frame(xmin=nbins/4+1, xmax=3*nbins/4+1, ymin=-Inf, ymax=Inf)
    
    plotTheme = theme_bw(base_size = 20) + theme(legend.position = "bottom", legend.title=element_blank())
    plotTADrect = geom_rect(data=TADrect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),color="grey20", alpha=0.25, inherit.aes = FALSE)
    plotColor = scale_color_manual(values=colValues) 
    plotFill = scale_fill_manual(values=colValues)
    
    p <- ggplot(NULL, aes( bins, values))  + 
    geom_area(aes(color = myGroups, fill= myGroups), alpha=.5, position = 'stack') + scale_x_discrete(labels=xlabs) + labs(x="", y=ylab, aesthetic="") + plotTheme + plotTADrect + plotColor + plotFill

    ggsave(p, file=outFile, width=7, height=7)

}


#-----------------------------------------------------------------------
# Plot ranges 
# (Source: http://master.bioconductor.org/help/course-materials/2009/GenentechNov2009/Module7/IRanges.R)
#-----------------------------------------------------------------------
plotRanges <- function(x, xlim = x, main = deparse(substitute(x)), col = "black", sep = 0.5, ...) {
  height <- 1
  if (is(xlim, "Ranges"))
    xlim <- c(min(start(xlim)), max(end(xlim)))
  bins <- disjointBins(IRanges(start(x), end(x) + 1))
  plot.new()
  plot.window(xlim, c(0, max(bins)*(height + sep)))
  ybottom <- bins * (sep + height) - height
  rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col = col, ...)
  title(main)
  axis(1)
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
dotplotWithDensityLogXY <- function(plotDF, xvar, yvar, zvar, COL=c("orange", "purple"), ALPHA=.5){

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
    scatter <- ggplot(plotDF,aes_string(xvar, yvar)) + 
      geom_point(aes_string(color=zvar), alpha=ALPHA) + scale_x_log10() + scale_y_log10() + 
      scale_color_manual(values = COL) + theme_bw() + 
      theme(legend.position=c(1,1),legend.justification=c(-.1,-.1), plot.margin =rep(grid::unit(c(0,0,1,1), "cm"),4))  
    
    #marginal density of x - plot on top
    plot_top <- ggplot(plotDF, aes_string(xvar, fill=zvar)) + 
      geom_density(alpha=.5) + scale_x_log10() + theme_bw() + 
      scale_fill_manual(values = COL) + 
      theme(legend.position = "none", axis.title.x = element_blank(), plot.margin = rep(grid::unit(c(1,0,0,1), "cm"),4))
      
    #marginal density of y - plot on the right
    plot_right <- ggplot(plotDF, aes_string(yvar, fill=zvar)) + 
      geom_density(alpha=.5) + scale_x_log10() + coord_flip() + theme_bw() + 
      scale_fill_manual(values = COL) + 
      theme(legend.position = "none", axis.title.y = element_blank(), plot.margin = rep(grid::unit(c(0,1,1,0), "cm"),4))
    
    #arrange the plots together, with appropriate height and width for each row and column
#~     grid.arrange(plot_top, empty, scatter, plot_right, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
    arrangeGrob(plot_top, empty, scatter, plot_right, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
}
#grid::grid.draw(dotplotWithDensityLogXY(plotDFcloseHiC, "dist", "HiCraw", "group", COL))

#-----------------------------------------------------------------------
# save the return value of arrangeGrob 
#-----------------------------------------------------------------------
grobSave <- ggplot2::ggsave; body(ggsave) <- body(ggplot2::ggsave)[-2]

#-----------------------------------------------------------------------
# A very simple heatmap plot of matirx like objects
#-----------------------------------------------------------------------
simple.heatmap <- function (x, Rowv=FALSE, Colv=FALSE, 
        revC=FALSE, dendrogram="none", trace="none", key=TRUE, ClabSide=3, 
        RlabSide=2, margins=c(1,1), col=colorRampPalette(brewer.pal(9,"Blues")),
          ...) {
    my.heatmap.2(
        x, Rowv=Rowv, Colv=Colv, revC=revC, dendrogram=dendrogram, 
        trace=trace, key=key, ClabSide=ClabSide, RlabSide=RlabSide, 
        margins=margins, col=col, ...)
}
 

#-----------------------------------------------------------------------
# Custom changes to the heatmap.2 function from the 
#-----------------------------------------------------------------------
my.heatmap.2 <- function (x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE, 
    distfun = dist, hclustfun = hclust, dendrogram = c("both", 
        "row", "column", "none"), reorderfun = function(d, w) reorder(d, 
        w), symm = FALSE, scale = c("none", "row", "column"), 
    na.rm = TRUE, revC = identical(Colv, "Rowv"), add.expr, breaks, 
    symbreaks = any(x < 0, na.rm = TRUE) || scale != "none", 
    col = "heat.colors", colsep, rowsep, sepcolor = "white", 
    sepwidth = c(0.05, 0.05), cellnote, notecex = 1, notecol = "cyan", 
    na.color = par("bg"), trace = c("column", "row", "both", 
        "none"), tracecol = "cyan", hline = median(breaks), vline = median(breaks), 
    linecol = tracecol, margins = c(5, 5), ColSideColors, RowSideColors, 
    cexRow = 0.2 + 1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL, 
    labCol = NULL, srtRow = NULL, srtCol = NULL, adjRow = c(0, 
        NA), adjCol = c(NA, 0), offsetRow = 0.5, offsetCol = 0.5, 
    key = TRUE, keysize = 1.5, density.info = c("histogram", 
        "density", "none"), denscol = tracecol, symkey = any(x < 
        0, na.rm = TRUE) || symbreaks, densadj = 0.25, key.title = NULL, 
    key.xlab = NULL, key.ylab = NULL, key.xtickfun = NULL, key.ytickfun = NULL, 
    key.par = list(), main = NULL, xlab = NULL, ylab = NULL, 
    lmat = NULL, lhei = NULL, lwid = NULL, extrafun = NULL, ClabColor = "black", RlabColor = "black", ClabSide=1, RlabSide=4, ...) 
{
    scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
    }
    retval <- list()
    scale <- if (symm && missing(scale)) 
        "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col)) 
        col <- get(col, mode = "function")
    if (!missing(breaks) && (scale != "none")) 
        warning("Using scale=\"row\" or scale=\"column\" when breaks are", 
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || is.na(Rowv)) 
        Rowv <- FALSE
    if (is.null(Colv) || is.na(Colv)) 
        Colv <- FALSE
    else if (all(Colv == "Rowv")) 
        Colv <- Rowv
    if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
        stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1) 
        stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2) 
        stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote)) 
        cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
        if (((is.logical(Rowv) && !isTRUE(Rowv)) || (is.null(Rowv))) && 
            (dendrogram %in% c("both", "row"))) {
            if (is.logical(Colv) && (Colv)) 
                dendrogram <- "column"
            else dendrogram <- "none"
            warning("Discrepancy: Rowv is FALSE, while dendrogram is `", 
                dendrogram, "'. Omitting row dendogram.")
        }
    }
    if (!inherits(Colv, "dendrogram")) {
        if (((is.logical(Colv) && !isTRUE(Colv)) || (is.null(Colv))) && 
            (dendrogram %in% c("both", "column"))) {
            if (is.logical(Rowv) && (Rowv)) 
                dendrogram <- "row"
            else dendrogram <- "none"
            warning("Discrepancy: Colv is FALSE, while dendrogram is `", 
                dendrogram, "'. Omitting column dendogram.")
        }
    }
    if (inherits(Rowv, "dendrogram")) {
        ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
        if (length(rowInd) > nr || any(rowInd < 1 | rowInd > 
            nr)) 
            stop("Rowv dendrogram doesn't match size of x")
    }
    else if (is.integer(Rowv)) {
        browser()
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorderfun(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd)) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
        Rowv <- rowMeans(x, na.rm = na.rm)
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        ddr <- reorderfun(ddr, Rowv)
        rowInd <- order.dendrogram(ddr)
        if (nr != length(rowInd)) 
            stop("row dendrogram ordering gave index of wrong length")
    }
    else {
        rowInd <- nr:1
    }
    if (inherits(Colv, "dendrogram")) {
        ddc <- Colv
        colInd <- order.dendrogram(ddc)
        if (length(colInd) > nc || any(colInd < 1 | colInd > 
            nc)) 
            stop("Colv dendrogram doesn't match size of x")
    }
    else if (identical(Colv, "Rowv")) {
        if (nr != nc) 
            stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        if (exists("ddr")) {
            ddc <- ddr
            colInd <- order.dendrogram(ddc)
        }
        else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
        hcc <- hclustfun(distfun(if (symm) 
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorderfun(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd)) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
        Colv <- colMeans(x, na.rm = na.rm)
        hcc <- hclustfun(distfun(if (symm) 
            x
        else t(x)))
        ddc <- as.dendrogram(hcc)
        ddc <- reorderfun(ddc, Colv)
        colInd <- order.dendrogram(ddc)
        if (nc != length(colInd)) 
            stop("column dendrogram ordering gave index of wrong length")
    }
    else {
        colInd <- 1:nc
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow)){
        labRow <- if (is.null(rownames(x))) 
            (1:nr)[rowInd]
        else rownames(x)
        if(is.vector(RlabColor)) {
            RlabColor=RlabColor[rowInd]
        }
        if(is.vector(ClabColor)) {
            ClabColor=ClabColor[colInd]
        }
    }else labRow <- labRow[rowInd]
    if (is.null(labCol)) 
        labCol <- if (is.null(colnames(x))) 
            (1:nc)[colInd]
        else colnames(x)
    else labCol <- labCol[colInd]
    if (scale == "row") {
        retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
        x <- sweep(x, 1, rm)
        retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
        x <- sweep(x, 2, rm)
        retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 
        1) {
        if (missing(col) || is.function(col)) 
            breaks <- 16
        else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
        if (!symbreaks) 
            breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm), 
                length = breaks)
        else {
            extreme <- max(abs(x), na.rm = TRUE)
            breaks <- seq(-extreme, extreme, length = breaks)
        }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (class(col) == "function") 
        col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei)) 
        lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid)) 
        lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
        lmat <- rbind(4:3, 2:1)
        if (!missing(ColSideColors)) {
            if (!is.character(ColSideColors) || length(ColSideColors) != 
                nc) 
                stop("'ColSideColors' must be a character vector of length ncol(x)")
            lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 
                1)
            lhei <- c(lhei[1], 0.2, lhei[2])
        }
        if (!missing(RowSideColors)) {
            if (!is.character(RowSideColors) || length(RowSideColors) != 
                nr) 
                stop("'RowSideColors' must be a character vector of length nrow(x)")
            lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 
                1), 1), lmat[, 2] + 1)
            lwid <- c(lwid[1], 0.2, lwid[2])
        }
        lmat[is.na(lmat)] <- 0
    }
    if (length(lhei) != nrow(lmat)) 
        stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat)) 
        stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    if (!missing(RowSideColors)) {
        par(mar = c(margins[1], 0, 0, 0.5))
        image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    }
    if (!missing(ColSideColors)) {
        par(mar = c(0.5, 0, 0, margins[2]))
        image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    }
    par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
        iy <- nr:1
        if (exists("ddr")) 
            ddr <- rev(ddr)
        x <- x[, iy]
        cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
        c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
        breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr")) 
        retval$rowDendrogram <- ddr
    if (exists("ddc")) 
        retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!gtools::invalid(na.color) & any(is.na(x))) {
        mmat <- ifelse(is.na(x), 1, NA)
        image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "", 
            col = na.color, add = TRUE)
    }
    if (is.null(srtCol)){ 
        #(1=bottom, 2=left, 3=top, 4=right)
        mtext(side = ClabSide, text = labCol, at = 1:nc, las = 2, line = 0.5,col = ClabColor, cex = cexCol)

    }else {
        if (is.numeric(srtCol)) {
            if (missing(adjCol) || is.null(adjCol)) 
                adjCol = c(1, NA)
            xpd.orig <- par("xpd")
            par(xpd = NA)
            xpos <- axis(1, 1:nc, labels = rep("", nc), las = 2, 
                tick = 0)
            text(x = xpos, y = par("usr")[3] - (1 + offsetCol) * 
                strheight("M"), labels = labCol, adj = adjCol, 
                cex = cexCol, srt = srtCol)
            par(xpd = xpd.orig)
        }
        else warning("Invalid value for srtCol ignored.")
    }
    if (is.null(srtRow)) {
        mtext(side = RlabSide, text = labRow, at = iy, las = 2, line = 0.5,col = RlabColor, cex = cexCol)
    }
    else {
        if (is.numeric(srtRow)) {
            xpd.orig <- par("xpd")
            par(xpd = NA)
            ypos <- axis(4, iy, labels = rep("", nr), las = 2, 
                line = -0.5, tick = 0)
            text(x = par("usr")[2] + (1 + offsetRow) * strwidth("M"), 
                y = ypos, labels = labRow, adj = adjRow, cex = cexRow, 
                srt = srtRow)
            par(xpd = xpd.orig)
        }
        else warning("Invalid value for srtRow ignored.")
    }
    if (!is.null(xlab)) 
        mtext(xlab, side = 1, line = margins[1] - 1.25)
    if (!is.null(ylab)) 
        mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr)) 
        eval(substitute(add.expr))
    if (!missing(colsep)) 
        for (csep in colsep) rect(xleft = csep + 0.5, ybottom = 0, 
            xright = csep + 0.5 + sepwidth[1], ytop = ncol(x) + 
                1, lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    if (!missing(rowsep)) 
        for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 
            1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 
            1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, 
            col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
        retval$vline <- vline
        vline.vals <- scale01(vline, min.scale, max.scale)
        for (i in colInd) {
            if (!is.null(vline)) {
                abline(v = i - 0.5 + vline.vals, col = linecol, 
                  lty = 2)
            }
            xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
            xv <- c(xv[1], xv)
            yv <- 1:length(xv) - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (trace %in% c("both", "row")) {
        retval$hline <- hline
        hline.vals <- scale01(hline, min.scale, max.scale)
        for (i in rowInd) {
            if (!is.null(hline)) {
                abline(h = i - 0.5 + hline.vals, col = linecol, 
                  lty = 2)
            }
            yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
            yv <- rev(c(yv[1], yv))
            xv <- length(yv):1 - 0.5
            lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
        }
    }
    if (!missing(cellnote)) 
        text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote), 
            col = notecol, cex = notecex)
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
        flag <- try(stats:::plot.dendrogram(ddr, horiz = TRUE, axes = FALSE, 
            yaxs = "i", leaflab = "none"))
        if ("try-error" %in% class(flag)) {
            cond <- attr(flag, "condition")
            if (!is.null(cond) && conditionMessage(cond) == "evaluation nested too deeply: infinite recursion / options(expressions=)?") 
                stop("Row dendrogram too deeply nested, recursion limit exceeded.  Try increasing option(\"expressions\"=...).")
        }
    }
    else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
        flag <- try(stats:::plot.dendrogram(ddc, axes = FALSE, xaxs = "i", 
            leaflab = "none"))
        if ("try-error" %in% class(flag)) {
            cond <- attr(flag, "condition")
            if (!is.null(cond) && conditionMessage(cond) == "evaluation nested too deeply: infinite recursion / options(expressions=)?") 
                stop("Column dendrogram too deeply nested, recursion limit exceeded.  Try increasing option(\"expressions\"=...).")
        }
    }
    else plot.new()
    if (!is.null(main)) 
        title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
        mar <- c(5, 4, 2, 1)
        if (!is.null(key.xlab) && is.na(key.xlab)) 
            mar[1] <- 2
        if (!is.null(key.ylab) && is.na(key.ylab)) 
            mar[2] <- 2
        if (!is.null(key.title) && is.na(key.title)) 
            mar[3] <- 1
        par(mar = mar, cex = 0.75, mgp = c(2, 1, 0))
        if (length(key.par) > 0) 
            do.call(par, key.par)
        tmpbreaks <- breaks
        if (symkey) {
            max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
            min.raw <- -max.raw
            tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
            tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
        }
        else {
            min.raw <- min(x, na.rm = TRUE)
            max.raw <- max(x, na.rm = TRUE)
        }
        z <- seq(min.raw, max.raw, length = length(col))
        image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks, 
            xaxt = "n", yaxt = "n")
        par(usr = c(0, 1, 0, 1))
        if (is.null(key.xtickfun)) {
            lv <- pretty(breaks)
            xv <- scale01(as.numeric(lv), min.raw, max.raw)
            xargs <- list(at = xv, labels = lv)
        }
        else {
            xargs <- key.xtickfun()
        }
        xargs$side <- 1
        do.call(axis, xargs)
        if (is.null(key.xlab)) {
            if (scale == "row") 
                key.xlab <- "Row Z-Score"
            else if (scale == "column") 
                key.xlab <- "Column Z-Score"
            else key.xlab <- "Value"
        }
        if (!is.na(key.xlab)) {
            mtext(side = 1, key.xlab, line = par("mgp")[1], padj = 0.5)
        }
        if (density.info == "density") {
            dens <- density(x, adjust = densadj, na.rm = TRUE)
            omit <- dens$x < min(breaks) | dens$x > max(breaks)
            dens$x <- dens$x[-omit]
            dens$y <- dens$y[-omit]
            dens$x <- scale01(dens$x, min.raw, max.raw)
            lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol, 
                lwd = 1)
            if (is.null(key.ytickfun)) {
                yargs <- list(at = pretty(dens$y)/max(dens$y) * 
                  0.95, labels = pretty(dens$y))
            }
            else {
                yargs <- key.ytickfun()
            }
            yargs$side <- 2
            do.call(axis, yargs)
            if (is.null(key.title)) 
                key.title <- "Color Key\nand Density Plot"
            if (!is.na(key.title)) 
                title(key.title)
            par(cex = 0.5)
            if (is.null(key.ylab)) 
                key.ylab <- "Density"
            if (!is.na(key.ylab)) 
                mtext(side = 2, key.ylab, line = par("mgp")[1], 
                  padj = 0.5)
        }
        else if (density.info == "histogram") {
            h <- hist(x, plot = FALSE, breaks = breaks)
            hx <- scale01(breaks, min.raw, max.raw)
            hy <- c(h$counts, h$counts[length(h$counts)])
            lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s", 
                col = denscol)
            if (is.null(key.ytickfun)) {
                yargs <- list(at = pretty(hy)/max(hy) * 0.95, 
                  labels = pretty(hy))
            }
            else {
                yargs <- key.ytickfun()
            }
            yargs$side <- 2
            do.call(axis, yargs)
            if (is.null(key.title)) 
                key.title <- "Color Key\nand Histogram"
            if (!is.na(key.title)) 
                title(key.title)
            par(cex = 0.5)
            if (is.null(key.ylab)) 
                key.ylab <- "Count"
            if (!is.na(key.ylab)) 
                mtext(side = 2, key.ylab, line = par("mgp")[1], 
                  padj = 0.5)
        }
        else if (is.null(key.title)) 
            title("Color Key")
        if (trace %in% c("both", "column")) {
            vline.vals <- scale01(vline, min.raw, max.raw)
            if (!is.null(vline)) {
                abline(v = vline.vals, col = linecol, lty = 2)
            }
        }
        if (trace %in% c("both", "row")) {
            hline.vals <- scale01(hline, min.raw, max.raw)
            if (!is.null(hline)) {
                abline(v = hline.vals, col = linecol, lty = 2)
            }
        }
    }
    else plot.new()
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)], 
        high = retval$breaks[-1], color = retval$col)
    if (!is.null(extrafun)) 
        extrafun()
    invisible(retval)
}

