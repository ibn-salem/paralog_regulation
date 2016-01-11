#!/usr/bin/env Rscript
########################################################################
#
# Chip-exo read plotting script.
# This script provies function to analyse ChIP-exo coverage profiles 
# around transcription factor binding sites (TFBS). It has been originally
# developed for the analysis in Starick and Ibn-Salem et al. 2015 Genome Research.
# 
########################################################################

require(gdata)

########################################################################
# some global parameters
########################################################################
# maximal 5'coverage count value for visualization in heatmaps

MAX_COUNT = 100
# pseudocount to 5'coverage values, needed for proper log transformation in the heatmaps
PSEUDOCOUNT = 1
PLOT_FORMATS=c("pdf", "png")

# read argument from commandline:
#~ args=(commandArgs(TRUE))
#~ prefix = args[1]


# converts a sequence into a vector of colors
seq2color <- function(s) {
    #return(sapply(strsplit(s, "")[[1]], switch , 
    return(sapply(s, switch , 
        "A"="green", "a"="green", 
        "C"="blue", "c"="blue",
        "G"="gold", "g"="gold",
        "T"="red", "t"="red", 
        "N"="gray", "n"="gray", 
        "darkgray"))
}

equal_dist = function(a,b,m) { return (1 -  sum(a==b)/m ) }

custom.dist <- function(input.matrix, my.function) {
    n <- nrow(input.matrix)
    m <- ncol(input.matrix)
    mat <- matrix(0, ncol = n, nrow = n)
    colnames(mat) <- rownames(mat) <- rownames(input.matrix)
    for(i in 1:nrow(mat)) {
        for(j in 1:ncol(mat)) {
            mat[i,j] <- my.function(input.matrix[i,],input.matrix[j,], m)
    }}
    return(as.dist(mat))
}

parse_profile <- function(prefix){

    # read 5'coverage count data
    upDF = read.table(paste(prefix, ".up_counts.tab", sep=""), header=FALSE)
    up_counts <- as.matrix(upDF[,2:ncol(upDF)])
    
    downDF = read.table(paste(prefix, ".down_counts.tab", sep=""), header=FALSE)
    down_counts <- as.matrix(downDF[,2:ncol(downDF)])
    
    
    # read consensus sequence and convert string in vector of chars
    consensus_str = readLines(paste(prefix, ".consensus.txt", sep=""))
    consensus = substring(consensus_str, seq(1,nchar(consensus_str),1), seq(1,nchar(consensus_str),1))
    
    # get dimensions
    n = dim(up_counts)[1] # number of regions/sites
    m = dim(up_counts)[2] # number of positions
    
    # x-axis values for plotting:
    xs = (-m/2):(m/2-1)
    
    return(list(
        up_counts=up_counts,
        down_counts=down_counts,
        consensus=consensus,
        n=n, m=m, xs=xs
    ))
}

plot_two_footprints <- function(p1, p2, prefix, consensus=NA){


    COL = brewer.pal(12, "Paired")[c(2,6,1,5)]
    IDX=10

#~     up1 = colSums(p1$up_counts[IDX,])
#~     up2 = colSums(p2$up_counts[IDX,])
#~     down1 = colSums(p1$down_counts[IDX,])
#~     down2 = colSums(p2$down_counts[IDX,])

    up1 = p1$up_counts[IDX,]
    up2 = p2$up_counts[IDX,]
    down1 = p1$down_counts[IDX,]
    down2 = p2$down_counts[IDX,]
    
    n = ifelse(is.matrix(up1), nrow(up1), 1)
    m = ifelse(is.matrix(up1), ncol(up1), length(up1))
    xs = (-m/2):(m/2-1)

    if (all(is.na(consensus))) {
        consensus = paste(rep(" ", m), collapse="")
    }
    
    
    ymax = max(c(up1, up2, down1, down2)) 
    DELTA = .5*ymax
    ymax = ymax + DELTA
    
#~     for ( f in PLOT_FORMATS ) {
#~         get(f)(file=paste(prefix,".profile.", IDX, ".", f, sep=""))

            plot(0,0,type="n", xlim=c(-m/2, m/2),ylim=c(-ymax,ymax),
                xlab="Distance from motif center (bp)", xaxt="n", 
                ylab=paste("5' coverage (", n, "sites )"),
                main=paste("Site:", IDX, ", File:\n", basename(prefix)))
            axis(1, line=1)
            #mtext(paste("sum 5' coverage (", n, "sites )"),2,line=2)
            mtext(strsplit(consensus, ""),side = 1,line = ,at = xs, col=seq2color(consensus), cex=0.8)
    
            # add grid
            abline(v=xs, col="gray", lty="dotted")
            nlines <- ceiling(m / 5 / 2)
            abline(v=-nlines:nlines * 5, col="black", lty=2)
    
            points(xs, up1, col=COL[1],type="h",lwd=3, cex=.6, pch=19)
            points(xs, -1*down1, col=COL[2],type="h",lwd=3, cex=.6, pch=19)

            points(xs, up2+DELTA, col=COL[3],type="o",lwd=3, cex=.6, pch=19)
            points(xs, -(1*down2)+DELTA, col=COL[4],type="o",lwd=3, cex=.6, pch=19)
                
            legend("topright", c("Anc1 for", "Anc1 rev", "Anc2 for", "Anc2 rev"), bg="white", 
                lty=1, pt.cex=.6, pch=19, col=COL, lwd=3)

#~         graphics.off()
#~     }
}

plot_anchor_footprints <- function(p1, p2, outFile){

    n = p1$n
    m = p1$m
    xs = (-1*m/2):(m/2-1)

    plotDF = data.frame(
        pos = rep(rep(xs, 2), 2*n),
        strand = rep(c(rep('+', m), rep('-',m)), 2*n),
        reg=rep(rep(1:n, each=2*m), 2),
        anchor=rep(factor(c("up", "down"), c("up", "down")), each=2*n*m),
        counts = c(
            as.vector(t(cbind(p1$up_counts, -1*p1$down_counts))),
            as.vector(t(cbind(p2$up_counts, -1*p2$down_counts)))
            )
        )

    
#~     plotDF = cbind(plotDF,t(cbind(p1$up_counts, p1$down_counts)))

#~         geom_ribbon(aes(ymin=0, ymax=counts)) + 
#~         geom_bar(stat='identity', aes(y=counts)) + 

    p = ggplot(plotDF, aes(x=pos, fill=strand)) +
        geom_rect(aes(ymin=0, ymax=counts, xmin=pos-.6, xmax=pos+.6)) + 
        facet_grid(anchor ~ reg ) +
        ylim(-50, 50) + theme_bw() + 
        labs(y="ChIP-exo 5' coverage", x="Positions around motif") +
        scale_fill_manual(values=COL[c(2, 1)])

    # save ggplot as pdf
    ggsave(p, file=outFile, w=14, h=3.5)

    
#~         geom_bar(stat='identity') +
#~         geom_bar(stat='identity') + facet_grid(reg~.) 
#~     p
 
#~     scale_color_manual(values=COL, guide_legend(title = "")) +
#~     theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + 
#~     guides(fill=guide_legend(title="")) +
#~     labs(y="Hi-C counts", x="Linear distance bin [kb]")
}

plot_exo_profile <- function(profileList, prefix){
        
    up_counts = profileList[["up_counts"]]
    down_counts = profileList[["down_counts"]]
    consensus = profileList[["consensus"]]
    n = profileList[["n"]]
    m = profileList[["m"]]
    
    ########################################################################
    # 5'coverage profile plot
    ########################################################################
    up_sum = colSums(up_counts) # /n
    down_sum = colSums(down_counts) # / n
    ymax = max(c(up_sum, down_sum))
    
    for ( f in PLOT_FORMATS ) {
        get(f)(file=paste(prefix,".profile.", f, sep=""))
            plot(0,0,type="n", xlim=c(-m/2, m/2),ylim=c(0,ymax),
                xlab="Distance from motif center (bp)", xaxt="n", 
                ylab=paste("5' coverage (", n, "sites )"))
            axis(1, line=1)
            #mtext(paste("sum 5' coverage (", n, "sites )"),2,line=2)
            mtext(strsplit(consensus, ""),side = 1,line = ,at = xs, col=seq2color(consensus), cex=0.8)
    
            # add grid
            abline(v=xs, col="gray", lty="dotted")
            nlines <- ceiling(m / 5 / 2)
            abline(v=-nlines:nlines * 5, col="black", lty=2)
    
            points(xs, up_sum, col="blue",type="o",lwd=3, cex=.6, pch=19)
            points(xs, down_sum, col="red",type="o",lwd=3, cex=.6, pch=19)
                
            legend("topright", c("Forward", "Reverse"), bg="white", 
                lty=c(1,1), pt.cex=.6, pch=19, col=c("blue", "red"), lwd=3)
    
        graphics.off()
    }

}

plot_exo_heatmaps <- function(profileList, prefix){
    
    up_counts = profileList[["up_counts"]]
    down_counts = profileList[["down_counts"]]
    consensus = profileList[["consensus"]]
    n = profileList[["n"]]
    m = profileList[["m"]]
    
    
    ########################################################################
    # Heatmap like plot with max of forward and reverse counts per cell:
    ########################################################################
    
    # scale data logarithmically for better visualisation of low count sites.
    #~ heat.val.up = log10(up_counts + PSEUDOCOUNT)/log10(max(up_counts) + PSEUDOCOUNT)
    #~ heat.val.down = log10(down_counts+ PSEUDOCOUNT)/log10(max(down_counts) + PSEUDOCOUNT)
    heat.val.up = log10(up_counts + PSEUDOCOUNT)/log10(MAX_COUNT + PSEUDOCOUNT)
    heat.val.down = log10(down_counts+ PSEUDOCOUNT)/log10(MAX_COUNT + PSEUDOCOUNT)
    
    # take maximum of up and down counts
    heat.val = ifelse(heat.val.up >= heat.val.down, heat.val.up, -heat.val.down )
    
    # force heat map values between -1 and 1
    # In this way, higher values than MAX_COUNT will have the color intensity coresponding to MAX_COUNT
    heat.val[heat.val > 1] = 1
    heat.val[heat.val < -1] = -1
    
    
    for ( f in PLOT_FORMATS ) {
        get(f)(file=paste(prefix,".maxstrand_heatmap.", f, sep=""))
    
            # make tow plots in one figure above each other
            #par(mar= 'c(bottom, left, top, right)', default is c(5,4,4,2)
    #~         def.par <- par(no.readonly = TRUE,mar=c(5, 3, 0, 2))
            def.par <- par(no.readonly = TRUE,mar=c(5, 2, 4, 1))
            nf <- layout(matrix(c(1,2), 1,2),  width = c(.9, .1), T)
            #layout.show(nf)
    
            # plot heatmap
            image(x=xs, z=t(heat.val), xlim=c(-m/2, m/2), zlim=c(-1,1), 
                col=colorRampPalette(c("red","white","blue"))(100), 
                xaxt="n", yaxt="n",xlab="Distance from motif (bp)")
            
            # add axis, axis labels and consensus seq
            axis(1, line=1, xpd=TRUE, xaxs="i")
            mtext(paste("Binding sites ( n =", n, ")"),2,line=1)
            mtext(strsplit(consensus, ""),side = 1,line = 0,at = xs, col=seq2color(consensus), cex=0.8)
            
            #===============================================================
            # Legend with color scales from white to blue/red
            #===============================================================
            r = cbind(seq(0, 1, 0.01), -seq(0, 1, 0.01) )
            # par(mar= 'c(bottom, left, top, right)', default is c(5,4,4,2)
            op = par(mar=c(12, 0, 4, 2) )
            image(x=c(1,2), y=seq(0, 1, 0.01), z=t(r), col=colorRampPalette(c("red","white","blue"))(100), ylab="", xlab="", xaxt="n", yaxt='n')
            par(op)
            
            mtext(c("Forward", "Reverse", "5' coverage"), side=4, line=c(-3, -2, 0), at=-0.075, adj=1)
            
            #pow=c(0,1)
            pow=0:(log10(MAX_COUNT)-1)
            ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
            myticks = log10(ticksat)/log10(max(ticksat))
    
    
            xlabs = 10^(0:log10(MAX_COUNT)) # 1, 10, 100
    
            axis(side=4, at=0:log10(MAX_COUNT)/log10(MAX_COUNT), labels=xlabs, line=-1)
            # add ticks:
            axis(side=4, at=myticks, labels=NA, tcl=-0.25, lwd=0, lwd.ticks=1, line=-1)
            par(def.par)
        dev.off()
    }


    ########################################################################
    # color map of centered sequences and clustering of sites by sequence
    ########################################################################
    # check if sequence matrix file exists:
    if (file.exists(paste(prefix, ".seq_matrix.tab", sep="")) ){
    #~     seq_matrix = read.table(paste(prefix, ".seq_matrix.tab", sep=""), header=FALSE, row.names=1)
        seq_matrix = read.table(paste(prefix, ".seq_matrix.tab", sep=""), header=FALSE)[,2:(m+1)]
        
        # plot color map
        for ( f in PLOT_FORMATS ) {
        get(f)(file=paste(prefix,"_seq.", f, sep=""))
            image(x=xs,y=1:n, z=t(seq_matrix), xlim=c(-m/2, m/2), 
                col=c("green", "blue", "yellow", "red"),
                xaxt="n", yaxt="n",xlab="distance relative to center", ylab="")
            axis(1, line=1)
            mtext(strsplit(consensus, ""),side = 1,line = 0,at = xs, col=seq2color(consensus), cex=0.8)
            mtext(paste("Binding sites ( n =", n, ")"),2,line=1)
            #,xaxt="n",yaxt="n"
        graphics.off()
        }
        
        # remove labels
        seq_matrix = as.matrix(seq_matrix)
        rownames(seq_matrix) = rep("", nrow(seq_matrix))
        # take only a subset of maximal 100 sites
        top.sites = min(n, 200) #500
        
        # coniser only the inner +- BASES nucleotides for the sequence clustering
        BASES=10
        consider = ((m/2)-BASES):((m/2)+BASES-1)
    
        # calcualte paire-wise distances
        d = custom.dist(seq_matrix[1:top.sites,consider], equal_dist)
        # cluster:
        hc = hclust(d)
        den = as.dendrogram(hc)
        
        pdf(file=paste(prefix,".seq_cluster.pdf", sep=""), width=14, height=7)
    #~     pdf(file=paste(prefix,".seq_cluster.pdf", sep=""), width=7, height=7)
        
            def.par <- par(no.readonly = TRUE,mar=c(5,2,4,1)) #c(5, 3, 0, 2))
            nf <- layout(matrix(1:4, 1,4),  width = c(.15, .425, 0.3825, 0.0425), T)
    #~         nf <- layout(matrix(c(1,2), 1,2),  width = c(.3, .7), T)
            dpar = par(mar=c(4, 4, 3, 0) )
            plot(den, horiz=TRUE, cex.axis=1.5)
            mtext("Distance" ,side = 1,line = 2)
            mtext(paste("Clustering by sequence (",length(hc$labels), "sites )"), side=2,line=2)
            
            # par(mar= 'c(bottom, left, top, right)', default is c(5,4,4,2)
            op = par(mar=c(5.5, 1, 4.5, 1) )
            image(x=xs,y=1:top.sites, z=t(seq_matrix[hc$order,]), xlim=c(-m/2, m/2), 
                col=c("green", "blue", "yellow", "red"),
                xaxt="n", yaxt="n", ylab="", xlab="")
            axis(1, line=1.5, xpd=TRUE, xaxs="i", cex.axis=1.5)
            mtext(strsplit(consensus, ""),side = 1,line = .5,at = xs, col=seq2color(consensus), cex=0.8)
            mtext("Distance from motif center (bp)", side=1, line=4)
            par(op)
    
            # par(mar= 'c(bottom, left, top, right)', default is c(5,4,4,2)
            op = par(mar=c(5.5, 1, 4.5, 1) )
            image(x=xs, z=t(heat.val[hc$order,]), xlim=c(-m/2, m/2), 
                col=colorRampPalette(c("red","white","blue"))(100), xaxt="n", yaxt="n", xlab="")
            axis(1, line=1.5, xpd=TRUE, xaxs="i", cex.axis=1.5)
            mtext(strsplit(consensus, ""),side = 1,line = .5,at = xs, col=seq2color(consensus), cex=0.8)
            mtext("Distance from motif center (bp)", side=1, line=4)
            par(op)
           
            r = cbind(seq(0, 1, 0.01), -seq(0, 1, 0.01) )
            # par(mar= 'c(bottom, left, top, right)', default is c(5,4,4,2)
            op = par(mar=c(12, 0, 4.5, 3) )
            image(x=c(1,2), y=seq(0, 1, 0.01), z=t(r), col=colorRampPalette(c("red","white","blue"))(100), ylab="", xlab="", xaxt="n", yaxt='n')
            par(op)
            mtext(c("Forward", "Reverse", "5' coverage"), side=4, line=c(-5, -4, -2), at=-0.075, adj=1)        
            pow=c(0,1)
            ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
            myticks = log10(ticksat)/log10(max(ticksat))
            xlabs = expression(1, 10, 100)
            axis(side=4, at=c(0, .5, 1), labels=xlabs, line=-3, cex.axis=1.5)
            axis(side=4, at=myticks, labels=NA, tcl=-0.25, lwd=0, lwd.ticks=1, line=-3, cex.axis=1.5)        
    
            par(def.par)
            
        dev.off()
    
    }

    ########################################################################
    # write heatmap values (just for transfering data to Paris...)
    ########################################################################
    #~ write.table(heat.val, file=paste(prefix,".heatmap_values.tab", sep=""), quote=FALSE, sep="\t", row.names=TRUE, col.names=FALSE)

    
    
    ########################################################################
    # hierarchical clustering of sites
    ########################################################################
    #~ print("INFO: Begin to clustering analysis.")
    
    top.sites = min(n, 500) #500
    pseudo.count = 0.01
    
    #mat = cbind(log10(up_counts+pseudo.count), log10(down_counts+pseudo.count))
    #mat = cbind(up_counts,down_counts)
    
    if (n>1) {
        maxCount = max(up_counts[1:top.sites,], down_counts[1:top.sites,]) +1
    #~     normedCountsUp = log10(up_counts[1:top.sites,]+1)/log10(maxCount)
    #~     normedCountsDown = log10(down_counts[1:top.sites,]+1)/log10(maxCount)
        normedCountsUp = log10(up_counts[1:top.sites,]+1)/log10(rowSums(up_counts[1:top.sites,]+1))
        normedCountsDown = log10(down_counts[1:top.sites,]+1)/log10(rowSums(down_counts[1:top.sites,]+1))
        mat = cbind(t(apply(normedCountsUp,1, smooth)), t(apply(normedCountsDown,1, smooth)) )
        
        #mat = t(apply(mat[1:top.sites,],1, smooth))
        #mat = cbind(up_counts, down_counts)[1:top.sites,]
        #mat = t(apply(normedCounts[1:top.sites,] / rowSums(normedCounts[1:top.sites,]),1, smooth))
        rownames(mat) = rep("", top.sites)
        
        hc = hclust(dist(mat))
        #hc = hclust(dist(mat), method="ward.D2")
        #hc = hclust(dist(mat, method="mahalanobis"))
        #hc = hclust(mahalanobis(mat))
        #hc = hclust(dist(mat))
        #hc$labels = NULL
        den = as.dendrogram(hc)
        
        pdf(file=paste(prefix,".cluster.pdf", sep=""))
        
            #par(mar= 'c(bottom, left, top, right)', default is c(5,4,4,2)
            def.par <- par(no.readonly = TRUE,mar=c(5,2,4,1)) #c(5, 3, 0, 2))
            nf <- layout(matrix(c(1,2), 1,2),  width = c(.3, .7), T)
        
            dpar = par(mar=c(4, 4, 3, 0) )
            plot(den, horiz=TRUE, 
            ylab=paste("Hierarchical clustering of n =",length(hc$labels), "sites"))
            mtext("Distance" ,side = 1,line = 2)
            par(dpar)
            
            # par(mar= 'c(bottom, left, top, right)', default is c(5,4,4,2)
            op = par(mar=c(5, 1, 4, 2) )
            image(x=xs, z=t(heat.val[hc$order,]), xlim=c(-m/2, m/2), col=colorRampPalette(c("red","white","blue"))(100), xaxt="n", yaxt="n",xlab="Distance from motif (bp)")
            axis(1, line=1, xpd=TRUE, xaxs="i")
            mtext(strsplit(consensus, ""),side = 1,line = 0,at = xs, col=seq2color(consensus), cex=0.8)
            par(op)
            
            par(def.par)
            
        dev.off()
    }
    #plot(as.dendrogram(hclust(dist(m)) leaflab="none"), horiz=TRUE)


    ########################################################################
    # K-means clustering
    ########################################################################
    # number of clusters for K-means clustering
    K=4 
    
    #library(RColorBrewer)
    #cols_paired = brewer.pal(2*K, "Paired")
    cols_paired = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00")
    cols_up = cols_paired[seq(1,2*K, 2)]
    cols_down = cols_paired[seq(2,2*K+1, 2)]
    
    # run clustering on normalized and smoothed counts
    km = kmeans(mat, K)
    
    # get cluster centers
    centers_up = km$centers[,seq(1,m)]
    centers_down = km$centers[,seq(m+1, 2*m)]
    
    # plot the profile of all centers
    pdf(file=paste(prefix,".kmeans_centers.pdf", sep=""))
    
        ymax = max(km$centers)
        plot(0,0,type="n", xlim=c(-m/2, m/2),ylim=c(0,ymax),
            xlab="Distance from motif center (bp)", xaxt="n", 
            ylab=paste("Normalized and smoothed 5' coverage"), main="K-means cluster centers")
        axis(1, line=1)
        mtext(strsplit(consensus, ""),side = 1,line = ,at = xs, col=seq2color(consensus), cex=0.8)
    
        # add grid
        abline(v=xs, col="gray", lty="dotted")
        nlines <- ceiling(m / 5 / 2)
        abline(v=-nlines:nlines * 5, col="black", lty=2)
    
        # plot all k centers separately on the same plot
        for ( k in 1:K) {
            points(xs, centers_up[k,], col=cols_up[k], type="o",lwd=3, cex=.6, pch=19)
            points(xs, centers_down[k,], col=cols_down[k], type="o",lwd=3, cex=.6, pch=19)
        }
        
        legend_labs = paste("cluster", rep(1:k, each=2), rep(c("For", "Rev"), K))
        legend("topright", legend_labs, bg="white", 
            lty=c(1,1), pt.cex=.6, pch=19, col=cols_paired, lwd=3)
    
    dev.off()
    
    # get mean profile for each cluster:
    cluster_mean_up = rbind()
    cluster_mean_down = rbind()
    for (k in 1:K){
        cluster_mean_up = rbind(cluster_mean_up, colMeans(up_counts[1:top.sites,][km$cluster == k,]))
        cluster_mean_down = rbind(cluster_mean_down, colMeans(down_counts[1:top.sites,][km$cluster == k,]))
    }
    
    # plot the profile of all clusters
    pdf(file=paste(prefix,".kmeans_cluster_profile.pdf", sep=""))
    
        ymax = max(cluster_mean_up, cluster_mean_down)
        plot(0,0,type="n", xlim=c(-m/2, m/2),ylim=c(0,ymax),
            xlab="Distance from motif center (bp)", xaxt="n", 
            ylab=paste("Average 5' coverage"),
            main=paste("K-means clustering of the top n =", top.sites, "sites"))
        axis(1, line=1)
        mtext(strsplit(consensus, ""),side = 1,line = ,at = xs, col=seq2color(consensus), cex=0.8)
    
        # add grid
        abline(v=xs, col="gray", lty="dotted")
        nlines <- ceiling(m / 5 / 2)
        abline(v=-nlines:nlines * 5, col="black", lty=2)
    
        # plot all k centers separately on the same plot
        for ( k in 1:K) {
            points(xs, cluster_mean_up[k,], col=cols_up[k], type="o",lwd=3, cex=.6, pch=19)
            points(xs, cluster_mean_down[k,], col=cols_down[k], type="o",lwd=3, cex=.6, pch=19)
        }
        
        legend_labs = paste("cluster", rep(1:k, each=2), rep(c("For", "Rev"), K))
        legend("topright", legend_labs, bg="white", 
            lty=c(1,1), pt.cex=.6, pch=19, col=cols_paired, lwd=3)
    
    dev.off()
    
    # plot the heat map ordered by clusters
    pdf(file=paste(prefix,".kmeans_cluster_heatmap.pdf", sep=""))
    
            # make tow plots in one figure above each other
            def.par <- par(no.readonly = TRUE,mar=c(5, 2, 4, 1))
            nf <- layout(matrix(c(1,2), 1,2),  width = c(.9, .1), T)
    
            # plot heatmap
            image(x=xs, z=t(heat.val[order(km$cluster, decreasing=TRUE),]), xlim=c(-m/2, m/2), zlim=c(-1,1), 
                col=colorRampPalette(c("red","white","blue"))(100), 
                xaxt="n", yaxt="n",xlab="Distance from motif (bp)",
                main=paste("K-means clustering of the top n =", top.sites, "sites"))
    
            # add cluster annotation
            rasterImage(as.raster(cols_down)[sort(km$cluster)], -31, 0, -29, 1, interpolate=FALSE)
            
            # add annotation to clusters
            n_cluster = table(sort(sort(km$cluster)))
            verticla_midpoints = sapply(1:K, function(k){mean(seq(1, top.sites)[sort(km$cluster, decreasing=TRUE) == k])}) / top.sites
            mtext(paste("cluster ", 1:K, "\nn=", n_cluster, sep=""), side=2, at=verticla_midpoints)
            
            # add axis, axis labels and consensus seq
            axis(1, line=1, xpd=TRUE, xaxs="i")
            mtext(strsplit(consensus, ""),side = 1,line = 0,at = xs, col=seq2color(consensus), cex=0.8)
            
            #===============================================================
            # Legend with color scales from white to blue/red
            #===============================================================
            r = cbind(seq(0, 1, 0.01), -seq(0, 1, 0.01) )
            # par(mar= 'c(bottom, left, top, right)', default is c(5,4,4,2)
            op = par(mar=c(12, 0, 4, 2) )
            image(x=c(1,2), y=seq(0, 1, 0.01), z=t(r), col=colorRampPalette(c("red","white","blue"))(100), ylab="", xlab="", xaxt="n", yaxt='n')
            par(op)
            
            mtext(c("Forward", "Reverse", "5' coverage"), side=4, line=c(-3, -2, 0), at=-0.075, adj=1)
            
            pow=0:(log10(MAX_COUNT)-1)
            ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
            myticks = log10(ticksat)/log10(max(ticksat))
            xlabs = 10^(0:log10(MAX_COUNT)) # 1, 10, 100
    
            axis(side=4, at=0:log10(MAX_COUNT)/log10(MAX_COUNT), labels=xlabs, line=-1)
            # add ticks:
            axis(side=4, at=myticks, labels=NA, tcl=-0.25, lwd=0, lwd.ticks=1, line=-1)
    
            #axis(1, at=10^c(0,2,4,6), labels=expression(1, 10^2,10^4, 10^6))
            par(def.par)
        
    dev.off()
    
    getConsensus <- function(seq_matrix){
        freq = apply(seq_matrix/nrow(seq_matrix), 2, table) 
        max_freq = lapply(counts, which.max)
    }
    
    # write count tables for each cluster
    for (k in 1:K) {
        write.table(up_counts[1:top.sites,][km$cluster == k,], 
            file=paste(prefix,".kmeans_cluster_", k, ".up_counts.tab", sep=""),
            quote=FALSE, sep="\t", row.names=TRUE, col.names=FALSE)
        write.table(down_counts[1:top.sites,][km$cluster == k,], 
            file=paste(prefix,".kmeans_cluster_", k, ".down_counts.tab", sep=""),
            quote=FALSE, sep="\t", row.names=TRUE, col.names=FALSE)
    
        # TODO: write consensus sequence for clusters and write BED files
        #if (file.exists(paste(prefix, ".seq_matrix.tab", sep="")) ){
        #   subset_consensus = getConsensus(seq_matrix[1:top.sites,][km$cluster == k,])
        #}
    }
}

# TODO: use different weight function:
# - see: http://stats.stackexchange.com/questions/15289/when-to-use-weighted-euclidean-distance-and-how-to-determine-the-weights-to-use
# - http://en.wikipedia.org/wiki/Mahalanobis_distance

########################################################################

