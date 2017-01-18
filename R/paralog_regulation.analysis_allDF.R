########################################################################
#
# A script to analyse the co-regulation of paralog genes 
#
########################################################################

require(stringr)        # for some string functionality
require(RColorBrewer)   # for nice colors
require(colorspace)     # for some more colors
require(GenomicRanges)  # for genomic intervals and overlap calculation
require(plyr)           # count() function
require(ggplot2)        # for nice plots
require(scales)         # for proper logarithmic scales in ggplot

#-----------------------------------------------------------------------
# Load parameters from external script
#-----------------------------------------------------------------------
# read parameter script location from command-line argument
#~ args <- commandArgs(trailingOnly = TRUE)
#~ PARAM_SCRIPT=args[1]

PARAM_SCRIPT="R/paralog_regulation.param.v16.R"
source(PARAM_SCRIPT)

#-----------------------------------------------------------------------
# Options for parallel computation
#-----------------------------------------------------------------------
# use all available cores but generate random number streams on each worker
multicorParam <- MulticoreParam(RNGseed=RANDOM_SEED)   
# set options
register(multicorParam)  
# bpparam() # to print current options


HiCcolumns = c(c("HiCRaw", "HiC", "HiCobsExp", "captureC_raw", "captureC_ObsExp"), paste0(c("HiCRaw", "HiC", "HiCobsExp", "captureC_raw", "captureC_ObsExp"), "NoZero"))
HiClabels = c(
HiCRaw="Hi-C contacts",
HiC="Normalized Hi-C contacts",
HiCobsExp="Normalized Hi-C [obs/exp]",
captureC_raw="Capture Hi-C contacts",
captureC_ObsExp="Normalized Capture Hi-C [obs/exp]",
HiCRawNoZero="Hi-C contacts",
HiCNoZero="Normalized Hi-C contacts",
HiCobsExpNoZero="Normalized Hi-C [obs/exp]",
captureC_rawNoZero="Capture Hi-C contacts",
captureC_ObsExpNoZero="Normalized Capture Hi-C [obs/exp]"
)

subTADlevels <- c("no TAD", "diff TAD", "diff sub TAD", "same sub TAD")
subTADlabels <- paste0("\n\n\n\n", subTADlevels)
names(subTADlabels) <- subTADlevels

#-----------------------------------------------------------------------
# load some custom functions
#-----------------------------------------------------------------------
source("R/functions.plot.R")
source("R/functions.GRanges.R")
source("R/functions.genePairs.R")
source("R/functions.genePairs.randomization.R")
source("R/functions.genePairs.paralog_analysis.R")
source("R/parseHiC.R")


#=======================================================================
# 1.) Load data from exported data.frames
#=======================================================================

#~ INallDF <- read.delim(paste0(outPrefix, ".allDF.csv"), header=TRUE, sep="\t")

# load aP and allDF
message(paste("INFO: Start loading data from file:", paste0(WORKIMAGE_FILE, ".aP_allDF.Rdata")))
load(paste0(WORKIMAGE_FILE, ".aP_allDF.Rdata"))

# define some often used subsetting vectors for the allDF
sSamp <- allDF$sampType %in% c("paralogs", "sampled by distance")
sSampEh <- allDF$sampType %in% c("paralogs", "sampled by distance\nand enhancer")
sTAD <- allDF$tadSource == "Rao_IMR90"
sExp <- allDF$expSource == "GTEx"
sSpec <- allDF$species == "mouse"

#=======================================================================
# 2.) Enhancer sharing for different sampling approaches
#=======================================================================

for (sampName in levels(allDF$sampType)){
    
    #select subset based on samplingType
	sampNameStr <- gsub("\n| ", "_", sampName)
    
    subSet <- sTAD & sExp &sSpec & allDF$sampType %in% c("paralogs", sampName)
    ehDF <- allDF[subSet,]
    
    
    # get subset for the size 10-1000kb
    subDF <- ehDF[ehDF$distTadBin == "10-1000kb",]
    
    #-------------------------------------------------------------------
    # Shared Enhancer vs. groups
    #-------------------------------------------------------------------
    dc <- ddply(subDF, .(group, replicate), summarize, n=length(eh), count=sum(eh), percent=sum(eh)/length(eh)*100)
    freqDF <- ddply(dc, .(group), summarize, 
        avgCount=mean(count, na.rm=TRUE), 
        sdCount=sd(count, na.rm=TRUE), 
        avgPercent=mean(percent), 
        sdPercent=sd(percent), 
        avgN=mean(n), 
        sdN=sd(n)
    )
    
    # calculate p-values
    pvals <- fisher.test(subDF$eh, subDF$group)$p.value 
    
    pvalDF <- data.frame(
        p = pvals,
        ypos = max(freqDF$avgCount)+10,
        yposPercent = max(freqDF$avgPercent)+2,
        group=NA
    )
    
    p <- ggplot(freqDF, aes(x=group, y=avgCount, fill=group)) +
        geom_errorbar(aes(ymax = avgCount + sdCount , ymin=avgCount - sdCount), width=.25) +
        geom_bar(stat="identity", colour="black") + 
        scale_fill_manual(values=COL, guide_legend(title = "")) +
        theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="none") + labs(y="Gene pairs with shared enhancer", x="") + 
        geom_text(aes(label=avgCount, y= avgCount), vjust=1.25, size=5) +
        geom_text(aes(label=paste0("p=",signif(p,2)), y=ypos, x=1.5), size=5, data=pvalDF)
    ggsave(paste0(outPrefix, ".", sampNameStr, ".10_1000kb.eh.by_group.barplot.pdf"), w=3.5)
    
    p <- ggplot(freqDF, aes(x=group, y=avgPercent, fill=group)) +
        geom_errorbar(aes(ymax = avgPercent + sdPercent , ymin=avgPercent - sdPercent), width=.25) +
        geom_bar(stat="identity", colour="black") + 
        scale_fill_manual(values=COL, guide_legend(title = "")) +
        theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="none") + labs(y="Percent of pairs with shared enhancer  [%]", x="") + 
        geom_text(aes(label=signif(avgPercent,3), y= avgPercent), vjust=1.25, size=5) +
        geom_text(aes(label=paste0("p=",signif(p,2)), y=yposPercent, x=1.5), size=5, data=pvalDF)
    ggsave(paste0(outPrefix, ".", sampNameStr, ".10_1000kb.ehPercent.by_group.barplot.pdf"), w=3.5)
    
    #-------------------------------------------------------------------
    # Shared Enhancer vs. groups and subTad structure
    #-------------------------------------------------------------------
    dc <- ddply(subDF, .(group, subTAD, replicate), summarize, n=length(eh), count=sum(eh), percent=sum(eh)/length(eh)*100)
    freqDF <- ddply(dc, .(group, subTAD), summarize, 
        avgCount=mean(count, na.rm=TRUE), 
        sdCount=sd(count, na.rm=TRUE), 
        avgPercent=mean(percent), 
        sdPercent=sd(percent), 
        avgN=mean(n), 
        sdN=sd(n)
    )
    
    # calculate p-values
    pvals <- sapply(levels(subDF$subTAD), function(sT) fisher.test(table(subDF[subDF$subTAD == sT,"eh"], subDF[subDF$subTAD == sT,"group"]))$p.value) 
    
    pvalDF <- data.frame(
        subTAD=names(pvals),
        p = pvals,
        ypos = ddply(freqDF, .(subTAD), summarize, y=max(avgCount)+10)$y,
        yposPercent = ddply(freqDF, .(subTAD), summarize, y=max(avgPercent)+2)$y,
        group=NA
    )
    
    p <- ggplot(freqDF, aes(x=group, y=avgCount, fill=group)) +
        geom_bar(stat="identity", colour="black") + 
        geom_errorbar(aes(ymax = avgCount + sdCount , ymin=avgCount - sdCount), width=.25) +
        facet_grid(.~subTAD) +
        scale_fill_manual(values=COL, guide_legend(title = "")) +
        theme_bw() + theme(text = element_text(size=20), axis.text.x=element_blank(), legend.position="bottom") + labs(y="Gene pairs with shared enhancer", x="") + theme(strip.text.x = element_text(size=4, angle=90)) +
        geom_text(aes(label=avgCount, y= avgCount), vjust=1.25, size=5) +
        geom_text(aes(label=paste0("p=",signif(p,2)), y=ypos, x=1.5), size=5, data=pvalDF)
    
    g <- addPictureLabels(p, subTADfigPaths)
    pdf(paste0(outPrefix, ".", sampNameStr, ".10_1000kb.eh.by_subTAD_and_group.barplot.pdf"))
        grid.draw(g)
    dev.off()
    
    # same with percent values
    p <- ggplot(freqDF, aes(x=group, y=avgPercent, fill=group)) +
        geom_bar(stat="identity", colour="black") + 
        geom_errorbar(aes(ymax = avgPercent + sdPercent , ymin=avgPercent - sdPercent), width=.25) +
        facet_grid(.~subTAD) +
        scale_fill_manual(values=COL, guide_legend(title = "")) +
        theme_bw() + theme(text = element_text(size=20), axis.text.x=element_blank(), legend.position="bottom") + labs(y="Percent of pairs with shared enhancer  [%]", x="") + theme(strip.text.x = element_text(size=4, angle=90)) +
        geom_text(aes(label=signif(avgPercent,3), y= avgPercent), vjust=1.25, size=5) +
        geom_text(aes(label=paste0("p=",signif(p,2)), y=yposPercent, x=1.5), size=5, data=pvalDF)
    
    g <- addPictureLabels(p, subTADfigPaths)
    pdf(paste0(outPrefix, ".", sampNameStr, ".10_1000kb.ehPercent.by_subTAD_and_group.barplot.pdf"))
        grid.draw(g)
    dev.off()

}

#-----------------------------------------------------------------------
# plot all shared enhancer by sampling type
#-----------------------------------------------------------------------
subSet <- sTAD & sExp &sSpec
ehDF <- allDF[subSet,]
# get subset for the size 10-1000kb
subDF <- ehDF[ehDF$distTadBin == "10-1000kb",]

subDF <- subset(allDF, sTAD & sExp & sSpec & distTadBin == "10-1000kb")

#-------------------------------------------------------------------
# Shared Enhancer vs. samp Groups
#-------------------------------------------------------------------
dc <- ddply(subDF, .(sampType, replicate), summarize, n=length(eh), count=sum(eh), percent=sum(eh)/length(eh)*100)
freqDF <- ddply(dc, .(sampType), summarize, 
    avgCount=mean(count, na.rm=TRUE), 
    sdCount=sd(count, na.rm=TRUE), 
    avgPercent=mean(percent), 
    sdPercent=sd(percent), 
    avgN=mean(n), 
    sdN=sd(n)
)

# calculate p-values
pvals <- sapply(levels(subDF$sampType)[-1], function(st) fisher.test(subDF[subDF$sampType %in% c("paralogs", st),"eh"], subDF[subDF$sampType %in% c("paralogs", st), "sampType"])$p.value )


pvalDF <- data.frame(
    p = pvals,
    ypos = max(freqDF$avgCount)+10,
    yposPercent = max(freqDF$avgPercent)+2,
    sampType=levels(subDF$sampType)[-1]
)

p <- ggplot(freqDF, aes(x=sampType, y=avgCount, fill=sampType)) +
    geom_errorbar(aes(ymax = avgCount + sdCount , ymin=avgCount - sdCount), width=.25) +
    geom_bar(stat="identity", colour="black") + 
    scale_fill_manual(values=COL[c(1, rep(2, length(levels(subDF$sampType)[-1])))], guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="none") + labs(y="Gene pairs with shared enhancer", x="") + 
    geom_text(aes(label=avgCount, y= avgCount), vjust=1.25, size=5) +
    geom_text(aes(label=paste0("p=",signif(p,2)), y=ypos), size=5, data=pvalDF)
ggsave(paste0(outPrefix, ".10_1000kb.eh.by_sampType_and_group.barplot.pdf"))

p <- ggplot(freqDF, aes(x=sampType, y=avgPercent, fill=sampType)) +
    geom_errorbar(aes(ymax = avgPercent + sdPercent , ymin=avgPercent - sdPercent), width=.25) +
    geom_bar(stat="identity", colour="black") + 
    scale_fill_manual(values=COL[c(1, rep(2, length(levels(subDF$sampType)[-1])))], guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="none") + labs(y="Pairs with shared enhancer  [%]", x="") + 
    geom_text(aes(label=signif(avgPercent,3), y= avgPercent), vjust=1.25, size=5) +
    geom_text(aes(label=paste0("p=",signif(p,2)), y=yposPercent), size=5, data=pvalDF)
ggsave(paste0(outPrefix, ".10_1000kb.ehPercent.by_sampType_and_group.barplot.pdf"))


#-------------------------------------------------------------------
# Shared Enhancer vs. samp Groups for close paralogs
#-------------------------------------------------------------------

subDF <- subset(allDF, sTAD & sExp & sSpec & distGroup=="close")

dc <- ddply(subDF, .(sampType, replicate), summarize, n=length(eh), count=sum(eh), percent=sum(eh)/length(eh)*100)
freqDF <- ddply(dc, .(sampType), summarize, 
    avgCount=mean(count, na.rm=TRUE), 
    sdCount=sd(count, na.rm=TRUE), 
    avgPercent=mean(percent), 
    sdPercent=sd(percent), 
    avgN=mean(n), 
    sdN=sd(n)
)

# calculate p-values
pvals <- sapply(levels(subDF$sampType)[-1], function(st) fisher.test(subDF[subDF$sampType %in% c("paralogs", st),"eh"], subDF[subDF$sampType %in% c("paralogs", st), "sampType"])$p.value )

pvalDF <- data.frame(
    p = pvals,
    ypos = max(freqDF$avgCount)+10,
    yposPercent = max(freqDF$avgPercent)+2,
    sampType=levels(subDF$sampType)[-1]
)

p <- ggplot(freqDF, aes(x=sampType, y=avgCount, fill=sampType)) +
    geom_errorbar(aes(ymax = avgCount + sdCount , ymin=avgCount - sdCount), width=.25) +
    geom_bar(stat="identity", colour="black") + 
    scale_fill_manual(values=COL[c(1, rep(2, length(levels(subDF$sampType)[-1])))], guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="none") + labs(y="Gene pairs with shared enhancer", x="") + 
    geom_text(aes(label=avgCount, y= avgCount), vjust=1.25, size=5) +
    geom_text(aes(label=paste0("p=",signif(p,2)), y=ypos), size=5, data=pvalDF)
ggsave(paste0(outPrefix, ".close.eh.by_sampType_and_group.barplot.pdf"))

p <- ggplot(freqDF, aes(x=sampType, y=avgPercent, fill=sampType)) +
    geom_errorbar(aes(ymax = avgPercent + sdPercent , ymin=avgPercent - sdPercent), width=.25) +
    geom_bar(stat="identity", colour="black") + 
    scale_fill_manual(values=COL[c(1, rep(2, length(levels(subDF$sampType)[-1])))], guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="none") + labs(y="Pairs with shared enhancer  [%]", x="") + 
    geom_text(aes(label=signif(avgPercent,3), y= avgPercent), vjust=1.25, size=5) +
    geom_text(aes(label=paste0("p=",signif(p,2)), y=yposPercent), size=5, data=pvalDF)
ggsave(paste0(outPrefix, ".close.ehPercent.by_sampType_and_group.barplot.pdf"))


#=======================================================================
# Distance Bins and sub-TAD structure
#=======================================================================

# reduce to sampled by distance and only one source for TAD/exp/species
subDF <- allDF[sSamp & sTAD & sExp & sSpec,]

#-----------------------------------------------------------------------
# barplot percent of pairs by group, inTAD, distTadBin
#-----------------------------------------------------------------------

# count gene pairs by combinations
dc <- ddply(subDF, .(group, replicate, distTadBin, inTAD), summarize, count=length(inTAD), .drop=FALSE)
# flag get existing replicate/group combination
exitingComb <- dc$replicate <= 1 | dc$group == "sampled"
# sum up counts for sub-combination to get percentage values
summedCounts <- rep(ddply(dc, .(group, replicate, distTadBin), summarize, s=sum(count))$s, each=length(levels(dc$inTAD)))
# add percentage values and counts by marking non-existing combinations with NA
freqRepDF <- cbind(dc, cbind(summedCounts, existingCount=ifelse(exitingComb, dc$count, NA), percent=dc$count * 100 / summedCounts))

# combine replicates by taking average and sd of counts and percentages
freqDF <- ddply(freqRepDF, .(group, distTadBin,inTAD), summarize, avgCount=mean(existingCount, na.rm=TRUE), sdCount=sd(existingCount, na.rm=TRUE), avgPercent=mean(percent, na.rm=TRUE), sdPercent=sd(percent, na.rm=TRUE))

p <- ggplot(freqDF, aes(x=inTAD, y=avgCount, fill=group)) +
    geom_bar(stat="identity", position="dodge", colour="black") + 
    facet_grid(.~distTadBin) + 
    geom_errorbar(aes(ymax = avgCount + sdCount , ymin=avgCount - sdCount), position=position_dodge(width=0.9), width=.25) +
    geom_text(aes(label=round(avgCount), y= avgCount), position=position_dodge(width=1), vjust=1.25, size=5) +
    scale_fill_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + labs(y="Gene pairs [n]", x="")
ggsave(p, file=paste0(outPrefix, ".gene_pairs.by_group_inTAD_distTadBin.barplot.pdf"), w=7, h=7)

p <- ggplot(freqDF, aes(x=inTAD, y=avgPercent, fill=group)) +
    geom_bar(stat="identity", position="dodge", colour="black") + 
    facet_grid(.~distTadBin) + 
    geom_errorbar(aes(ymax = avgPercent + sdPercent , ymin=avgPercent - sdPercent), position=position_dodge(width=0.9), width=.25) +
    geom_text(aes(label=round(avgPercent), y= avgPercent), position=position_dodge(width=1), vjust=1.25, size=5) +
    scale_fill_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + labs(y="Gene pairs [n]", x="")
ggsave(p, file=paste0(outPrefix, ".percent_gene_pairs.by_group_in_TAD_distTadBin.barplot.pdf"), w=7, h=7)


#-----------------------------------------------------------------------
# percent in TAD by group and TADsource for [10-1000kb]
#-----------------------------------------------------------------------

# build data frame with all tadSources
subDF <- allDF[sSamp & sExp & sSpec & allDF$distTadBin == "10-1000kb",]


# count gene pairs by combinations
dc <- ddply(subDF, .(group, replicate, tadSource, inTAD), summarize, count=length(inTAD), .drop=FALSE)
# flag get existing replicate/group combination
existingComb <- dc$replicate <= 1 | dc$group == "sampled"
# sum up counts for sub-combination to get percentage values
summedCounts <- rep(ddply(dc, .(group, replicate, tadSource), summarize, s=sum(count))$s, each=length(levels(dc$inTAD)))
# add percentage values and counts by marking non-existing combinations with NA
freqRepDF <- cbind(dc, cbind(summedCounts, existingCount=ifelse(existingComb, dc$count, NA), percent=dc$count * 100 / summedCounts))

# combine replicates by taking average and sd of counts and percentages
freqDF <- ddply(freqRepDF, .(group, tadSource,inTAD), summarize, avgCount=mean(existingCount, na.rm=TRUE), sdCount=sd(existingCount, na.rm=TRUE), avgPercent=mean(percent, na.rm=TRUE), sdPercent=sd(percent, na.rm=TRUE))

# calculate p-values
pvalDF <- ddply(subDF, .(tadSource), summarize, p=fisher.test(group, inTAD)$p.value)

# combine p-values with plotting coordinates (max of values per group)
pvalDF <- merge(pvalDF, ddply(freqDF[freqDF$inTAD=="same TAD",], .(tadSource), summarize, avgCount=max(avgCount)+50, avgPercent=max(avgPercent)+5))
pvalDF[,"group"] <- "paralog"

p <- ggplot(freqDF[freqDF$inTAD=="same TAD",], aes(x=tadSource, y=avgCount, fill=group)) +
    geom_errorbar(aes(ymax = avgCount + sdCount , ymin=avgCount - sdCount), position=position_dodge(width=0.9), width=.25) +
    geom_bar(stat="identity", position="dodge", colour="black") + 
    geom_text(aes(label=round(avgCount), y= avgCount), position=position_dodge(width=1), vjust=1.25, size=5) +
    geom_text(aes(label=paste0("p=", signif(p,2))), data=pvalDF, size=5) + 
    scale_fill_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + labs(y="Gene pairs in same TAD", x="")
ggsave(p, file=paste0(outPrefix, ".10-1000kb.gene_pairs.by_group_in_TAD_tadSource.barplot.pdf"), w=14, h=7)

p <- ggplot(freqDF[freqDF$inTAD=="same TAD",], aes(x=tadSource, y=avgPercent, fill=group)) +
    geom_errorbar(aes(ymax = avgPercent + sdPercent , ymin=avgPercent - sdPercent), position=position_dodge(width=0.9), width=.25) +
    geom_bar(stat="identity", position="dodge", colour="black") + 
    geom_text(aes(label=signif(avgPercent,3), y= avgPercent), position=position_dodge(width=1), vjust=1.25, size=5) +
    geom_text(aes(label=paste0("p=", signif(p,2))), data=pvalDF, size=5) + 
    scale_fill_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + labs(y="Gene pairs in same TAD [%]", x="")
ggsave(p, file=paste0(outPrefix, ".10-1000kb.percent_gene_pairs.by_group_in_TAD_TADsource.barplot.pdf"), w=14, h=7)


#-----------------------------------------------------------------------
# percent in TAD by group and TADsource for close pairs only
#-----------------------------------------------------------------------

# build data frame with all tadSources
subDF <- allDF[sSamp & sExp & sSpec & allDF$distGroup == "close",]


# count gene pairs by combinations
dc <- ddply(subDF, .(group, replicate, tadSource, inTAD), summarize, count=length(inTAD), .drop=FALSE)
# flag get existing replicate/group combination
existingComb <- dc$replicate <= 1 | dc$group == "sampled"
# sum up counts for sub-combination to get percentage values
summedCounts <- rep(ddply(dc, .(group, replicate, tadSource), summarize, s=sum(count))$s, each=length(levels(dc$inTAD)))
# add percentage values and counts by marking non-existing combinations with NA
freqRepDF <- cbind(dc, cbind(summedCounts, existingCount=ifelse(existingComb, dc$count, NA), percent=dc$count * 100 / summedCounts))

# combine replicates by taking average and sd of counts and percentages
freqDF <- ddply(freqRepDF, .(group, tadSource,inTAD), summarize, avgCount=mean(existingCount, na.rm=TRUE), sdCount=sd(existingCount, na.rm=TRUE), avgPercent=mean(percent, na.rm=TRUE), sdPercent=sd(percent, na.rm=TRUE))

# calculate p-values
pvalDF <- ddply(subDF, .(tadSource), summarize, p=fisher.test(group, inTAD)$p.value)

# combine p-values with plotting coordinates (max of values per group)
pvalDF <- merge(pvalDF, ddply(freqDF[freqDF$inTAD=="same TAD",], .(tadSource), summarize, avgCount=max(avgCount)+50, avgPercent=max(avgPercent)+5))
pvalDF[,"group"] <- "paralog"

p <- ggplot(freqDF[freqDF$inTAD=="same TAD",], aes(x=tadSource, y=avgCount, fill=group)) +
    geom_errorbar(aes(ymax = avgCount + sdCount , ymin=avgCount - sdCount), position=position_dodge(width=0.9), width=.25) +
    geom_bar(stat="identity", position="dodge", colour="black") + 
    geom_text(aes(label=round(avgCount), y= avgCount), position=position_dodge(width=1), vjust=1.25, size=5) +
    geom_text(aes(label=paste0("p=", signif(p,2))), data=pvalDF, size=5) + 
    scale_fill_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + labs(y="Gene pairs in same TAD", x="")
ggsave(p, file=paste0(outPrefix, ".close.gene_pairs.by_group_in_TAD_tadSource.barplot.pdf"), w=14, h=7)

p <- ggplot(freqDF[freqDF$inTAD=="same TAD",], aes(x=tadSource, y=avgPercent, fill=group)) +
    geom_errorbar(aes(ymax = avgPercent + sdPercent , ymin=avgPercent - sdPercent), position=position_dodge(width=0.9), width=.25) +
    geom_bar(stat="identity", position="dodge", colour="black") + 
    geom_text(aes(label=signif(avgPercent,3), y= avgPercent), position=position_dodge(width=1), vjust=1.25, size=5) +
    geom_text(aes(label=paste0("p=", signif(p,2))), data=pvalDF, size=5) + 
    scale_fill_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + labs(y="Gene pairs in same TAD [%]", x="")
ggsave(p, file=paste0(outPrefix, ".close.percent_gene_pairs.by_group_in_TAD_TADsource.barplot.pdf"), w=14, h=7)


#------------------------------------------------------------------------
# HiC by group and inTADclose
#------------------------------------------------------------------------

# reduce to sampled by distance and only one source for TAD/exp/species
subDF <- allDF[sSamp & sTAD & sExp & sSpec,]


# iterate over all contact measurements (in HiCcolumns):
#~ for (HiCcol in c("HiC", "captureC_raw", "HiCNoZero", "captureC_rawNoZero")){
for (HiCcol in HiCcolumns){
    
    message(paste("INFO: Plotting:", HiCcol))
    HiClab <- HiClabels[HiCcol]
    
    nDF <- ddply(subDF, .(inTADclose, group), summarize, n=sum(!is.na(get(HiCcol))))
    
    # calculate p-values
    pvalDF <- ddply(subDF, .(inTADclose), summarize, p=wilcox.test(as.formula(paste(HiCcol, "~ group")))$p.value)

    p = ggplot(subDF, aes_string(x="group", y=HiCcol)) + 
        geom_boxplot(aes(colour = group), lwd=1.5) + scale_y_log10() +
        facet_grid(.~inTADclose) + 
        scale_color_manual(values=COL, guide_legend(title = "")) +
        theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="none") + 
        guides(fill=guide_legend(title="")) +
        labs(x="",y=HiClab) + 
        geom_text(aes(label=paste0("n=",n),  y=0.25), data=nDF) +
        geom_text(aes(label=paste0("p=", signif(p,2)), x=1.5, y=250), data=pvalDF, size=5)   
    ggsave(p, file=paste0(outPrefix, ".", HiCcol, ".by_inTADclose.boxplot.pdf"), w=7, h=7)

}

#-----------------------------------------------------------------------
# number of gene pairs by group and subTAD combination
#-----------------------------------------------------------------------

# restrict to the close pairs with distance within 10-1000kb
subDF <- allDF[sSamp & sTAD & sExp & sSpec & allDF$distTadBin == "10-1000kb",]

dc <- ddply(subDF, .(group, subTAD, replicate), summarize, count=length(subTAD))
summedCounts <- rep(ddply(dc, .(group, replicate), summarize, s=sum(count))$s, each=length(levels(dc$subTAD)))
freqRepDF <- cbind(dc, cbind(summedCounts, percent=dc$count * 100 / summedCounts))
freqDF <- ddply(freqRepDF, .(group, subTAD), summarize, avgCount=mean(count, na.rm=TRUE), sdCount=sd(count, na.rm=TRUE), avgPercent=mean(percent), sdPercent=sd(percent))

# calculate p-values
pvals <- sapply(levels(subDF$subTAD), function(sT) fisher.test(table(subDF$subTAD == sT, subDF$group))$p.value) 
pvalDF <- data.frame(
    subTAD=levels(subDF$subTAD),
    p = pvals,
    ypos = ddply(freqDF, .(subTAD), summarize, y=max(avgCount))$y + 50,
    group=NA,
    dupAgeGroup=NA
)

p <- ggplot(freqDF, aes(x=group, y=avgCount, fill=group)) +
    geom_errorbar(aes(ymax = avgCount + sdCount , ymin=avgCount - sdCount), position=position_dodge(width=0.9), width=.25) +
    geom_bar(stat="identity", position="dodge", colour="black") + 
    facet_grid(.~subTAD, labeller=as_labeller(subTADlabels)) +
    scale_fill_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_blank(), legend.position="bottom") + labs(y="Gene pairs", x="") +
    geom_text(aes(label=avgCount, y= avgCount), position=position_dodge(width=1), vjust=1.25, size=5) +
    geom_text(aes(label=paste0("p=",signif(p,2)), y=ypos, x=1.5), size=5, data=pvalDF)

g <- addPictureLabels(p, subTADfigPaths)
pdf(paste0(outPrefix, ".10_1000kb.subTAD.byGroup.barplot.pdf"))
    grid.draw(g)
dev.off()

#-----------------------------------------------------------------------
# percent of gene pairs by group and subTAD combination
#-----------------------------------------------------------------------
p <- ggplot(freqDF, aes(x=group, y=avgPercent, fill=group)) +
    geom_errorbar(aes(ymax = avgPercent + sdPercent , ymin=avgPercent - sdPercent), position=position_dodge(width=0.9), width=.25) +
    geom_bar(stat="identity", position="dodge", colour="black") +     
    facet_grid(.~subTAD, labeller=as_labeller(subTADlabels)) +
    scale_fill_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_blank(), legend.position="bottom") + labs(y="Gene pairs [%]", x="") +
    geom_text(aes(label=signif(avgPercent,3)), position=position_dodge(width=1), vjust=1.25, size=5) +
    geom_text(aes(label=paste0("p=",signif(p,2)), y=1.1*max(freqDF$avgPercent), x=1.5), size=5, data=pvalDF)

g <- addPictureLabels(p, subTADfigPaths)
pdf(paste0(outPrefix, ".10_1000kb.percent_gene_pairs.subTAD.byGroup.barplot.pdf"))
    grid.draw(g)
dev.off()

#-----------------------------------------------------------------------
# number of gene pairs by dupAgeGroup, group and subTAD combination
#-----------------------------------------------------------------------
# make dupAge for sampled pairs as NA
subDF[subDF$group == "sampled", "dupAgeGroup"] <- NA

dcAgeGroup <- ddply(subDF, .(group, dupAgeGroup, replicate), function(d) {data.frame(
    table(d$subTAD)
    )})
names(dcAgeGroup)[c(4,5)] <- c("subTAD", "count")
freqAgeGroupDF <- ddply(dcAgeGroup, .(group, dupAgeGroup, subTAD), summarize, avgCount=mean(count, na.rm=TRUE), sdCount=sd(count, na.rm=TRUE))
# convert to factor that do not exclude NAs for plotting
freqAgeGroupDF$dupAgeGroup <- factor(freqAgeGroupDF$dupAgeGroup, exclude = NULL)

p <- ggplot(freqAgeGroupDF, aes(x=group, y=avgCount, fill=dupAgeGroup)) +
    geom_bar(stat="identity", colour="black") + 
    geom_errorbar(aes(ymax = avgCount + sdCount , ymin=avgCount - sdCount), width=.25) +
    facet_grid(.~subTAD, labeller=as_labeller(subTADlabels)) +
    scale_fill_manual(values=c(COL_AGE, COL[2]), guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + labs(y="Gene pairs", x="") +
    geom_text(aes(label=avgCount, y= avgCount), vjust=1.25, size=5) +
    geom_text(aes(label=paste0("p=",signif(p,2)), y=ypos, x=1.5), size=5, data=pvalDF)

g <- addPictureLabels(p, subTADfigPaths)

pdf(paste0(outPrefix, ".10_1000kb.subTAD.by_dupAgeGroup_and_group.barplot.pdf"))
    grid.draw(g)
dev.off()

#-----------------------------------------------------------------------
# number of gene pairs with shared enhancer by group and subTAD combination
#-----------------------------------------------------------------------

# get subset for sampled by distance and enhancer and only the distance 10-1000kb
subDF <- allDF[sSampEh & sTAD & sExp & sSpec & allDF$distTadBin == "10-1000kb",]

dc <- ddply(subDF, .(group, subTAD, replicate), summarize, n=length(eh), count=sum(eh), percent=sum(eh)/length(eh)*100)
freqDF <- ddply(dc, .(group, subTAD), summarize, 
    avgCount=mean(count, na.rm=TRUE), 
    sdCount=sd(count, na.rm=TRUE), 
    avgPercent=mean(percent), 
    sdPercent=sd(percent), 
    avgN=mean(n), 
    sdN=sd(n)
)

# calculate p-values
pvals <- sapply(levels(subDF$subTAD), function(sT) fisher.test(table(subDF[subDF$subTAD == sT,"eh"], subDF[subDF$subTAD == sT,"group"]))$p.value) 

pvalDF <- data.frame(
    subTAD=names(pvals),
    p = pvals,
    ypos = ddply(freqDF, .(subTAD), summarize, y=max(avgCount)+10)$y,
    yposPercent = ddply(freqDF, .(subTAD), summarize, y=max(avgPercent)+2)$y,
    group=NA
)

p <- ggplot(freqDF, aes(x=group, y=avgCount, fill=group)) +
    geom_bar(stat="identity", colour="black") + 
    geom_errorbar(aes(ymax = avgCount + sdCount , ymin=avgCount - sdCount), width=.25) +
    facet_grid(.~subTAD) +
    scale_fill_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_blank(), legend.position="bottom") + labs(y="Gene pairs with shared enhancer", x="") + theme(strip.text.x = element_text(size=4, angle=90)) +
    geom_text(aes(label=avgCount, y= avgCount), vjust=1.25, size=5) +
    geom_text(aes(label=paste0("p=",signif(p,2)), y=ypos, x=1.5), size=5, data=pvalDF)

g <- addPictureLabels(p, subTADfigPaths)
pdf(paste0(outPrefix, ".10_1000kb.eh.by_subTAD_and_group.barplot.pdf"))
    grid.draw(g)
dev.off()

# same with percent values
p <- ggplot(freqDF, aes(x=group, y=avgPercent, fill=group)) +
    geom_bar(stat="identity", colour="black") + 
    geom_errorbar(aes(ymax = avgPercent + sdPercent , ymin=avgPercent - sdPercent), width=.25) +
    facet_grid(.~subTAD) +
    scale_fill_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_blank(), legend.position="bottom") + labs(y="Percent of pairs with shared enhancer  [%]", x="") + theme(strip.text.x = element_text(size=4, angle=90)) +
    geom_text(aes(label=signif(avgPercent,3), y= avgPercent), vjust=1.25, size=5) +
    geom_text(aes(label=paste0("p=",signif(p,2)), y=yposPercent, x=1.5), size=5, data=pvalDF)

g <- addPictureLabels(p, subTADfigPaths)
pdf(paste0(outPrefix, ".10_1000kb.ehPercent.by_subTAD_and_group.barplot.pdf"))
    grid.draw(g)
dev.off()

#-----------------------------------------------------------------------
# subTAD by group
#-----------------------------------------------------------------------

# get subset on all distances
subDF <- allDF[sSamp & sTAD & sExp & sSpec,]


dc <- ddply(subDF, .(group, subTAD, replicate), summarize, count=length(subTAD))
summedCounts <- rep(ddply(dc, .(group, replicate), summarize, s=sum(count))$s, each=length(levels(dc$subTAD)))
freqRepDF <- cbind(dc, cbind(summedCounts, percent=dc$count * 100 / summedCounts))
freqDF <- ddply(freqRepDF, .(group, subTAD), summarize, avgPercent=mean(percent), sd=sd(percent))

p <- ggplot(freqDF, aes(x=subTAD, y=avgPercent, fill=group)) +
    geom_bar(stat="identity", position="dodge", colour="black") + 
    geom_errorbar(aes(ymax = avgPercent + sd , ymin=avgPercent - sd), position=position_dodge(width=0.9), width=.25) +
    geom_text(aes(label=signif(avgPercent, 2), y= avgPercent), position=position_dodge(width=1), vjust=1.25, size=5) +
    scale_fill_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + labs(y="Gene pairs [%]", x="")
ggsave(p, file=paste0(outPrefix, ".subTAD.byGroup.barplot.pdf"), w=3.5, h=7)


#-----------------------------------------------------------------------
# HiC by group and subTAD combination
#-----------------------------------------------------------------------
subDF <- allDF[sSamp & sTAD & sExp & sSpec & allDF$distTadBin == "10-1000kb",]

# iterate over all contact measurements (in HiCcolumns):
#~ for (HiCcol in HiCcolumns){
for (HiCcol in c("HiC", "captureC_raw", "HiCNoZero", "captureC_rawNoZero")){
    
    message(paste("INFO: Plotting:", HiCcol))
    HiClab <- HiClabels[HiCcol]

    nDF <- ddply(subDF, .(subTAD, group), summarize, n=sum(!is.na(get(HiCcol))))
    
    # calculate p-values
    pvalDF <- ddply(subDF, .(subTAD), summarize, p=wilcox.test(as.formula(paste(HiCcol, "~ group")))$p.value)


    p = ggplot(subDF, aes_string(x="group", y=HiCcol)) + 
        geom_boxplot(aes(colour = group), lwd=1.5) + scale_y_log10() +
        facet_grid(.~subTAD, labeller=as_labeller(subTADlabels)) + 
        scale_color_manual(values=COL, guide_legend(title = "")) +
        theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="none") + 
        guides(fill=guide_legend(title="")) +
        labs(x="", y=HiClab) + 
        geom_text(aes(label=paste0("n=",n),  y=0.25), data=nDF) +
        geom_text(aes(label=paste0("p=", signif(p,2)), x=1.5, y=250), data=pvalDF, size=5) 
          
    g <- addPictureLabels(p, subTADfigPaths)
    pdf(paste0(outPrefix, ".10_1000kb.", HiCcol, ".by_subTAD.boxplot.pdf"))
        grid.draw(g)
    dev.off()
    
}


#-----------------------------------------------------------------------
# SUBTAD WITH 3 GROUPS: number of gene pairs by group and sub3TAD combination
#-----------------------------------------------------------------------
subDF <- allDF[sSamp & sTAD & sExp & sSpec & allDF$distTadBin == "10-1000kb",]

dc <- ddply(subDF, .(group, sub3TAD, replicate), summarize, count=length(sub3TAD))
summedCounts <- rep(ddply(dc, .(group, replicate), summarize, s=sum(count))$s, each=length(levels(dc$sub3TAD)))
freqRepDF <- cbind(dc, cbind(summedCounts, percent=dc$count * 100 / summedCounts))
freqDF <- ddply(freqRepDF, .(group, sub3TAD), summarize, avgCount=mean(count, na.rm=TRUE), sdCount=sd(count, na.rm=TRUE), avgPercent=mean(percent), sd=sd(percent))

# calculate p-values
pvals <- sapply(levels(subDF$sub3TAD), function(sT) fisher.test(table(subDF$sub3TAD == sT, subDF$group))$p.value) 
pvalDF <- data.frame(
    sub3TAD=levels(subDF$sub3TAD),
    p = pvals,
    ypos = ddply(freqDF, .(sub3TAD), summarize, y=max(avgCount))$y + 50,
    group=NA,
    dupAgeGroup=NA
)


p <- ggplot(freqDF, aes(x=group, y=avgCount, fill=group)) +
    geom_errorbar(aes(ymax = avgCount + sdCount , ymin=avgCount - sdCount), position=
    position_dodge(width=0.9), width=.25) +
    geom_bar(stat="identity", position="dodge", colour="black") + 
    facet_grid(.~sub3TAD) +
    scale_fill_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_blank(), legend.position="bottom") + labs(y="Gene pairs", x="") + theme(strip.text.x = element_text(size=5, angle=90)) +
    geom_text(aes(label=avgCount, y= avgCount), position=position_dodge(width=1), vjust=1.25, size=5) +
    geom_text(aes(label=paste0("p=",signif(p,2)), y=ypos, x=1.5), size=5, data=pvalDF)


g <- addPictureLabels(p, sub3TADfigPaths)
pdf(paste0(outPrefix, ".10_1000kb.sub3TAD.byGroup.barplot.pdf"))
    grid.draw(g)
dev.off()

#-----------------------------------------------------------------------
# number of gene pairs by dupAgeGroup, group and sub3TAD combination
#-----------------------------------------------------------------------
# make dupAge for sampled pairs as NA
subDF[subDF$group == "sampled", "dupAgeGroup"] <- NA

dcAgeGroup <- ddply(subDF, .(group, dupAgeGroup, replicate), function(d) {data.frame(
    table(d$sub3TAD)
    )})
names(dcAgeGroup)[c(4,5)] <- c("sub3TAD", "count")
freqAgeGroupDF <- ddply(dcAgeGroup, .(group, dupAgeGroup, sub3TAD), summarize, avgCount=mean(count, na.rm=TRUE), sdCount=sd(count, na.rm=TRUE))
# convert to factor that do not exclude NAs for plotting
freqAgeGroupDF$dupAgeGroup <- factor(freqAgeGroupDF$dupAgeGroup, exclude = NULL)

p <- ggplot(freqAgeGroupDF, aes(x=group, y=avgCount, fill=dupAgeGroup)) +
    geom_errorbar(aes(ymax = avgCount + sdCount , ymin=avgCount - sdCount), width=.25) +
    geom_bar(stat="identity", colour="black") + 
    facet_grid(.~sub3TAD) +
    scale_fill_manual(values=c(COL_AGE, COL[2]), guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + labs(y="Gene pairs", x="") + theme(strip.text.x = element_text(size=5, angle=90)) +
    geom_text(aes(label=avgCount, y= avgCount), vjust=1.25, size=5) +
    geom_text(aes(label=paste0("p=",signif(p,2)), y=ypos, x=1.5), size=5, data=pvalDF)

g <- addPictureLabels(p, sub3TADfigPaths)

pdf(paste0(outPrefix, ".10_1000kb.sub3TAD.by_dupAgeGroup_and_group.barplot.pdf"))
    grid.draw(g)
dev.off()

#-----------------------------------------------------------------------
# number of gene pairs with shared enhancer by group and sub3TAD combination
#-----------------------------------------------------------------------

# get subset for sampled by distance and enhancer and only the distance 10-1000kb
subDF <- allDF[sSampEh & sTAD & sExp & sSpec & allDF$distTadBin == "10-1000kb",]

dc <- ddply(subDF, .(group, sub3TAD, replicate), summarize, n=length(eh), count=sum(eh), percent=sum(eh)/length(eh)*100)
freqDF <- ddply(dc, .(group, sub3TAD), summarize, 
    avgCount=mean(count, na.rm=TRUE), 
    sdCount=sd(count, na.rm=TRUE), 
    avgPercent=mean(percent), 
    sdPercent=sd(percent), 
    avgN=mean(n), 
    sdN=sd(n)
)

# calculate p-values
pvals <- sapply(levels(subDF$sub3TAD), function(sT) fisher.test(table(subDF[subDF$sub3TAD == sT,"eh"], subDF[subDF$sub3TAD == sT,"group"]))$p.value) 

pvalDF <- data.frame(
    sub3TAD=names(pvals),
    p = pvals,
    ypos = ddply(freqDF, .(sub3TAD), summarize, y=max(avgCount)+10)$y,
    yposPercent = ddply(freqDF, .(sub3TAD), summarize, y=max(avgPercent)+2)$y,
    group=NA
)

p <- ggplot(freqDF, aes(x=group, y=avgCount, fill=group)) +
    geom_bar(stat="identity", colour="black") + 
    geom_errorbar(aes(ymax = avgCount + sdCount , ymin=avgCount - sdCount), width=.25) +
    facet_grid(.~sub3TAD) +
    scale_fill_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_blank(), legend.position="bottom") + labs(y="Gene pairs with shared enhancer", x="") + theme(strip.text.x = element_text(size=5, angle=90)) +
    geom_text(aes(label=avgCount, y= avgCount), vjust=1.25, size=5) +
    geom_text(aes(label=paste0("p=",signif(p,2)), y=ypos, x=1.5), size=5, data=pvalDF)

g <- addPictureLabels(p, sub3TADfigPaths)
pdf(paste0(outPrefix, ".10_1000kb.eh.by_sub3TAD_and_group.barplot.pdf"))
    grid.draw(g)
dev.off()

#-----------------------------------------------------------------------
# same with percent values
#-----------------------------------------------------------------------
p <- ggplot(freqDF, aes(x=group, y=avgPercent, fill=group)) +
    geom_bar(stat="identity", colour="black") + 
    geom_errorbar(aes(ymax = avgPercent + sdPercent , ymin=avgPercent - sdPercent), width=.25) +
    facet_grid(.~sub3TAD) +
    scale_fill_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_blank(), legend.position="bottom") + labs(y="Percent of pairs with shared enhancer  [%]", x="") + theme(strip.text.x = element_text(size=5, angle=90)) +
    geom_text(aes(label=signif(avgPercent,3), y= avgPercent), vjust=1.25, size=5) +
    geom_text(aes(label=paste0("p=",signif(p,2)), y=yposPercent, x=1.5), size=5, data=pvalDF)

g <- addPictureLabels(p, sub3TADfigPaths)
pdf(paste0(outPrefix, ".10_1000kb.ehPercent.by_sub3TAD_and_group.barplot.pdf"))
    grid.draw(g)
dev.off()

#-----------------------------------------------------------------------
# sub3TAD by group
#-----------------------------------------------------------------------
subDF <- allDF[sSamp & sTAD & sExp & sSpec,]

dc <- ddply(subDF, .(group, sub3TAD, replicate), summarize, count=length(sub3TAD))
summedCounts <- rep(ddply(dc, .(group, replicate), summarize, s=sum(count))$s, each=length(levels(dc$sub3TAD)))
freqRepDF <- cbind(dc, cbind(summedCounts, percent=dc$count * 100 / summedCounts))
freqDF <- ddply(freqRepDF, .(group, sub3TAD), summarize, avgPercent=mean(percent), sd=sd(percent))

p <- ggplot(freqDF, aes(x=sub3TAD, y=avgPercent, fill=group)) +
    geom_bar(stat="identity", position="dodge", colour="black") + 
    geom_errorbar(aes(ymax = avgPercent + sd , ymin=avgPercent - sd), position=position_dodge(width=0.9), width=.25) +
    geom_text(aes(label=signif(avgPercent, 2), y= avgPercent), position=position_dodge(width=1), vjust=1.25, size=5) +
    scale_fill_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + labs(y="Gene pairs [%]", x="")
ggsave(p, file=paste0(outPrefix, ".sub3TAD.byGroup.barplot.pdf"), w=3.5, h=7)


#-----------------------------------------------------------------------
# HiC by group and sub3TAD combination
#-----------------------------------------------------------------------
subDF <- allDF[sSamp & sTAD & sExp & sSpec & allDF$distTadBin == "10-1000kb",]

# iterate over all contact measurements (in HiCcolumns):
#~ for (HiCcol in HiCcolumns){
for (HiCcol in c("HiC", "captureC_raw", "HiCNoZero", "captureC_rawNoZero")){
    
    message(paste("INFO: Plotting:", HiCcol))
    HiClab <- HiClabels[HiCcol]
    
    nDF <- ddply(subDF, .(sub3TAD, group), summarize, n=sum(!is.na(get(HiCcol))))
    
    # calculate p-values
    pvalDF <- ddply(subDF, .(sub3TAD), summarize, p=wilcox.test(as.formula(paste(HiCcol, "~ group")))$p.value)
    p = ggplot(subDF, aes_string(x="group", y=HiCcol)) + 
        geom_boxplot(aes(colour = group), lwd=1.5) + scale_y_log10() +
        facet_grid(.~sub3TAD) + 
        scale_color_manual(values=COL, guide_legend(title = "")) +
        theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="none") + 
        guides(fill=guide_legend(title="")) +
        labs(x="", y=HiClab) + 
        geom_text(aes(label=paste0("n=",n),  y=0.25), data=nDF) +
        geom_text(aes(label=paste0("p=", signif(p,2)), x=1.5, y=250), data=pvalDF, size=5)   
    ggsave(p, file=paste0(outPrefix, ".10_1000kb.", HiCcol, ".by_sub3TAD.boxplot.pdf"), w=7, h=7)
}

#-----------------------------------------------------------------------
# Density of average expression of pairs in IMR90
#-----------------------------------------------------------------------
subDF <- allDF[sSamp & sTAD & sExp & sSpec,]

p <- ggplot(subDF, aes(avgExp+1, ..density.., fill = group, color=group)) +
   geom_density(alpha=.5) + scale_x_log10() + 
   scale_color_manual(values=COL, guide_legend(title = "")) + scale_fill_manual(values=COL, guide_legend(title = "")) +
   theme_bw() + theme(text = element_text(size=20), legend.justification=c(1,1), legend.position=c(1,1))
ggsave(p, file=paste0(outPrefix, ".avgExp_IMR90.by_group.density.pdf"), w=3.5, h=3.5)    

#-----------------------------------------------------------------------
# expression of pairs in IMR90 by HiC for close pairs
#-----------------------------------------------------------------------
subDF <- allDF[sSamp & sTAD & sExp & sSpec & allDF$distGroup == "close",]

p <- dotplotWithDensityLogXY(revDF(subDF[subDF$dist > 1,]), "avgExp+1", "HiC", "group", COL, fit=TRUE, ALPHA=.2)
pdf(paste0(outPrefix, ".close.HiC.vs_avgExp.dotplot.pdf"))
    grid::grid.draw(p)
dev.off()

#-----------------------------------------------------------------------
# expression of pairs in IMR90 by HiC
#-----------------------------------------------------------------------
subDF <- allDF[sSamp & sTAD & sExp & sSpec,]

p <- dotplotWithDensityLogXY(revDF(subDF[subDF$dist > 1,]), "avgExp+1", "HiC", "group", COL, fit=TRUE, ALPHA=.2)
pdf(paste0(outPrefix, ".HiC.vs_avgExp.dotplot.pdf"))
    grid::grid.draw(p)
dev.off()

#-----------------------------------------------------------------------
# expression correlation by TAD and TAD distance bins
#-----------------------------------------------------------------------
subDF <- allDF[sSamp & sTAD & sExp & sSpec,]

nDF <- data.frame(table(subDF[, c("group", "distTadBin", "inTAD")], useNA="ifany"))
p = ggplot(subDF, aes(x=inTAD, y=expCor^2)) + 
    geom_boxplot(aes(colour = group), lwd=1.5) + 
    facet_grid(.~distTadBin) + 
    scale_color_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + 
    guides(fill=guide_legend(title="")) +
    labs(y="Expression Correlation [R^2]", x="Linear distance bin [kb]") + 
    geom_text(aes(label=paste0("n=",Freq), vjust=(as.numeric(group)-1)*2-1,  y=0, angle=90), data=nDF)
ggsave(p, file=paste0(outPrefix, ".expCor.by_TADdist_byTAD.boxplot.pdf"), w=7, h=7)

#-----------------------------------------------------------------------
# expression correlation by subTAD and TAD distance bins
#-----------------------------------------------------------------------
nDF <- data.frame(table(subDF[, c("group", "distTadBin", "subTAD")], useNA="ifany"))
p = ggplot(subDF, aes(x=subTAD, y=expCor^2)) + 
    geom_boxplot(aes(colour = group), lwd=1.5) + 
    facet_grid(.~distTadBin) + 
    scale_color_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + 
    guides(fill=guide_legend(title="")) +
    labs(y="Expression Correlation [R^2]", x="Linear distance bin [kb]") + 
    geom_text(aes(label=paste0("n=",Freq), vjust=(as.numeric(group)-1)*2-1,  y=0, angle=90), data=nDF)    
ggsave(p, file=paste0(outPrefix, ".expCor.by_TADdist_subTAD.boxplot.pdf"), w=7, h=7)

#-----------------------------------------------------------------------
# expression correlation by subTAD and group in 10-1000kb
#-----------------------------------------------------------------------
subDF <- allDF[sSamp & sTAD & sExp & sSpec & allDF$distTadBin == "10-1000kb",]

nDF <- data.frame(table(subDF[, c("group", "subTAD")], useNA="ifany"))
p = ggplot(subDF, aes(x=group, y=expCor^2)) + 
    geom_boxplot(aes(colour = group), lwd=1.5) + 
    facet_grid(.~subTAD) + 
    scale_color_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="none") + 
    guides(fill=guide_legend(title="")) +
    labs(y="Expression Correlation [R^2]", x="") + 
    geom_text(aes(label=paste0("n=",Freq),  y=1.05), data=nDF)    
ggsave(p, file=paste0(outPrefix, ".10_1000kb.expCor.by_subTAD.boxplot.pdf"), w=7, h=7)

#=======================================================================
# Protein protein interaction from HIPPIE
#=======================================================================

#-----------------------------------------------------------------------
# HIPPIE by TAD and distTadBin
#-----------------------------------------------------------------------
subDF <- allDF[sSamp & sTAD & sExp & sSpec,]

nDF <- ddply(subDF, .(group, inTAD, distTadBin), summarize, Freq=sum(!is.na(HIPPIE)))
p = ggplot(subDF, aes(x=group, y=HIPPIE)) + 
    geom_boxplot(aes(colour = group), lwd=1.5) +
    facet_grid(.~distTadBin*inTAD) + 
    scale_color_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="none") + 
    guides(fill=guide_legend(title="")) +
    labs(x="") + 
    geom_text(aes(label=paste0("n=",Freq),  y=1), data=nDF) 
ggsave(p, file=paste0(outPrefix, ".HIPPIEscore.by_TAD_and_distTadBin.boxplot.pdf"), w=7, h=7)

#-----------------------------------------------------------------------
# HIPPIE by TAD and groups [10kb - 1000kb]
#-----------------------------------------------------------------------
subDF <- allDF[sSamp & sTAD & sExp & sSpec & allDF$distTadBin == "10-1000kb",]

pvalDF <- ddply(subDF, .(inTAD), summarize, p=wilcox.test(HIPPIE ~ group)$p.value)

nDF <- ddply(subDF, .(group, inTAD), summarize, Freq=sum(!is.na(HIPPIE)))
p = ggplot(subDF, aes(x=group, y=HIPPIE)) + 
    geom_boxplot(aes(colour = group), lwd=1.5) +
    facet_grid(.~inTAD) + 
    scale_color_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="none") + 
    guides(fill=guide_legend(title="")) +
    labs(x="") + 
    geom_text(aes(label=paste0("n=",Freq),  y=-0.05), data=nDF) + 
    geom_text(aes(label=paste0("p=", signif(p,3)), y=1, x=1.5), data=pvalDF, size=5)
ggsave(p, file=paste0(outPrefix, ".10_1000kb.HIPPIEscore.by_group_and_TAD.boxplot.pdf"), w=3.5, h=7)

#-----------------------------------------------------------------------
# HIPPIE by distTadBin and group
#-----------------------------------------------------------------------
subDF <- allDF[sSamp & sTAD & sExp & sSpec,]
nDF <- ddply(subDF, .(group, distTadBin), summarize, Freq=sum(!is.na(HIPPIE)))

p <- ggplot(subDF, aes(x=group, y=HIPPIE)) + 
    geom_boxplot(aes(colour = group), lwd=1.5) +
    facet_grid(.~distTadBin) + 
    scale_color_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="none") + 
    guides(fill=guide_legend(title="")) +
    labs(x="") + 
    geom_text(aes(label=paste0("n=",Freq),  y=1), data=nDF)    
ggsave(p, file=paste0(outPrefix, ".HIPPIEscore.by_distTadBin_and_group.boxplot.pdf"), w=3.5, h=7)

#-----------------------------------------------------------------------
# PPI by TAD and group [10kb - 1000kb]
#-----------------------------------------------------------------------
subDF <- subset(allDF, sSamp & sTAD & sExp & sSpec & distTadBin == "10-1000kb")
#~ subDF <- subset(allDF, sSamp & sTAD & sExp & sSpec & distGroup == "close")

dc <- ddply(subDF, .(group, inTAD, replicate, PPI), summarize, count=length(PPI), .drop=FALSE)
summedCounts <- rep(ddply(dc, .(group, inTAD, replicate), summarize, s=sum(count))$s, each=length(levels(dc$PPI))+1)
freqRepDF <- cbind(dc, cbind(summedCounts, percent=dc$count * 100 / summedCounts))
freqDF <- ddply(freqRepDF, .(group, inTAD, PPI), summarize, avg=sum(count), avgPercent=mean(percent,na.rm=TRUE), sd=sd(percent,na.rm=TRUE))

p <- ggplot(freqDF, aes(x=PPI, y=avgPercent, fill=group)) +
    geom_errorbar(aes(ymax = avgPercent + sd , ymin=avgPercent - sd), position=position_dodge(width=0.9), width=.25) +
    geom_bar(stat="identity", position="dodge", colour="black") + 
    facet_grid(.~inTAD) +
    geom_text(aes(label=signif(avgPercent, 2), y= avgPercent), position=position_dodge(width=1), hjust=-0.25, size=5, angle=90) +
    scale_fill_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + labs(y="Gene pairs [%]", x="")
ggsave(p, file=paste0(outPrefix, ".10_1000kb.PPI.by_TAD_and_group.barplot.pdf"), w=3.5, h=7)

# PPI by TAD and group [10kb - 1000kb] only non NA
p <- ggplot(freqDF[freqDF$PPI == "PPI" & !is.na(freqDF$PPI),], aes(x=group, y=avgPercent, fill=group)) +
    geom_errorbar(aes(ymax = avgPercent + sd , ymin=avgPercent - sd), position=position_dodge(width=0.9), width=.25) +
    geom_bar(stat="identity", position="dodge", colour="black") + 
    facet_grid(.~inTAD) +
    geom_text(aes(label=signif(avgPercent, 2), y= avgPercent), position=position_dodge(width=1), vjust=1.25, size=5) +
    scale_fill_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + labs(y="Gene pairs with PPI  [%]", x="")
ggsave(p, file=paste0(outPrefix, ".10_1000kb.PPI_noNA.by_TAD_and_group.barplot.pdf"), w=3.5, h=7)

#-----------------------------------------------------------------------
# PPI by TAD and group for close pairs
#-----------------------------------------------------------------------
#~ subDF <- subset(allDF, sSamp & sTAD & sExp & sSpec & distTadBin == "10-1000kb")
subDF <- subset(allDF, sSamp & sTAD & sExp & sSpec & distGroup == "close")

dc <- ddply(subDF, .(group, inTAD, replicate, PPI), summarize, count=length(PPI), .drop=FALSE)
summedCounts <- rep(ddply(dc, .(group, inTAD, replicate), summarize, s=sum(count))$s, each=length(levels(dc$PPI))+1)
freqRepDF <- cbind(dc, cbind(summedCounts, percent=dc$count * 100 / summedCounts))
freqDF <- ddply(freqRepDF, .(group, inTAD, PPI), summarize, avg=mean(count), avgPercent=mean(percent,na.rm=TRUE), sd=sd(percent,na.rm=TRUE))

pvalDF <- ddply(subDF, .(group), summarize, p=fisher.test(PPI=="PPI", inTAD)$p.value)

p <- ggplot(freqDF, aes(x=PPI, y=avgPercent, fill=group)) +
    geom_errorbar(aes(ymax = avgPercent + sd , ymin=avgPercent - sd), position=position_dodge(width=0.9), width=.25) +
    geom_bar(stat="identity", position="dodge", colour="black") + 
    facet_grid(.~inTAD) +
    geom_text(aes(label=signif(avgPercent, 2), y= avgPercent), position=position_dodge(width=1), hjust=-0.25, size=5, angle=90) +
    scale_fill_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + labs(y="Gene pairs [%]", x="")
ggsave(p, file=paste0(outPrefix, "close.PPI.by_TAD_and_group.barplot.pdf"), w=3.5, h=7)

# PPI by TAD and group close only non NA
p <- ggplot(freqDF[freqDF$PPI == "PPI" & !is.na(freqDF$PPI),], aes(x=group, y=avgPercent, fill=group)) +
    geom_errorbar(aes(ymax = avgPercent + sd , ymin=avgPercent - sd), position=position_dodge(width=0.9), width=.25) +
    geom_bar(stat="identity", position="dodge", colour="black") + 
    facet_grid(.~inTAD) +
    geom_text(aes(label=signif(avgPercent, 2), y= avgPercent), position=position_dodge(width=1), vjust=1.25, size=5) +
    scale_fill_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + labs(y="Gene pairs with PPI  [%]", x="")
ggsave(p, file=paste0(outPrefix, ".close.PPI_noNA.by_TAD_and_group.barplot.pdf"), w=3.5, h=7)


# PPI by TAD and group close only non NA
p <- ggplot(freqDF[freqDF$PPI == "PPI" & !is.na(freqDF$PPI),], aes(x=inTAD, y=avgPercent, fill=group)) +
    geom_errorbar(aes(ymax = avgPercent + sd , ymin=avgPercent - sd), position=position_dodge(width=0.9), width=.25) +
    geom_bar(stat="identity", position="dodge", colour="black") + 
    facet_grid(.~group) +
    geom_text(aes(label=signif(avgPercent, 2), y= avgPercent), position=position_dodge(width=1), vjust=1.25, size=5) +
    scale_fill_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + labs(y="Gene pairs with PPI  [%]", x="")
ggsave(p, file=paste0(outPrefix, ".close.PPI_noNA.by_group_and_inTAD.barplot.pdf"), w=3.5, h=7)


#-----------------------------------------------------------------------
# Any PPI by TAD and group for close pairs
#-----------------------------------------------------------------------
subDF <- subset(allDF, sSamp & sTAD & sExp & sSpec & distGroup == "close")

# redefine PPI by allowing any score in HIPPIE
subDF$anyPPI <- factor(!is.na(subDF$HIPPIE), levels=c(TRUE, FALSE), labels=c("PPI", "no PPI"))

dc <- ddply(subDF, .(group, inTAD, replicate), summarize, n=length(anyPPI), count=sum(anyPPI=="PPI"), percent=sum(anyPPI=="PPI")/length(anyPPI)*100)

freqDF <- ddply(dc, .(group, inTAD), summarize, 
    avgCount=mean(count, na.rm=TRUE), 
    sdCount=sd(count, na.rm=TRUE), 
    avgPercent=mean(percent), 
    sdPercent=sd(percent), 
    avgN=mean(n), 
    sdN=sd(n)
)

# calculate p-values

pvals <- sapply(levels(subDF$group), function(gr) fisher.test(table(subDF[subDF$group == gr,"anyPPI"], subDF[subDF$group == gr,"inTAD"]))$p.value) 

pvalDF <- ddply(subDF, .(group), summarize, 
    p=fisher.test(anyPPI=="PPI", inTAD)$p.value,
    ypos = 1.1*max(freqDF$avgCount),
    yposPercent = 1.1*max(freqDF$avgPercent),
    inTAD=NA
)

# anyPPI by TAD and group close only non NA
p <- ggplot(freqDF, aes(x=inTAD, y=avgCount, fill=group)) +
    geom_errorbar(aes(ymax = avgCount + sdCount , ymin=avgCount - sdCount), position=position_dodge(width=0.9), width=.25) +
    geom_bar(stat="identity", position="dodge", colour="black") + 
    facet_grid(.~group) +
    geom_text(aes(label=signif(avgCount, 2)), position=position_dodge(width=1), vjust=1.25, size=5) +
    geom_text(aes(label=paste0("p=", signif(p,3)), x=1.5, y=ypos), size=5, data=pvalDF) + 
    scale_fill_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="none") + labs(y="Gene pairs with PPI", x="")
ggsave(p, file=paste0(outPrefix, ".close.anyPPI_noNA.by_group_and_inTAD.barplot.pdf"), w=3.5, h=7)


# anyPPI by TAD and group close only non NA
p <- ggplot(freqDF, aes(x=inTAD, y=avgPercent, fill=group)) +
    geom_errorbar(aes(ymax = avgPercent + sdPercent , ymin=avgPercent - sdPercent), position=position_dodge(width=0.9), width=.25) +
    geom_bar(stat="identity", position="dodge", colour="black") + 
    facet_grid(.~group) +
    geom_text(aes(label=signif(avgPercent, 2)), position=position_dodge(width=1), vjust=1.25, size=5) +
    geom_text(aes(label=paste0("p=", signif(p,3)), x=1.5, y=yposPercent), size=5, data=pvalDF) + 
    scale_fill_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="none") + labs(y="Gene pairs with PPI [%]", x="")
ggsave(p, file=paste0(outPrefix, ".close.anyPPI_noNA_percent.by_group_and_inTAD.barplot.pdf"), w=3.5, h=7)

#-----------------------------------------------------------------------
# PPI by subTAD and group [10kb - 1000kb]
#-----------------------------------------------------------------------
subDF <- allDF[sSamp & sTAD & sExp & sSpec & allDF$distTadBin == "10-1000kb",]

dc <- ddply(subDF, .(group, subTAD, replicate, PPI), summarize, count=length(PPI), .drop=FALSE)
summedCounts <- rep(ddply(dc, .(group, subTAD, replicate), summarize, s=sum(count))$s, each=length(levels(dc$PPI))+1)
freqRepDF <- cbind(dc, cbind(summedCounts, percent=dc$count * 100 / summedCounts))
freqDF <- ddply(freqRepDF, .(group, subTAD, PPI), summarize, avg=sum(count), avgPercent=mean(percent,na.rm=TRUE), sd=sd(percent,na.rm=TRUE))

p <- ggplot(freqDF, aes(x=PPI, y=avgPercent, fill=group)) +
    geom_bar(stat="identity", position="dodge", colour="black") + 
    geom_errorbar(aes(ymax = avgPercent + sd , ymin=avgPercent - sd), position=position_dodge(width=0.9), width=.25) +
    facet_grid(.~subTAD) +
    geom_text(aes(label=signif(avgPercent, 2), y= avgPercent), position=position_dodge(width=1), hjust=-0.25, size=5, angle=90) +
    scale_fill_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + labs(y="Gene pairs [%]", x="")
ggsave(p, file=paste0(outPrefix, ".10_1000kb.PPI.by_subTAD_and_group.barplot.pdf"), w=7, h=7)

# PPI by subTAD and group [10kb - 1000kb] only non NA
p <- ggplot(freqDF[freqDF$PPI == "PPI" & !is.na(freqDF$PPI),], aes(x=PPI, y=avgPercent, fill=group)) +
    geom_bar(stat="identity", position="dodge", colour="black") + 
    geom_errorbar(aes(ymax = avgPercent + sd , ymin=avgPercent - sd), position=position_dodge(width=0.9), width=.25) +
    facet_grid(.~subTAD) +
    geom_text(aes(label=signif(avgPercent, 2), y= avgPercent), position=position_dodge(width=1), vjust=1.25, size=5) +
    scale_fill_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + labs(y="Gene pairs with PPI [%]", x="")
ggsave(p, file=paste0(outPrefix, ".10_1000kb.PPI_noNA.by_subTAD_and_group.barplot.pdf"), w=7, h=7)

#-----------------------------------------------------------------------
# PPI by subTAD and distTadBin
#-----------------------------------------------------------------------
subDF <- allDF[sSamp & sTAD & sExp & sSpec,]

dc <- ddply(subDF, .(group, subTAD, distTadBin, replicate, PPI), summarize, count=length(PPI), .drop=FALSE)
summedCounts <- rep(ddply(dc, .(group, subTAD, distTadBin, replicate), summarize, s=sum(count))$s, each=length(levels(dc$PPI))+1)
freqRepDF <- cbind(dc, cbind(summedCounts, percent=dc$count * 100 / summedCounts))
freqDF <- ddply(freqRepDF, .(group, subTAD, distTadBin, PPI), summarize, avg=sum(count), avgPercent=mean(percent,na.rm=TRUE), sd=sd(percent,na.rm=TRUE))

p <- ggplot(freqDF, aes(x=PPI, y=avgPercent, fill=group)) +
    geom_bar(stat="identity", position="dodge", colour="black") + 
    geom_errorbar(aes(ymax = avgPercent + sd , ymin=avgPercent - sd), position=position_dodge(width=0.9), width=.25) +
    facet_grid(.~distTadBin*subTAD) +
    geom_text(aes(label=signif(avgPercent, 2), y= avgPercent), position=position_dodge(width=1), hjust=-0.25, size=5, angle=90) +
    scale_fill_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + labs(y="Gene pairs [%]", x="")
ggsave(p, file=paste0(outPrefix, ".PPI.by_subTAD_and_group.barplot.pdf"), w=14, h=7)

# PPI by subTAD and group [10kb - 1000kb] only non NA
p <- ggplot(freqDF[freqDF$PPI == "PPI" & !is.na(freqDF$PPI),], aes(x=PPI, y=avgPercent, fill=group)) +
    geom_bar(stat="identity", position="dodge", colour="black") + 
    geom_errorbar(aes(ymax = avgPercent + sd , ymin=avgPercent - sd), position=position_dodge(width=0.9), width=.25) +
    facet_grid(.~distTadBin*subTAD) +
    geom_text(aes(label=signif(avgPercent, 2), y= avgPercent), position=position_dodge(width=1), vjust=1.25, size=5) +
    scale_fill_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + labs(y="Gene pairs with PPI [%]", x="")
ggsave(p, file=paste0(outPrefix, ".PPI_noNA.by_subTAD_and_group.barplot.pdf"), w=7, h=7)

#-----------------------------------------------------------------------
# plot contact fequencies against different features
#-----------------------------------------------------------------------

# create a list of subset for use in loop
subDFlist <- list(
"close"=subset(allDF, sSamp & sTAD & sExp & sSpec & allDF$distGroup == "close"),
"10-1000kb"=subset(allDF, sSamp & sTAD & sExp & sSpec & allDF$distTadBin == "10-1000kb"),
"distal"=subset(allDF, sSamp & sTAD & sExp & sSpec & allDF$distGroup == "distal"),
"allCis"=subset(allDF, sSamp & sTAD & sExp & sSpec)
)

# iterate over all contact measurements (in HiCcolumns):
#~ for (HiCcol in HiCcolumns){
for (HiCcol in c("HiC", "HiCNoZero", "captureC_raw", "captureC_rawNoZero")){
    
    message(paste("INFO: Plotting:", HiCcol))
    HiClab <- HiClabels[HiCcol]
    
    subDF <- allDF[sSamp & sTAD & sExp & sSpec,]

    # contacts by TAD and TAD distance bins
    nDF <- data.frame(table(subDF[, c("group", "distTadBin", "inTAD")], useNA="ifany"))
    p = ggplot(subDF, aes_string(x="inTAD", y=HiCcol)) + 
        scale_y_log10() + annotation_logticks(sides="l") +
        geom_boxplot(aes(colour = group), lwd=1.5) + 
        facet_grid(.~distTadBin) + 
        scale_color_manual(values=COL, guide_legend(title = "")) +
        theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + 
        guides(fill=guide_legend(title="")) +
        labs(y=HiClab, x="Linear distance bin [kb]") + 
        geom_text(aes(label=paste0("n=",Freq), vjust=(as.numeric(group)-1)*2-1,  y=10^-1, angle=90), data=nDF)
    ggsave(p, file=paste0(outPrefix, ".", HiCcol, ".by_TADdist_byTAD.boxplot.pdf"), w=7, h=7)
    
    
    #------------------------------------------------------------------------
    # iterate over close, distal, and all:
    for(subName in names(subDFlist)){
        
        subDF <- subDFlist[[subName]]
        
        # contacts by group        
        xlabels <- paste0(c("Paralogs", "Sampled"), 
            "\n n=", applyToSubset(subDF, function(v) sum(!is.na(v)), HiCcol, "group"), 
            "\n med=", signif(applyToSubset(subDF, median, HiCcol, "group", na.rm=TRUE), 3),
            "\n avg=", signif(applyToSubset(subDF, mean, HiCcol, "group", na.rm=TRUE), 3)
            )
        ws.test <- wilcox.test(subDF[,HiCcol] ~ subDF[,"group"])
        
        p = ggplot(subDF, aes_string(x="group", y=HiCcol, color="group")) + 
            geom_boxplot(lwd=1.5) + scale_y_log10() + annotation_logticks(sides="l") +
            scale_color_manual(values=COL, guide_legend(title = "")) +
            theme_bw() + theme(text = element_text(size=20), legend.position="none", axis.text.x=element_text(angle = 45, hjust = 1)) + 
            guides(fill=guide_legend(title="")) +
            labs(y=HiClab, x="", title=paste0("p = ", signif(ws.test$p.value, 3))) 
        ggsave(p, file=paste0(outPrefix, ".", subName, ".", HiCcol, ".ggboxplot.noLables.pdf"), width=3.5)    
		p <- p + scale_x_discrete(labels=xlabels) + theme(axis.text.x=element_text(angle = 0, hjust = .5))
        ggsave(p, file=paste0(outPrefix, ".", subName, ".", HiCcol, ".ggboxplot.withLables.pdf"), width=3.5)    
    
        # contacts by group but with at least 100kb distance pairs:
        for (D in c(10, 100)){
            subDFMinDist <- subDF[subDF$dist >= D,]
            xlabels=paste0(c("Paralogs", "Sampled"), 
                "\n n=", applyToSubset(subDFMinDist, function(v) sum(!is.na(v)), HiCcol, "group"), 
                "\n med=", signif(applyToSubset(subDFMinDist, median, HiCcol, "group", na.rm=TRUE), 3),
                "\n avg=", signif(applyToSubset(subDFMinDist, mean, HiCcol, "group", na.rm=TRUE), 3)
                )
            ws.test = wilcox.test(subDFMinDist[,HiCcol] ~ subDFMinDist[,"group"])
        
            p = ggplot(subDFMinDist, aes_string(x="group", y=HiCcol, color="group")) + 
                geom_boxplot(lwd=1.5)  + scale_y_log10() + annotation_logticks(sides="l") +
                scale_color_manual(values=COL, guide_legend(title = "")) +
                theme_bw() + theme(text = element_text(size=20), legend.position="none") + 
                guides(fill=guide_legend(title="")) +
                labs(y=HiClab, x="", title=paste0("p = ", signif(ws.test$p.value, 3))) +
                scale_x_discrete(labels=xlabels )
            ggsave(p, file=paste0(outPrefix, ".close.min", D, "kb.", HiCcol, ".ggboxplot.pdf"), width=3.5)    
        }
        
        # contacts by group and distance bin
        nDF <- data.frame(table(subDF[, c("group", "distBin")]))
        p = ggplot(subDF, aes_string(x="group", y=HiCcol)) + 
            scale_y_log10() + annotation_logticks(sides="l") +
            geom_boxplot(aes(x=group, colour = group), lwd=1.5) + 
            #geom_jitter(aes(x=group, colour = group), alpha=0.1, width=.5) + 
            facet_grid(.~distBin) + 
            scale_color_manual(values=COL, guide_legend(title = "")) +
            theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="none") + 
            guides(fill=guide_legend(title="")) +
            labs(y=HiClab, x="Linear distance bin [kb]") + 
            geom_text(aes(label=paste0("n=",Freq), y=max(subDF[,HiCcol], na.rm=TRUE)), data=nDF)
            
        ggsave(p, file=paste0(outPrefix, ".", subName, ".", HiCcol, ".by_dist.boxplot.pdf"), w=7)
    
        # contacts by group, inTAD and distance bin
        p = ggplot(subDF, aes_string(x="inTAD", y=HiCcol)) + 
            scale_y_log10() + annotation_logticks(sides="l") +
            geom_boxplot(aes(colour = group), lwd=1.5) + 
            facet_grid(.~distBin) + 
            scale_color_manual(values=COL, guide_legend(title = "")) +
            theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + 
            guides(fill=guide_legend(title="")) +
            labs(y=HiClab, x="Linear distance bin [kb]")
        ggsave(p, file=paste0(outPrefix, ".", subName, ".", HiCcol, ".by_dist_byTAD.boxplot.pdf"), w=7, h=7)
    
        # contacts by group, strand and distance bin    
        nDF <- data.frame(table(subDF[, c("group", "distBin", "sameStrand")]))
        p = ggplot(subDF, aes_string(x="sameStrand", y=HiCcol)) + 
            scale_y_log10() + annotation_logticks(sides="l") +
            geom_boxplot(aes(colour = group), lwd=1.5) + 
            facet_grid(.~distBin) + 
            scale_color_manual(values=COL, guide_legend(title = "")) +
            theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + 
            guides(fill=guide_legend(title="")) +
            labs(y=HiClab, x="Linear distance bin [kb]") +
            geom_text(aes(label=paste0("n=",Freq), vjust=(as.numeric(group)-1)*2-1,  y=max(subDF[,HiCcol], na.rm=TRUE), hadjust=1, angle=90), data=nDF)
        ggsave(p, file=paste0(outPrefix, ".", subName, ".", HiCcol, ".by_dist_byStrand.boxplot.pdf"), w=7, h=7)
        
        # contacts by group and expression bin
        #nDF <- data.frame(table(subDF[, c("group", "expBin")], useNA="ifany"))
        nDF <- ddply(subDF, .(expBin, group), summarize, n=sum(!is.na(get(HiCcol))))
        pvalDF <- ddply(subDF, .(expBin), summarize, p=ifelse(min(table(!is.na(get(HiCcol)), group))>1, wilcox.test(get(HiCcol) ~ group)$p.value, NA))
        
        p = ggplot(subDF, aes_string(x="group", y=HiCcol)) + 
            scale_y_log10() + annotation_logticks(sides="l") +
            geom_boxplot(aes(x=group, colour = group), lwd=1.5) + 
            #geom_jitter(aes(x=group, colour = group), alpha=0.1, width=.5) + 
            facet_grid(.~expBin) + 
            scale_color_manual(values=COL, guide_legend(title = "")) +
            theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="none") + 
            guides(fill=guide_legend(title="")) +
            labs(y=HiClab, x="Average Expression") + 
            geom_text(aes(label=paste0("p=",signif(p, 3)), x=1.5, y=max(subDF[,HiCcol], na.rm=TRUE)), size=5, data=pvalDF) +
            geom_text(aes(label=paste0("n=",n), y=0.75*min(subDF[,HiCcol], na.rm=TRUE)), data=nDF)
        ggsave(p, file=paste0(outPrefix, ".", subName, ".", HiCcol, ".by_expBin.boxplot.pdf"))
    
        # contacts by group and average expression as dotplot
        p <- dotplotWithDensityLogXY(revDF(subDF), "avgExp+1", HiCcol, "group", COL, fit=TRUE, ALPHA=.2, ylab=HiClab, xlab="Mean expression of gene pair")
        pdf(paste0(outPrefix, ".", subName, ".", HiCcol, ".vs_avgExp.dotplot.pdf"))
            grid::grid.draw(p)
        dev.off()
    
        # contacts by expression correlation bins
        nDF <- data.frame(table(subDF[, c("group", "expCorBin")], useNA="ifany"))
        p = ggplot(subDF, aes_string(x="group", y=HiCcol)) + 
            scale_y_log10() + annotation_logticks(sides="l") +
            geom_boxplot(aes(x=group, colour = group), lwd=1.5) + 
            facet_grid(.~expCorBin) + 
            scale_color_manual(values=COL, guide_legend(title = "")) +
            theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="none") + 
            guides(fill=guide_legend(title="")) +
            labs(y=HiClab, x="Expression Correlation") + 
            geom_text(aes(label=paste0("n=",Freq), y=max(subDF[,HiCcol], na.rm=TRUE)), data=nDF)
        ggsave(p, file=paste0(outPrefix, ".", subName, ".", HiCcol, ".by_expCorBin.boxplot.pdf"))
    
        p <- dotplotWithDensityLogXY(revDF(subDF), "expCor", HiCcol, "group", COL, fit=TRUE, ALPHA=.2, xlog=FALSE, ylog=TRUE, ylab=HiClab, xlab="Correlation of gene expression")
        pdf(paste0(outPrefix, ".", subName, ".", HiCcol, ".vs_expCor.dotplot.pdf"))
            grid::grid.draw(p)
        dev.off()

        # contacts by group and distance (dotplot)
        distSelcet <- subDF$dist >= 1 & subDF$dist <= 10^5
        p <- ggplot(revDF(subDF[distSelcet,]), aes(x=dist, y=get(HiCcol), color=group, fill=group)) + 
            geom_point(alpha=.5, size=.5, shape=20)  +
            scale_y_log10() + 
            scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +  
            annotation_logticks(sides="bl") +
            theme_bw() + theme(text = element_text(size=20), legend.position="bottom") +
            scale_color_manual(values = COL, guide_legend(title = "")) + scale_fill_manual(values = COL, guide_legend(title = "")) +
            xlab("Distance [kb]") + ylab(HiClab)
        
        p + geom_smooth(alpha=.5, fullrange=TRUE) + 
        ggsave(paste0(outPrefix, ".", subName, ".", HiCcol, ".vs_dist.dotplot_smoothFit.pdf"))
        
        p + geom_smooth(method="lm", fullrange=TRUE, alpha=.5) + 
        ggsave(paste0(outPrefix, ".", subName, ".", HiCcol, ".vs_dist.dotplot_linearFit.pdf"))

        # contacts by distance, group, and inTAD (dotplot)

        p <- ggplot(revDF(subDF[distSelcet,]), aes(x=dist, y=get(HiCcol), color=group, fill=group)) + 
            geom_point(alpha=.5, size=.5, shape=20)  +
            facet_grid(.~inTAD, scales="free_x") + 
            scale_y_log10() + 
            scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +  
            annotation_logticks(sides="bl") +
            theme_bw() + theme(text = element_text(size=20), legend.position="bottom") +
            scale_color_manual(values = COL, guide_legend(title = "")) + scale_fill_manual(values = COL, guide_legend(title = "")) +
            xlab("Distance [kb]") + ylab(HiClab)

        p + geom_smooth(alpha=.5, fullrange=TRUE) + 
        ggsave(paste0(outPrefix, ".", subName, ".", HiCcol, ".vs_dist_by_inTAD.dotplot_smoothFit.pdf"))
        
        p + geom_smooth(method="lm", fullrange=TRUE, alpha=.5) + 
        ggsave(paste0(outPrefix, ".", subName, ".", HiCcol, ".vs_dist_by_inTAD.dotplot_linearFit.pdf"))
    
    }    
}


#------------------------------------------------------------------------
# plot Hi-C contacts versus genomic distance
#-----------------------------------------------------------------------

#------------------------------------------------------------------------
# iterate over close, distal, and all:
for(subName in names(subDFlist)){
    
    subDF <- subDFlist[[subName]]
        
    p0a <- dotplotWithDensityLogXY(revDF(subDF), "dist", "HiCRaw", "group", COL, fit=TRUE, ALPHA=.2)
    p0b <- dotplotWithDensityLogXY(revDF(subDF), "dist", "HiCRawNoZero", "group", COL, fit=TRUE, ALPHA=.2)
    p1 <- dotplotWithDensityLogXY(revDF(subDF), "dist", "HiC", "group", COL, fit=TRUE, ALPHA=.2)
    p2 <- dotplotWithDensityLogXY(revDF(subDF), "dist", "HiCobsExp", "group", COL, fit=TRUE, ALPHA=.2)
    p3 <- dotplotWithDensityLogXY(revDF(subDF), "dist", "captureC_raw", "group", COL, fit=TRUE, ALPHA=.2)
    p4 <- dotplotWithDensityLogXY(revDF(subDF), "dist", "captureC_ObsExp", "group", COL, fit=TRUE, ALPHA=.2)
    
    pdf(paste0(outPrefix, ".", subName, ".Hi-C_vs_dist_all.ggboxplot.pdf"), w=12, h=18)
        do.call(grid.arrange, list(p0a, p0b, p1,p2,p3,p4)) 
    dev.off()

}

#------------------------------------------------------------------------
# plot expression correlation vs. genomic distance
#-----------------------------------------------------------------------
subDF <- allDF[sSamp & sTAD & sExp & sSpec,]

p <- dotplotWithDensityLogXY(revDF(subDF), "dist", "expCor", "group", COL, fit=TRUE, ALPHA=.2, xlog=TRUE, ylog=FALSE)

pdf(paste0(outPrefix, ".expCor_vs_dist.dotplot.pdf"))
    grid::grid.draw(p)
dev.off()

#=======================================================================
# Compare one2one orthologs of human paralogs in mouse and dog
#=======================================================================

#-----------------------------------------------------------------------
# Ortholog_TAD by group and species
#-----------------------------------------------------------------------
# select with ortholog data from both species but only close distance
subDF <- allDF[sSamp & sTAD & sExp & allDF$distGroup == "close" & allDF$ortholog_one2one,]

# count gene pairs by combinations
dc <- ddply(subDF, .(group, replicate, speciesFactor, ortholog_TAD), summarize, count=length(ortholog_TAD))

# sum up counts for sub-combination to get percentage values
summedCounts <- rep(ddply(dc, .(group, replicate, speciesFactor), summarize, s=sum(count))$s, each=length(levels(dc$ortholog_TAD)))
# add percentage values and counts by marking non-existing combinations with NA

freqRepDF <- cbind(dc, cbind(summedCounts, dc$count), percent=dc$count * 100 / summedCounts)

# combine replicates by taking average and sd of counts and percentages
freqDF <- ddply(freqRepDF, .(group, speciesFactor, ortholog_TAD), summarize, avgCount=mean(count, na.rm=TRUE), sdCount=sd(count, na.rm=TRUE), avgPercent=mean(percent, na.rm=TRUE), sdPercent=sd(percent, na.rm=TRUE))

# calculate p-values
pvalDF <- ddply(subDF, .(speciesFactor), summarize, p=fisher.test(group, ortholog_TAD)$p.value)

# combine p-values with plotting coordinates (max of values per group)
pvalDF <- merge(pvalDF, ddply(freqDF[freqDF$ortholog_TAD=="same TAD",], .(speciesFactor), summarize, avgCount=max(avgCount)+50, avgPercent=max(avgPercent)+5))
pvalDF[,"group"] <- "paralog"

p <- ggplot(freqDF[freqDF$ortholog_TAD=="same TAD",], aes(x=group, y=avgPercent, fill=group)) +
    facet_grid(.~speciesFactor) +
    geom_errorbar(aes(ymax = avgPercent + sdPercent , ymin=avgPercent - sdPercent), position=position_dodge(width=0.9), width=.25) +
    geom_bar(stat="identity", position="dodge", colour="black") + 
    geom_text(aes(label=round(avgPercent), y= avgPercent), position=position_dodge(width=1), vjust=1.25, size=5) +
    geom_text(aes(x=1.5, label=paste0("p=", signif(p,2))), data=pvalDF, size=5) + 
    scale_fill_manual(values=COL_ORTHO, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="none") + labs(y="Ortholog pairs in same TAD [%]", x="")
ggsave(p, file=paste0(outPrefix, ".close.orthologs.orthologs_TAD.by_group_and_species.barplot.pdf"), w=3.5)


#-----------------------------------------------------------------------
# inTAD AND Ortholog_TAD by group and species
#-----------------------------------------------------------------------
# select with ortholog data from both species but only close distance
subDF <- allDF[sSamp & sTAD & sExp & allDF$distGroup == "close" & allDF$ortholog_one2one & allDF$inTAD == "same TAD",]

# count gene pairs by combinations
dc <- ddply(subDF, .(group, replicate, speciesFactor, ortholog_TAD), summarize, count=length(ortholog_TAD))

# sum up counts for sub-combination to get percentage values
summedCounts <- rep(ddply(dc, .(group, replicate, speciesFactor), summarize, s=sum(count))$s, each=length(levels(dc$ortholog_TAD)))
# add percentage values and counts by marking non-existing combinations with NA

freqRepDF <- cbind(dc, cbind(summedCounts, dc$count), percent=dc$count * 100 / summedCounts)

# combine replicates by taking average and sd of counts and percentages
freqDF <- ddply(freqRepDF, .(group, speciesFactor, ortholog_TAD), summarize, avgCount=mean(count, na.rm=TRUE), sdCount=sd(count, na.rm=TRUE), avgPercent=mean(percent, na.rm=TRUE), sdPercent=sd(percent, na.rm=TRUE))

# calculate p-values
pvalDF <- ddply(subDF, .(speciesFactor), summarize, p=fisher.test(group, ortholog_TAD)$p.value)

# combine p-values with plotting coordinates (max of values per group)
pvalDF <- merge(pvalDF, ddply(freqDF[freqDF$ortholog_TAD=="same TAD",], .(speciesFactor), summarize, avgCount=max(avgCount)+50, avgPercent=max(avgPercent)+5))
pvalDF[,"group"] <- "paralog"


p <- ggplot(freqDF[freqDF$ortholog_TAD=="same TAD",], aes(x=group, y=avgPercent, fill=group)) +
    facet_grid(.~speciesFactor) +
    geom_errorbar(aes(ymax = avgPercent + sdPercent , ymin=avgPercent - sdPercent), position=position_dodge(width=0.9), width=.25) +
    geom_bar(stat="identity", position="dodge", colour="black") + 
    geom_text(aes(label=round(avgPercent), y= avgPercent), position=position_dodge(width=1), vjust=1.25, size=5) +
    geom_text(aes(x=1.5, label=paste0("p=", signif(p,2))), data=pvalDF, size=5) + 
    scale_fill_manual(values=COL_ORTHO, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="none") + labs(y="Pairs conserved in same TAD [%]", x="")
ggsave(p, file=paste0(outPrefix, ".close.orthologs.orthologs_conserved_TAD.by_group_and_species.barplot.pdf"), w=3.5)

#-----------------------------------------------------------------------
# Ortholog contacts by group and species
#-----------------------------------------------------------------------

# create a list of subset for use in loop with both species included
subDFspeciesList <- list(
"close"=allDF[sSamp & sTAD & sExp & allDF$distGroup == "close",],
"distal"=allDF[sSamp & sTAD & sExp & allDF$distGroup == "distal",],
"allCis"=allDF[sSamp & sTAD & sExp ,]
)


for (HiCcol in c("ortholog_HiC", "ortholog_HiCnorm", "ortholog_HiCNoZero", "ortholog_HiCnormNoZero")){
    
    message(paste("INFO: Plotting:", HiCcol))

    subDF <- allDF[sSamp & sTAD & sExp & allDF$ortholog_one2one,]

    #-----------------------------------------------------------------------
    # Ortholog contacts by group human TAD and human TAD distance bins
    #-----------------------------------------------------------------------
    nDF <- data.frame(table(subDF[, c("group", "speciesFactor", "distTadBin", "inTAD")], useNA="ifany"))
    p = ggplot(subDF, aes_string(x="inTAD", y=HiCcol)) + 
        scale_y_log10() + annotation_logticks(sides="l") +
        geom_boxplot(aes(colour = group), lwd=1.5) + 
        facet_grid(.~distTadBin*speciesFactor) + 
        scale_color_manual(values=COL_ORTHO, guide_legend(title = "")) +
        theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + 
        guides(fill=guide_legend(title="")) +
        labs(y=HiCcol, x="Human TAD") + 
        geom_text(aes(label=paste0("n=",Freq), vjust=(as.numeric(group)-1)*2-1,  y=10^-1, angle=90), data=nDF)
    ggsave(p, file=paste0(outPrefix, ".", HiCcol, ".by_TADdist_bySpecies_byTAD.boxplot.pdf"), w=14, h=7)
    
    #-----------------------------------------------------------------------
    # Ortholog contacts by group ortholog TAD and human TAD distance bins
    #-----------------------------------------------------------------------
    nDF <- data.frame(table(subDF[, c("group", "speciesFactor", "distTadBin", "ortholog_TAD")], useNA="ifany"))
    p = ggplot(subDF, aes_string(x="ortholog_TAD", y=HiCcol)) + 
        scale_y_log10() + annotation_logticks(sides="l") +
        geom_boxplot(aes(colour = group), lwd=1.5) + 
        facet_grid(.~distTadBin*speciesFactor) + 
        scale_color_manual(values=COL_ORTHO, guide_legend(title = "")) +
        theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + 
        guides(fill=guide_legend(title="")) +
        labs(y=HiCcol, x="Ortholog TAD") + 
        geom_text(aes(label=paste0("n=",Freq), vjust=(as.numeric(group)-1)*2-1,  y=10^-1, angle=90), data=nDF)
    ggsave(p, file=paste0(outPrefix, ".", HiCcol, ".by_TADdist_bySpecies_byTAD.boxplot.pdf"), w=14, h=7)    

    #------------------------------------------------------------------------
    # iterate over close, distal, and all:
    for(subName in names(subDFspeciesList)){
        
        subDF <- subDFspeciesList[[subName]]
        
        # select for those pairs with one-to-one orghologs
        subDF <- subDF[subDF$ortholog_one2one,]
    
        #-------------------------------------------------------------------
        # Compare Hi-C counts of close orthologs      
        #------------------------------------------------------------------------

        summaryDF <- ddply(subDF, .(speciesFactor, group), summarize, n=sum(!is.na(get(HiCcol))), med=median(get(HiCcol), na.rm=TRUE), avg=mean(get(HiCcol), na.rm=TRUE))
        
        pvalDF <- ddply(subDF, .(speciesFactor), summarize, p=wilcox.test(get(HiCcol) ~ group)$p.value)
        pvalDF$group="paralog"
    
        p = ggplot(subDF, aes(x=group, y=get(HiCcol))) +
            facet_grid(.~speciesFactor) +
            geom_boxplot(aes(colour=group), lwd=1.5) + scale_y_log10() + annotation_logticks(sides="l") +
            scale_color_manual(values=COL_ORTHO) +
            theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="none") + 
            labs(y="Hi-C contacts of orthologs", x="")  +
            geom_text(aes(label=paste0("p=",signif(p,2)), x=1.5, y=max(subDF[,HiCcol], na.rm=TRUE)), size=5, data=pvalDF)
            
        
        ggsave(p, file=paste0(outPrefix, ".", subName, ".", HiCcol, ".by_group.boxplot.withLables.pdf"), width=3.5, height=7)
        p <- p + geom_text(aes(label=paste0("n=",n, "\nmed=", signif(med,3), "\navg=", signif(avg,3)),  y=min(subDF[,HiCcol], na.rm=TRUE)), vjust = -.1, data=summaryDF)
        ggsave(p, file=paste0(outPrefix, ".", subName, ".", HiCcol, ".by_group.boxplot.noLables.pdf"), width=3.5, height=7)
    
    }
    
}


# iterate over all species:
for (orgName in orgStr2Name){
    
    #-------------------------------------------------------------------
    # number of orthologs by group
    #-------------------------------------------------------------------
    # select with ortholog data from the species
    subDF <- allDF[sSamp & sTAD & sExp & allDF$speciesFactor == orgName,]

    dc <- ddply(subDF, .(group, replicate), summarize, n=length(ortholog_one2one), count=sum(ortholog_one2one), percent=sum(ortholog_one2one)/length(ortholog_one2one)*100)
    freqDF <- ddply(dc, .(group), summarize, 
        avgCount=mean(count, na.rm=TRUE), 
        sdCount=sd(count, na.rm=TRUE), 
        avgPercent=mean(percent), 
        sdPercent=sd(percent), 
        avgN=mean(n), 
        sdN=sd(n)
    )
    
    # calculate p-values
    pvals <- fisher.test(subDF$ortholog_one2one, subDF$group)$p.value 
    
    pvalDF <- data.frame(
        p = pvals,
        ypos = max(freqDF$avgCount)+10,
        yposPercent = max(freqDF$avgPercent)+2,
        group=NA
    )

    p <- ggplot(freqDF, aes(x=group, y=avgPercent, fill=group)) +
        geom_errorbar(aes(ymax = avgPercent + sdPercent , ymin=avgPercent - sdPercent), width=.25) +
        geom_bar(stat="identity", colour="black") + 
        scale_fill_manual(values=COL, guide_legend(title = "")) +
        theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="none") + labs(y="Pairs with one-to-one orthologs [%]", x="") + 
        geom_text(aes(label=signif(avgPercent,3), y= avgPercent), vjust=1.25, size=5) +
        geom_text(aes(label=paste0("p=",signif(p,2)), y=yposPercent, x=1.5), size=5, data=pvalDF)
    ggsave(paste0(outPrefix, ".", orgName, ".One2One_orthologsPercent.by_group.barplot.pdf"), w=3.5)

    #-------------------------------------------------------------------
    # correlate linear distances of paralogs in human with mouse orthologs
    #-------------------------------------------------------------------
    subDF <- allDF[sSamp & sTAD & sExp & allDF$speciesFactor == orgName,]

    p <- dotplotWithDensityLogXY(revDF(subDF), "dist", "ortholog_dist", "group", COL, fit=TRUE, ALPHA=.2, xlog=TRUE, ylog=TRUE)
    
    pdf(paste0(outPrefix, ".", orgName, ".dist_vs_orgholog_dist.dotplot.pdf"))
        grid::grid.draw(p)
    dev.off()
}


#-----------------------------------------------------------------------
# Expression correlation by group, exp-data-set, and distGroup
#-----------------------------------------------------------------------
subDF <- subset(allDF, sSamp & sTAD & sSpec)

g <- ggplot(subDF, aes(x=group, y = expCor, fill=group)) + 
    geom_violin(aes(color = group), adjust = .4)  + geom_boxplot(fill=NA, width=.25, lwd=1, outlier.shape = NA) + 
    facet_grid(distGroup~expSource) + scale_fill_manual(values=COL) + scale_color_manual(values=COL) + 
    theme_bw()+ theme(legend.position = "none", axis.text.x=element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=14)) + 
    labs(y="Pearson correlation coefficient", x = "", title= "Co-expression correlation over tissues")  #+ theme(legend.position = "bottom") + guides(fill=guide_legend(title=""))

ggsave(paste0(outPrefix,".expCor.by_distGroup.violinplot.pdf"), g)

#-----------------------------------------------------------------------
# close paralogs Hi-C by inTAD
#-----------------------------------------------------------------------
subDF <- subset(allDF, distGroup=="close" & sSamp & sExp & sSpec)

p <- ggplot(subDF, aes(x=inTAD, fill=group, y=HiCobsExp)) + 
        geom_boxplot() + labs(y="Normalized Hi-C") + scale_y_log10()  + 
        facet_grid(.~tadSource) + theme_bw() + theme(axis.text.x=element_text(angle = 45, hjust = 1), legend.position = "bottom") + theme(text = element_text(size=15)) + scale_fill_manual(values=COL)
ggsave(paste0(outPrefix,".close.HiCobsExp.byTAD_and_tadSource.boxplot.pdf"), w=14,h=7)

p <- ggplot(subDF, aes(x=inTAD, fill=group, y=captureC_ObsExp)) + 
        geom_boxplot() + theme_bw() + labs(y="Captrue HiC (obs/exp)") + scale_y_log10()  + 
        facet_grid(.~tadSource) + theme(axis.text.x=element_text(angle = 45, hjust = 1), legend.position = "bottom") + theme(text = element_text(size=15)) + scale_fill_manual(values=COL) 
ggsave(paste0(outPrefix,".close.captureC_ObsExp.byTAD_and_tadSource.boxplot.pdf"), w=14,h=7)


#-----------------------------------------------------------------------
# common_ortholog by TAD
#-----------------------------------------------------------------------
subDF <- subset(allDF, sSamp & sExp & sSpec & group=="paralog" & distGroup=="close")

p <- ggplot(subDF, aes(x=ortholog_commonOrtholg, fill=inTAD)) +
    geom_bar(aes(y = ..count..), position="dodge") +
    facet_grid(~tadSource, scales="free_y") + theme_bw() + theme(legend.position = "bottom") + theme(text = element_text(size=15)) + scale_fill_manual(values=rev(COL_TAD))
ggsave(paste0(outPrefix,".close.age_commonOrth_by_TAD.barplot.pdf"), w=14,h=7)

#=======================================================================
# Strand and distance
#=======================================================================
#-----------------------------------------------------------------------
# same strand vs. distance in paralogs and sampled
#-----------------------------------------------------------------------
subDF <- subset(allDF, sSamp & sTAD & sExp & sSpec & distGroup=="close")
# add labeld factor 
subDF$sameStrand <- factor(subDF$sameStrand, levels=c(TRUE, FALSE), c("same strand", "opposite strand"))

pPara <- wilcox.test(dist ~ sameStrand, data=subDF[subDF$group=="paralog",])
pSamp <- wilcox.test(dist ~ sameStrand, data=subDF[subDF$group=="sampled",])

pVals <- data.frame(sameStrand=rep(1.5,2), dist=1.1*MAX_DIST/10^3, group=c("paralog", "sampled"), pVal=paste0("p = ", signif(c(pPara$p.value, pSamp$p.value),3))) 

p <- ggplot(subDF, aes(x=sameStrand, y=dist, fill=group)) +
    geom_boxplot() + ylim(0,1.1*MAX_DIST/10^3) + #scale_y_log10() +
    geom_text(data=pVals, aes(label=pVal), size=5) + 
    facet_grid(~group) + theme_bw() + theme(legend.position = "none", axis.text.x=element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=20)) + scale_fill_manual(values=COL) + labs(x="", y="Distance [kb]")
ggsave(paste0(outPrefix,".close.dist_by_sameStrand_and_group.boxplot.pdf"), w=3.5, h=7)


#=======================================================================
# Distance Bins and sub-TAD structure
#=======================================================================

# reduce to sampled by distance and only one source for TAD/exp/species
subDF <- allDF[sSamp & sTAD & sExp & sSpec,]

#-----------------------------
dc <- ddply(subDF, .(group, replicate), summarize, n=length(sameStrand), count=sum(sameStrand), percent=sum(sameStrand)/length(sameStrand)*100)
freqDF <- ddply(dc, .(group), summarize, 
    avgCount=mean(count, na.rm=TRUE), 
    sdCount=sd(count, na.rm=TRUE), 
    avgPercent=mean(percent), 
    sdPercent=sd(percent), 
    avgN=mean(n), 
    sdN=sd(n)
)

# calculate p-values
pVal <- fisher.test(subDF$sameStrand, subDF$group)$p.value

p <- ggplot(freqDF, aes(x=group, y=avgCount, fill=group)) +
    geom_errorbar(aes(ymax = avgCount + sdCount , ymin=avgCount - sdCount), width=.25) +
    geom_bar(stat="identity", colour="black") + 
    scale_fill_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="none") + labs(y="Transcribed from same strand", x="") + 
    geom_text(aes(label=avgCount, y=avgCount), vjust=1.25, size=5) +
    geom_text(aes(label=paste0("p=",signif(pVal,2)), y=1.1*max(freqDF$avgCount), x=1.5), size=5)
    
ggsave(paste0(outPrefix,".gene_pairs.sameStrand_by_group.barplot.pdf"), w=3.5, h=7)

p <- ggplot(freqDF, aes(x=group, y=avgPercent, fill=group)) +
    geom_errorbar(aes(ymax = avgPercent + sdPercent , ymin=avgPercent - sdPercent), width=.25) +
    geom_bar(stat="identity", colour="black") + 
    scale_fill_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="none") + labs(y="Transcribed from same strand [%]", x="") + 
    geom_text(aes(label=signif(avgPercent, 3), y=avgPercent), vjust=1.25, size=5) +
    geom_text(aes(label=paste0("p=",signif(pVal,2)), y=1.1*max(freqDF$avgPercent), x=1.5), size=5)
    
ggsave(paste0(outPrefix,".percent_gene_pairs.sameStrand_by_group.barplot.pdf"), w=3.5, h=7)

#-----------------------------------------------------------------------
# same strand for close pairs by group
#-----------------------------------------------------------------------

# reduce to sampled by distance and only one source for TAD/exp/species and close pairs
subDF <- subset(allDF, sSamp & sTAD & sExp & sSpec & distGroup=="close")

#-----------------------------
dc <- ddply(subDF, .(group, replicate), summarize, n=length(sameStrand), count=sum(sameStrand), percent=sum(sameStrand)/length(sameStrand)*100)
freqDF <- ddply(dc, .(group), summarize, 
    avgCount=mean(count, na.rm=TRUE), 
    sdCount=sd(count, na.rm=TRUE), 
    avgPercent=mean(percent), 
    sdPercent=sd(percent), 
    avgN=mean(n), 
    sdN=sd(n)
)

# calculate p-values
pVal <- fisher.test(subDF$sameStrand, subDF$group)$p.value

p <- ggplot(freqDF, aes(x=group, y=avgCount, fill=group)) +
    geom_errorbar(aes(ymax = avgCount + sdCount , ymin=avgCount - sdCount), width=.25) +
    geom_bar(stat="identity", colour="black") + 
    scale_fill_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="none") + labs(y="Transcribed from same strand", x="") + 
    geom_text(aes(label=avgCount, y=avgCount), vjust=1.25, size=5) +
    geom_text(aes(label=paste0("p=",signif(pVal,2)), y=1.1*max(freqDF$avgCount), x=1.5), size=5)
    
ggsave(paste0(outPrefix,".close.gene_pairs.sameStrand_by_group.barplot.pdf"), w=3.5, h=7)

p <- ggplot(freqDF, aes(x=group, y=avgPercent, fill=group)) +
    geom_errorbar(aes(ymax = avgPercent + sdPercent , ymin=avgPercent - sdPercent), width=.25) +
    geom_bar(stat="identity", colour="black") + 
    scale_fill_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="none") + labs(y="Transcribed from same strand [%]", x="") + 
    geom_text(aes(label=signif(avgPercent, 3), y=avgPercent), vjust=1.25, size=5) +
    geom_text(aes(label=paste0("p=",signif(pVal,2)), y=1.1*max(freqDF$avgPercent), x=1.5), size=5)
    
ggsave(paste0(outPrefix,".close.percent_gene_pairs.sameStrand_by_group.barplot.pdf"), w=3.5, h=7)

#=======================================================================
# In TAD vs. not in TAD 
#=======================================================================

#-----------------------------------------------------------------------
# do paralogs exp correlation vs in same TAD
#-----------------------------------------------------------------------
subDF <- subset(allDF, sSamp & sExp & sSpec & distGroup=="close")

pvalDF <- ddply(subDF, .(group, tadSource), summarize, p=wilcox.test(expCor^2 ~ inTAD)$p.value)

p <- ggplot(subDF, aes(x=inTAD, y=expCor^2)) +
    geom_boxplot(aes(fill=inTAD)) +
    facet_grid(group~tadSource) + theme_bw() + theme(legend.position = "bottom", axis.text.x=element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=15)) + scale_fill_manual(values=rev(COL_TAD)) +
    geom_text(aes(label=paste("p=", signif(p, 3)), y=1.1, x=1.5), data=pvalDF)
    
ggsave(paste0(outPrefix,".close.expCor2_by_inTAD_and_tadSource.boxplot.pdf"), w=14, h=7)

#-----------------------------------------------------------------------
# subset of pairs within same TAD: shared enhancers, expression correlation
#-----------------------------------------------------------------------
# do paralogs exp correlation vs in same TAD
#~ subDF <- subset(allDF, distGroup=="close" & tadSource=="stable_TADs"  & species=="mouse")
subDF <- subset(allDF, sSamp & sTAD & sSpec & distGroup=="close")

p <- ggplot(subDF, aes(x=group, y=expCor, fill=group)) +
    geom_boxplot() +
    facet_grid(inTAD~expSource) + theme_bw() + theme(legend.position = "none", axis.text.x=element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=15)) + scale_fill_manual(values=COL) + xlab("")
ggsave(paste0(outPrefix,".close.expCor_vs_group_and_inTAD.boxplot.pdf"), w=7, h=7)

p <- ggplot(subDF, aes(x=group, y=expCor, fill=group)) +
        geom_violin(aes(color = group), adjust=.5, lwd=1) + geom_boxplot(fill=NA, width=.25, lwd=1) +
        facet_grid(inTAD~expSource) + theme_bw() + theme(legend.position = "none", axis.text.x=element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=14)) + scale_fill_manual(values=COL) + scale_color_manual(values=COL) +
        ylab("Expression Correlation [Pearson R]") + xlab("")
ggsave(paste0(outPrefix,".close.expCor_vs_group_and_inTAD.violinplot.pdf"), w=7, h=7)

#-----------------------------------------------------------------------
# shared enhancer by group and inTAD
#-----------------------------------------------------------------------

subDFList <- list(
"close"=subset(allDF, sSampEh & sTAD & sExp & sSpec & distGroup=="close"),
"10-1000kb"=subset(allDF, sSampEh & sTAD & sExp & sSpec & distTadBin=="10-1000kb")
)

for(subName in names(subDFList)){

    subDF <- subDFList[[subName]]
    
    #-----------------------------
    dc <- ddply(subDF, .(group, inTAD, replicate), summarize, n=length(eh), count=sum(eh), percent=sum(eh)/length(eh)*100)
    freqDF <- ddply(dc, .(inTAD, group), summarize, 
        avgCount=mean(count, na.rm=TRUE), 
        sdCount=sd(count, na.rm=TRUE), 
        avgPercent=mean(percent), 
        sdPercent=sd(percent), 
        avgN=mean(n), 
        sdN=sd(n)
    )
    
    # calculate p-values
    pvalDF <- ddply(subDF, .(inTAD), summarize, p=fisher.test(eh, group)$p.value) 
    pvalDF$ypos <- max(freqDF$avgCount)+10
    pvalDF$yposPercent <- max(freqDF$avgPercent)+2
    pvalDF$group <- NA
    
    p <- ggplot(freqDF, aes(x=group, y=avgCount, fill=group)) +
        geom_errorbar(aes(ymax = avgCount + sdCount , ymin=avgCount - sdCount), width=.25) +
        geom_bar(stat="identity", colour="black") + 
        facet_grid(.~inTAD) +
        scale_fill_manual(values=COL, guide_legend(title = "")) +
        theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="none") + labs(y="Gene pairs with shared enhancer", x="") + 
        geom_text(aes(label=avgCount, y= avgCount), vjust=1.25, size=5) +
        geom_text(aes(label=paste0("p=",signif(p,2)), y=ypos, x=1.5), size=5, data=pvalDF)
    ggsave(paste0(outPrefix, ".", subName, ".eh_by_group_and_inTAD.barplot.pdf"), w=3.5, h=7)
    
    p <- ggplot(freqDF, aes(x=group, y=avgPercent, fill=group)) +
        geom_errorbar(aes(ymax = avgPercent + sdPercent , ymin=avgPercent - sdPercent), width=.25) +
        geom_bar(stat="identity", colour="black") + 
        facet_grid(.~inTAD) +
        scale_fill_manual(values=COL, guide_legend(title = "")) +
        theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="none") + labs(y="Percent of pairs with shared enhancer  [%]", x="") + 
        geom_text(aes(label=signif(avgPercent,3), y= avgPercent), vjust=1.25, size=5) +
        geom_text(aes(label=paste0("p=",signif(p,2)), y=yposPercent, x=1.5), size=5, data=pvalDF)
    ggsave(paste0(outPrefix,".", subName, ".ehPercent_by_group_and_inTAD.barplot.pdf"), w=3.5, h=7)

}

#=======================================================================
# Compartments and Sub-compartments
#=======================================================================

#-----------------------------------------------------------------------
# Percent of pairs with same compartment and same sub-compartment
#-----------------------------------------------------------------------
subDF <- subset(allDF, sSamp & sTAD & sExp & sSpec)

subDF$comp_combination[subDF$same_comp_region] <- NA
subDF$subcomp_combination[subDF$same_subcomp_region] <- NA

# mark common cmopartment as NA if combination is NA
subDF$common_comp[is.na(subDF$comp_combination)] <- NA
subDF$common_subcomp[is.na(subDF$subcomp_combination)] <- NA

dc <- ddply(subDF, .(group, replicate), summarize, 
    nComp=sum(!is.na(common_comp)), 
    countComp=sum(common_comp, na.rm=TRUE), 
    percentComp=sum(common_comp, na.rm=TRUE)/sum(!is.na(common_comp))*100,
    nSubComp=sum(!is.na(common_subcomp)), 
    countSubComp=sum(common_subcomp, na.rm=TRUE), 
    percentSubComp=sum(common_subcomp, na.rm=TRUE)/sum(!is.na(common_subcomp))*100
)

freqDF <- ddply(dc, .(group), summarize,
    meanComp=mean(countComp),
    sdComp=sd(countComp),
    meanPercentComp=mean(percentComp),
    sdPercentComp=sd(percentComp),
    meanSubComp=mean(countSubComp),
    sdSubComp=sd(countSubComp),
    meanPercentSubComp=mean(percentSubComp),
    sdPercentSubComp=sd(percentSubComp)
)
    

allCompP <- fisher.test(subDF$group, subDF$common_comp)$p.value
allSubCompP <- fisher.test(subDF$group, subDF$common_subcomp)$p.value


p <- ggplot(freqDF, aes(x=group, y=meanPercentComp, fill=group)) +
    geom_errorbar(aes(ymax = meanPercentComp + sdPercentComp , ymin=meanPercentComp - sdPercentComp), width=.25) +
    geom_bar(stat="identity", color="black") + 
    geom_text(aes(label=signif(meanPercentComp, 3)), vjust=1.5, size=7) + 
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position = "none") + 
    scale_fill_manual(values=COL) + ylab("Same A/B-compartment [%]") + xlab("") + ggtitle(paste0("p=", signif(allCompP,3)))    
ggsave(paste0(outPrefix,".common_compartment.noLables.pdf"), p, w=3.5)
ggsave(paste0(outPrefix,".common_compartment.pdf"), p, w=3.5)

p <- ggplot(freqDF, aes(x=group, y=meanPercentSubComp, fill=group)) +
    geom_errorbar(aes(ymax = meanPercentSubComp + sdPercentSubComp , ymin=meanPercentSubComp - sdPercentSubComp), width=.25) +
    geom_bar(stat="identity", color="black") + 
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position = "none") + 
    scale_fill_manual(values=COL) + ylab("Same sub-compartment [%]") + xlab("") + ggtitle(paste0("p=", signif(allSubCompP <- fisher.test(subDF$group, subDF$common_subcomp)$p.value
,3))) +
    geom_text(aes(label=signif(meanPercentSubComp, 3)), vjust=1.5, size=7)
ggsave(paste0(outPrefix,".common_subcompartment.noLables.pdf"), p, w=3.5)
ggsave(paste0(outPrefix,".common_subcompartment.pdf"), p, w=3.5)

#=======================================================================
# Housekeeping genes
#=======================================================================
subDF <- subset(allDF, sSamp & sTAD & sExp & sSpec)

dc <- ddply(subDF, .(group, housekeeping), summarize, count=length(housekeeping))
summedCounts <- rep(ddply(dc, .(group), summarize, s=sum(count))$s, each=length(levels(dc$housekeeping)))
subDFcomb <- cbind(dc, cbind(summedCounts, percent=dc$count * 100 / summedCounts))

p <- ggplot(subDFcomb, aes(x=housekeeping, fill = group, y = percent)) +
    geom_bar(stat="identity", position="dodge", color="black") + theme_bw() + 
    theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position = "bottom") + 
    scale_fill_manual("", values=COL) + ylab("Gene pairs [%]") + xlab("Housekeeping genes\nin pairs") +
    geom_text(aes(label = signif(percent,3)), position=position_dodge(width=1), vjust=1.25, size=4)
ggsave(paste0(outPrefix,".housekeeping_in_pair_by_group.barplot.pdf"), p, w=3.5)

#=======================================================================
# save workspace image
#=======================================================================
save.image(WORKIMAGE_FILE)

