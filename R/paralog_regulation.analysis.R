########################################################################
#
# A script to run parts of the paralog regulation analyse 
#
########################################################################

require(stringr)        # for some string functionality
require(RColorBrewer)   # for nice colors
require(colorspace)     # for some more colors
require(GenomicRanges)  # for genomic intervals and overlap calculation
require(rtracklayer)    # for import.bed
require(plyr)           # count() function
require(data.table)     # for data.table object
require(ggplot2)        # for nice plots
require(scales)         # for proper logarithmic scales in ggplot
require(BiocParallel)   # for parallel computing

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

#-----------------------------------------------------------------------
# load some custom functions
#-----------------------------------------------------------------------
source("R/functions.plot.R")
source("R/functions.genePairs.R")

#=======================================================================
# 1.) Load data from exported WORKIMAGES
#=======================================================================
message("INFO: Start loading WORKIMAGE_FILE ...")
load(paste0(WORKIMAGE_FILE, ".sampling_and_annotation.Rdata"))
message("INFO: Finished loading WORKIMAGE_FILE.")

#-----------------------------------------------------------------------
# separate close and distal pairs and create some combined data sets
#-----------------------------------------------------------------------
sampClosePairs <- lapply(sampCisPairs, function(gP) gP[abs(gP$dist) <= MAX_DIST,])
sampClosePairsCombined <- do.call("rbind", sampClosePairs)

sampDistalPairs <- lapply(sampCisPairs, function(gP) gP[abs(gP$dist) > MAX_DIST,])
sampDistalPairsCombined <- do.call("rbind", sampDistalPairs)

sampDistEhPairsClose <- lapply(sampDistEhPairs, function(gP) gP[abs(gP$dist) <= MAX_DIST,])
sampDistEhPairsCloseCombined <- do.call("rbind", sampDistEhPairsClose)

#=======================================================================
# 2.) Run basic analysis
#=======================================================================

#-----------------------------------------------------------------------
# Validate sampling by analysing distribution of distance, linked enhancer...
#-----------------------------------------------------------------------
pdf(paste0(outPrefix, ".sampling.Dist.pdf"), w=9, h=2.25)

    par(cex=1, lwd=1.5, mfrow=c(1,4))

    paraDist <- abs(allCisPairs[,"dist"] / 10^3)
    paraDist <- paraDist[paraDist <=1000]
    sampDist <- abs(sampCisPairsCombined[,"dist"] / 10^3)
    sampDist <- sampDist[sampDist <=1000]

    hist(paraDist, 20, col=COL[1],
    main="Paralogs", xlab="Distance [kb]")
    hist(sampDist, 20, col=COL[2],
    main="Sampled by distance", xlab="Distance [kb]")   
    qqplot(paraDist, sampDist, xlab="Paralogs distance [kb]", ylab="Sampled distance [kb]", main="QQ-Plot")
    abline(0,1, col="red")    
    qqplot(paraDist, sampDist, log="xy", xlab="Paralogs distance [kb]", ylab="Sampled distance [kb]", main="QQ-Plot\n(log scale)")
    abline(0,1, col="red")    

dev.off()

pdf(paste0(outPrefix, ".sampling.DistEh.pdf"), w=9, h=4.5)

    par(cex=1, lwd=1.5, mfrow=c(2,4))

    paraDist <- abs(allCisPairs[,"dist"] / 10^3)
    paraDist <- paraDist[paraDist <=1000]
    sampDist <- abs(sampDistEhPairsCombined[,"dist"] / 10^3)
    sampDist <- sampDist[sampDist <=1000]
    paraLinkedEnhancer = tssGR[c(allCisPairs[,1], allCisPairs[,2])]$linked_enhancer 
    sampEh = tssGR$linked_enhancer[c(sampDistEhPairsCombined[,1], sampDistEhPairsCombined[,2])]


    hist(paraDist, 20, col=COL[1],
    main="Paralogs", xlab="Distance [kb]")
    hist(sampDist, 20, col=COL[2],
    main="Sampled by distance\n and enhancer", xlab="Distance [kb]")   
    qqplot(paraDist, sampDist, xlab="Paralogs distance [kb]", ylab="Sampled distance [kb]", main="QQ-Plot")
    abline(0,1, col="red")    
    qqplot(paraDist, sampDist, log="xy", xlab="Paralogs distance [kb]", ylab="Sampled distance [kb]", main="QQ-Plot\n(log scale)")
    abline(0,1, col="red")    
    
    
    # number of linked enhancers
    hist(paraLinkedEnhancer[paraLinkedEnhancer<=50], 50, col=COL[1],
    main="Paralogs", xlab="Enhancers")
    hist(sampEh[sampEh<=50], 50, col=COL[2],
    main="Sampled by distance\n and enhancer", xlab="Enhancers")    
    
    qqplot(paraLinkedEnhancer, sampEh, xlab="Enhancers in paralog genes", ylab="Enhancers in sampled genes",  main="QQ-Plot")
    abline(0,1, col="red")
    qqplot(paraLinkedEnhancer, sampEh, xlab="Enhancers in paralog genes", log="xy", ylab="Enhancers in sampled genes",  main="QQ-Plot\n(log scale)")
    abline(0,1, col="red")
dev.off()


pdf(paste0(outPrefix, ".sampling.DistEhStrand.pdf"), w=9, h=6.75)

    par(cex=1, lwd=1.5, mfrow=c(3,4))
    
    paraDist <- abs(allCisPairs[,"dist"] / 10^3)
    paraDist <- paraDist[paraDist <=1000]
    sampDist <- abs(sampDistEhStrandPairsCombined[,"dist"] / 10^3)
    sampDist <- sampDist[sampDist <=1000]
    
    paraLinkedEnhancer = tssGR[allCisPairs[,1]]$linked_enhancer + tssGR[allCisPairs[,2]]$linked_enhancer    
    sampEh = tssGR[sampDistEhStrandPairsCombined[,1]]$linked_enhancer + tssGR[sampDistEhStrandPairsCombined[,2]]$linked_enhancer
    
    strandPara <- table(allCisPairs$sameStrand) / nrow(allCisPairs) * 100
    strandSamp <- table(sampDistEhStrandPairsCombined$sameStrand) / nrow(sampDistEhStrandPairsCombined) * 100

    # distance
    hist(paraDist, 20, col=COL[1],
    main="Paralogs", xlab="Distance [kb]")
    hist(sampDist, 20, col=COL[2],
    main="Sampled by distance,\n enhancer and strand", xlab="Distance [kb]")   
    qqplot(paraDist, sampDist, xlab="Paralogs distance [kb]", ylab="Sampled distance [kb]", main="QQ-Plot")
    abline(0,1, col="red")    
    qqplot(paraDist, sampDist, log="xy", xlab="Paralogs distance [kb]", ylab="Sampled distance [kb]", main="QQ-Plot\n(log scale)")
    abline(0,1, col="red")    
    
    # number of linked enhancers
    hist(paraLinkedEnhancer[paraLinkedEnhancer<=50], 50, col=COL[1],
    main="Paralogs", xlab="Enhancers", ylab="%")
    hist(sampEh[sampEh<=50], 50, col=COL[2],
    main="Sampled by distance,\n enhancer and strand", xlab="Enhancers", ylab="%")    
    
    qqplot(paraLinkedEnhancer, sampEh, xlab="Enhancers in paralog genes", ylab="Enhancers in sampled genes",  main="QQ-Plot")
    abline(0,1, col="red")
    qqplot(paraLinkedEnhancer, sampEh, xlab="Enhancers in paralog genes", log="xy", ylab="Enhancers in sampled genes",  main="QQ-Plot\n(log scale)")
    abline(0,1, col="red")    
    
    # satrand
    bp <- barplot(strandPara, col=COL[1], main="Paralogs", ylab="%", names.arg=c("Opposite\nstrand", "Same\nstrand"))
    text(x=bp, y=strandPara, label=signif(strandPara,3), pos=1)
    bp <- barplot(strandSamp, col=COL[2], main="Sampled by distance,\n enhancer and strand", ylab="%", names.arg=c("Opposite\nstrand", "Same\nstrand"))
    text(x=bp, y=strandSamp, label=signif(strandSamp,3), pos=1)

dev.off()


pdf(paste0(outPrefix, ".sampling.DistEhLen.pdf"), h=6.75, w=9)

    par(cex=1, lwd=1.5, mfrow=c(3,4))
    
    paraDist <- abs(allCisPairs[,"dist"] / 10^3)
    paraDist <- paraDist[paraDist <=1000]
    sampDist <- abs(sampDistEhLenPairsCombined[,"dist"] / 10^3)
    sampDist <- sampDist[sampDist <=1000]
    
    paraLinkedEnhancer = tssGR[allCisPairs[,1]]$linked_enhancer + tssGR[allCisPairs[,2]]$linked_enhancer    
    sampEh = tssGR[sampDistEhLenPairsCombined[,1]]$linked_enhancer + tssGR[sampDistEhLenPairsCombined[,2]]$linked_enhancer
    
    
    paraLen = c(tssGR[allCisPairs[,1]]$gene_size, tssGR[allCisPairs[,2]]$gene_size)
    sampLen <- c(tssGR[sampDistEhLenPairsCombined[,1]]$gene_size, tssGR[sampDistEhLenPairsCombined[,2]]$gene_size)
    
    # Distance
    # distance
    hist(paraDist, 20, col=COL[1],
    main="Paralogs", xlab="Distance [kb]")
    hist(sampDist, 20, col=COL[2],
    main="Sampled by distance,\n enhancer and gene length", xlab="Distance [kb]")   
    qqplot(paraDist, sampDist, xlab="Paralogs distance [kb]", ylab="Sampled distance [kb]", main="QQ-Plot")
    abline(0,1, col="red")    
    qqplot(paraDist, sampDist, log="xy", xlab="Paralogs distance [kb]", ylab="Sampled distance [kb]", main="QQ-Plot\n(log scale)")
    abline(0,1, col="red")    
    
    # number of linked enhancers
    hist(paraLinkedEnhancer[paraLinkedEnhancer<=50], 50, col=COL[1],
    main="Paralogs", xlab="Enhancers", ylab="%")
    hist(sampEh[sampEh<=50], 50, col=COL[2],
    main="Sampled by distance,\n enhancer and gene length", xlab="Enhancers", ylab="%")    

    qqplot(paraLinkedEnhancer, sampEh, xlab="Enhancers in paralog genes", ylab="Enhancers in sampled genes",  main="QQ-Plot")
    abline(0,1, col="red")
    qqplot(paraLinkedEnhancer, sampEh, xlab="Enhancers in paralog genes", log="xy", ylab="Enhancers in sampled genes",  main="QQ-Plot\n(log scale)")
    abline(0,1, col="red")    
    
    # gene length log
    hist(log10(paraLen), 20, col=COL[1],
    main="Gene length paralogs", xlab="Gene length [log_10 kb]")
    hist(log10(sampLen), 20, col=COL[2],
    main="Gene length sampled", xlab="Gene length [log_10 kb]")
    
    qqplot(paraLen, sampLen, xlab="Gene length in paralogs", ylab="Gene length in sampled genes",  main="QQ-Plot")
    abline(0,1, col="red")
    qqplot(log10(paraLen), log10(sampLen), xlab="Gene length in paralogs", log="xy", ylab="Gene length in sampled genes",  main="QQ-Plot\n(log scale)")
    abline(0,1, col="red")
dev.off()

#-----------------------------------------------------------------------
# Analyse how many sampled pairs are actually real paralog pairs
#-----------------------------------------------------------------------

paraPairList <- list(cisUniq=allCisPairs, all=paralogPairs)
sampledPairs <- list(
    "rand"=randPairsCombined,
    "cisByDist"=sampCisPairsCombined,
    "cisByDist_Eh"=sampDistEhPairsCombined,
    "cisByDist_Eh_Strand"=sampDistEhStrandPairsCombined,
    "cisByDist_Eh_GeneLength"=sampDistEhLenPairsCombined
    )
    
# calculate percent of paralogs in sampled gene pairs
percentPara <- cbind(sapply(names(paraPairList), function(L){
    para <- paraPairList[[L]]    
    sapply(sampledPairs, percentIncluded, para)
}))

pdf(paste0(outPrefix, ".sampling.percentParalogs.pdf"))
    bp <- my.barplot(t(percentPara), beside=TRUE, ylab="Paralogs in sampled pairs [%]", addValues=TRUE, names=rownames(percentPara), srt=20, col=brewer.pal(3, "Paired")[1:2])
    legend("topleft", c("Unique paralog pairs on same chrom", "All paralog pairs"), fill=brewer.pal(3, "Paired")[1:2])
    text(bp, t(percentPara), labels=signif(t(percentPara), 3), pos=3)
dev.off()

#-----------------------------------------------------------------------
# plot number of genes in each group
#-----------------------------------------------------------------------
nrGenes = sapply(list(paralogs, nonParalogs), length)
pdf(paste0(outPrefix, ".number_of_genes.pdf"), width=3.5)

    par(cex=1.5, mgp=c(3,1,0), mar = c(7, 4, 4, 2) + 0.1)
    barP = barplot(nrGenes, beside=TRUE, col=COL2, names.arg=NA, ylab="Number of genes", ylim=c(0, 1.25*max(nrGenes))) 
    text(barP, nrGenes, nrGenes, pos=3)
    axis(1, at=barP, labels = FALSE)
    labels = c("Paralog genes", "Non-paralog genes")
    text(barP, par("usr")[3] - 0.1*max(nrGenes), srt = 45, adj = 1, labels = labels, xpd = TRUE)

dev.off()

#-----------------------------------------------------------------------
# Plot distribution of paralog group size
#-----------------------------------------------------------------------
# count paralog partners for each gene
countGenes = table(paralogPairs[,1])
countSizes = table(as.data.frame(countGenes)[,2])
countSizes = c("0"=length(nonParalogs), countSizes)
groupSizes = as.numeric(names(countSizes))

pdf(paste0(outPrefix, ".paralog_group_size.pdf"))
    par(cex=1.5, lwd=3)
    plot(groupSizes, as.vector(countSizes), xlim=c(0,40), type="o", pch=19, xlab="Number of paralogs copies", ylab="Number of genes", main="Paralog group size distribution", col=COL[1])
dev.off()

#-----------------------------------------------------------------------
# Gene size distribution
#-----------------------------------------------------------------------
paralogGeneSize = mcols(tssGR[paralogs])$gene_size / 10^3
nonParalogGeneSize = mcols(tssGR[nonParalogs])$gene_size / 10^3

# Wilcoxon-rank-sum test
ws.test = wilcox.test(paralogGeneSize, nonParalogGeneSize)

pdf(paste0(outPrefix, ".Gene_length_boxplot.pdf"), width=3.5)
    par(cex=1.5, mgp=c(3,1,0), lwd=2)

    xlabs=paste(c("Paralogs", "Non-paralogs"), "\n n = ", c(length(paralogGeneSize), length(nonParalogGeneSize)))
    my.boxplot(list(paralogGeneSize, nonParalogGeneSize), offset=0.122,
        names=xlabs, log="y", border=COL2, ylab="Gene length [kb]", 
        main=paste0("p = ", signif(ws.test$p.value, 2)))

dev.off()

#-----------------------------------------------------------------------
# Housekeeping genes
#-----------------------------------------------------------------------

genesDF <- data.frame(
    ENSG = names(tssGR),
    paralog = factor(names(tssGR) %in% paralogs, c(TRUE, FALSE), c("paralog", "non-paralog")),
    housekeeping = factor(names(tssGR) %in% housekeepingGenes, c(TRUE, FALSE), c("housekeeping", "tissue-specific"))
)

d <- ddply(genesDF, .(paralog, housekeeping), summarize, count=length(paralog))

p <- ggplot(d, aes(x=housekeeping, y=count, fill=paralog)) + 
    geom_bar(stat="identity", position="dodge") + 
    facet_grid(.~paralog) + xlab("") +
    scale_fill_manual(values=COL2, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + 
    geom_text(aes(label=count), vjust=1.25, size=5)
ggsave(p, file=paste0(outPrefix, ".genes_by_pralog_and_housekeeping.pdf"), width=3.5)
    
#-----------------------------------------------------------------------
# Number of genes
#-----------------------------------------------------------------------
nrGenes = sapply(list(paralogs, nonParalogs), length)
pdf(paste0(outPrefix, ".number_of_genes.pdf"), width=3.5)

    par(cex=1.5, mgp=c(3,1,0), mar = c(7, 4, 4, 2) + 0.1)
    barP = barplot(nrGenes, beside=TRUE, col=COL2, names.arg=NA, ylab="Number of genes", ylim=c(0, 1.25*max(nrGenes))) 
    text(barP, nrGenes, nrGenes, pos=3)
    axis(1, at=barP, labels = FALSE)
    labels = c("Paralog genes", "Non-paralog genes")
    text(barP, par("usr")[3] - 0.1*max(nrGenes), srt = 45, adj = 1, labels = labels, xpd = TRUE)

dev.off()


#-----------------------------------------------------------------------
# fraction of paralog pairs on same chromosom
#-----------------------------------------------------------------------

# Plot the percent of paralog pairs on the same chromosome
pdf(paste0(outPrefix, ".paralogPairs_on_same_chrom.barplot.pdf"), width=3.5)

        paraSameChrom <- !is.na(paralogPairsUniq$dist)
        sampSameChrom <- lapply(randPairs, function(gP) !is.na(gP$dist))

        paraPercent <- percentTrue(paraSameChrom)
        randPercent <- sapply(sampSameChrom, percentTrue) 
        
        height=c(paraPercent, mean(randPercent))
        
        # create contingency table and run Fisher test
        contab = rbind(
            para=table(paraSameChrom),  
            rand=table(unlist(sampSameChrom))
        )
        pval <- fisher.test(contab)$p.value
        
        yMax = 1.3*max(height) + sd(randPercent)
        par(cex=1.5, lwd=2)
        
        bp = my.barplot(height, 
            names=paste0(c("Paralogs\n (n=", "Random pairs\n (n="), c(nrow(paralogPairsUniq), paste0(N_RAND, "x",nrow(paralogPairsUniq))), c(")", ")")), 
            addValues=TRUE, yMax=yMax, col=COL_RAND,
            main="Gene pairs on\n same chromosome", 
            ylab="Pairs on same chromosome [%]")
        error.bar(bp,height, c(NA, sd(randPercent)))

        # p-value matrix for only two group comparison
        pval_mat = matrix(c(NA, signif(pval, 3)), 1)
        add_pval_all_allPairs(bp, ymax=1.1*max(height)+ sd(randPercent), pval_mat, min_pval=10^-16)
        
dev.off()

#---------------------------------------------------------------
# fraction of paralog pairs with same strand on same chrom
#---------------------------------------------------------------
paraPercent = sapply(list(paralogPairsUniq, distalPairs, closePairs), function(gP) percentTrue(gP[,"sameStrand"]))
randPercent = lapply(list(
    randPairs, 
    lapply(sampCisPairs, function(gP) gP[abs(gP$dist) > 10^6,]), 
    lapply(sampCisPairs, function(gP) gP[abs(gP$dist) <= 10^6,])
    ), function(gpl) {
        sapply(gpl, function(gP) percentTrue(gP[,"sameStrand"]))
        })

height=rbind(paraPercent=paraPercent, randPercent=sapply(randPercent, mean))
pdf(paste0(outPrefix, ".paralogPairs_sameStrand.barplot.pdf"))
    
    par(cex=1.5, lwd=2)
    bp = my.barplot(height, beside=TRUE,
        names=c("All pairs", "Distal cis pairs", "Close cis pairs"), 
        addValues=TRUE, col=c(COL_RAND, COL, COL),
        main="Gene pairs on same strand", 
        ylab="%")
    
    error.bar(bp, height, rbind(rep(NA, 3), sapply(randPercent, sd)))
    legend("topleft", c("Paralog pairs", "Sampled controls"), fill=COL)
dev.off()

# plot distance according to same strand as boxplot
wstest = wilcox.test(abs(closePairs$dist[closePairs$sameStrand]), abs(closePairs$dist[!closePairs$sameStrand]))
pdf(paste0(outPrefix, ".paralogPairs_sameStrand_dist.boxplot.pdf"), width=3.5)
    par(cex=1.5, lwd=2)
    bp = my.boxplot(list(abs(closePairs$dist[closePairs$sameStrand])/10^3 ,abs(closePairs$dist[!closePairs$sameStrand])/10^3),
        returnBp=TRUE, border=COL[1], names=c("Same strand", "Opposite strand"),
        ylab="Distance [kb]")
    
    add_pval_all_allPairs(1:length(bp$n), 1.1*max(bp$out), pval_mat=cbind(NA, signif(wstest$p.value, 3)), xpd=TRUE)

dev.off()

#-----------------------------------------------------------------------
# 3) linear distance between paralogs
#-----------------------------------------------------------------------

# Show bias of equaly sampled genes to paralog genes

# get paralog distances <= MAX_DIST
paraDist = abs(closePairs[,"dist"])/10^3
# get cispairs and distance form uniformaly random genes
randDist <- abs(randPairsInCisCombined[abs(randPairsInCisCombined[,"dist"])<= MAX_DIST, "dist"])/10^3
# distance form sampled pairs
sampledDist = abs(sampCisPairsCombined[abs(sampCisPairsCombined[,"dist"]) <= MAX_DIST,"dist"])/10^3

randDistPval <- wilcox.test(paraDist, randDist)$p.value
sampleDistPval <- wilcox.test(paraDist, sampledDist)$p.value
    
pdf(paste0(outPrefix, ".samped_random_para_dist.hist.pdf"))
    par(lwd=2, mfrow=c(3,1))
    par(cex=1.5,  mar=c(3, 4.1, 1.5, 2.1))
    
    hp <- hist(paraDist, 50, col=COL_RAND[1],
    main="Paralog gene pairs", xlab="")
    hr <- hist(randDist, 50, col=COL_RAND[2],
    main=paste0("Random gene pairs (p=", signif(randDistPval, 2), ")"), xlab="")
    hs <- hist(sampledDist, 50, col=COL[2],
    main=paste0("Sampled gene pairs (p=", signif(sampleDistPval, 2), ")"), xlab="")
    mtext("Distance (kb)", side=1, line=2, cex=1.5)    

dev.off()


#-----------------------------------------------------------------------
# Plot size distribution of domains
#-----------------------------------------------------------------------
sizeList = lapply(c(allTADs, speciesTADs), function(gr) width(gr)/10^3)
numberList = lapply(c(allTADs, speciesTADs), length)
sizeListMean = sapply(sizeList, mean)
sizeListMedian = sapply(sizeList, median)

pdf(paste0(outPrefix, ".Compare_domain_size.boxplot.pdf"))
    par(cex=1.3, lwd=1.5)
    my.boxplot(sizeList, names=paste0(gsub("_", " ", names(sizeList)), " n=", numberList),
        main="Size distribution of TADs",
        ylab="TAD size [kb]", col=COL_DOMAIN)
dev.off()

pdf(paste0(outPrefix, ".Compare_domain_size.log10.boxplot.pdf"))
    par(cex=1.3, lwd=1.5)
    my.boxplot(lapply(sizeList, function(x){log10(x*10^3)}), names=paste0(gsub("_", " ", names(sizeList)), " (n=", numberList,")"),
        main="Size distribution of TADs",
        ylab="TAD size [log_10 bp]", col=COL_DOMAIN)
dev.off()


#-----------------------------------------------------------------------
# Co-occurance in same TAD
#-----------------------------------------------------------------------

# Plot fraction of paralog pairs within TAD
pvalues = c()
heightMat = c()
sdMat = c()

# repeat analysis for all TADs from Dixon et al and Rao et al    
for( LABEL in names(allTADs) ){
    
    # get boolean vectors
    paraInTAD = closePairs[, LABEL]
    sampInTAD = lapply(sampClosePairs, function(gP) gP[, LABEL])

    # create contingency table and run Fisher test
    contab = rbind(
        para=table(paraInTAD),  
        rand=table(unlist(sampInTAD))
    )
    fs.test = fisher.test(contab)
    
    # get percent values for barplot
    paraPercent = percentTrue(paraInTAD)
    sampPercent = sapply(sampInTAD, percentTrue) 
    heights = c(paraPercent, mean(sampPercent, na.rm=TRUE))
    sdSamp = sd(sampPercent, na.rm=TRUE)

    pvalues = c(pvalues, fs.test$p.value)
    heightMat = cbind(heightMat, heights)
    sdMat = c(sdMat, sdSamp)

    pdf(paste0(outPrefix, ".close.paralogPairs.within_", LABEL, ".barplot.pdf"), width=3.5)

        par(cex=1.5, lwd=2)
        yMax = 1.3*max(heights) + sdSamp
        
        bp = my.barplot(heights, addValues=TRUE, yMax=yMax,
            names=c("Paralogs", "Sampled pairs"),
            main=paste("Close pairs in\n", LABEL, "TAD"), ylab="Gene pairs with TSSs in same TAD [%]", col=COL)
        error.bar(bp, heights, c(NA,  sdSamp), lwd=2)
        
        # pvalue matrix for only two group comparision
        pval_mat = matrix(c(NA, signif(fs.test$p.value, 3)), 1)
        add_pval_all_allPairs(bp, ymax=1.1*max(heights)+sdSamp, pval_mat)

    dev.off()
    
}

# Plot for TADs from all cell-types together
pdf(paste0(outPrefix, ".close.paralogPairs.within_allTADs.barplot.pdf"), width=14)

    par(cex=1.5, lwd=2)
    bp = my.barplot(heightMat, addValues=TRUE, beside=TRUE, names=gsub("_", " ", names(allTADs)), col=COL, yMax=1.3*max(heightMat), digits=2, ylab="Gene pairs in same TAD [%]")
    error.bar(bp, heightMat, rbind(rep(NA, length(allTADs)), sdMat), length=0.05)

    add_pval_two_pairs(bp, heightMat, pvalues, offset=.15*max(heightMat)) #, digits=1

    legend("topleft", c("Paralogs", "Sampled pairs"), fill=COL, bty="n")
dev.off()


#-----------------------------------------------------------------------
# Number of enhancers linked to each gene in the groups
#-----------------------------------------------------------------------

# show bias of number of enhancers linked to each gene:
paraEhPercent = percentCounts(tssGR[paralogs]$linked_enhancer)
nonParaEhPercent = percentCounts(tssGR[nonParalogs]$linked_enhancer)
mergedCounts = merge(paraEhPercent, nonParaEhPercent, all=TRUE, by="count")

paraMeanEh = mean(tssGR[paralogs]$linked_enhancer)
nonParaMeanEh = mean(tssGR[nonParalogs]$linked_enhancer)

pVal <- wilcox.test(tssGR[paralogs]$linked_enhancer, tssGR[nonParalogs]$linked_enhancer)$p.value

# plot distribution of number of enhancers
pdf(paste0(outPrefix, ".number_of_enhancers.barplot.pdf"))
    # define max number of counts and add sum up all higher counts
    maxCount = 9
    countTab = mergedCounts[seq(maxCount+1),]
    countTab[maxCount+1,] = countTab[maxCount+1,] +  colSums(mergedCounts[seq(maxCount+2, nrow(mergedCounts)),], na.rm=TRUE)
    
    height=t(countTab[,seq(2, ncol(countTab))])
    names=countTab[,1]
    names[maxCount+1] = paste0(maxCount, "+")
    par(cex=1.5, lwd=1.5)
    bp <- barplot(height, names.arg=names,
        beside=TRUE, col=COL2, main="Number of enhancers\n linked to single genes",
        xlab="Number of enhancers", ylab="Genes [%]")
    text(x=mean(bp), y=.5*max(height), labels=paste0("p=", signif(pVal, 3)))
    legend("topright", paste0(legend=c("Paralog genes", "Non-Paralogs"), " (avg=", signif(c(paraMeanEh, nonParaMeanEh), 2), ")"), fill=COL2)
dev.off()

#-----------------------------------------------------------------------
# plot fraction of pairs without and with at least one shared enhancers
#-----------------------------------------------------------------------

# get enhancer-promotor distances 
tab = data.table(dist=abs(regMap$dist), tssIDX=regMap$tss)
enhancerDist = tab[,list(
    max=max(dist),
    mean=mean(dist),
    min=min(dist),
    n=length(dist)
        ),by=tssIDX]

paralogsDists = enhancerDist[names(tssGR[enhancerDist$tssIDX]) %in% paralogs,]
nonParalogsDists = enhancerDist[names(tssGR[enhancerDist$tssIDX]) %in% nonParalogs,]

# nearest enhancer
min.test = wilcox.test(paralogsDists$min, nonParalogsDists$min)

pdf(paste0(outPrefix, ".enhancer_dist_min.pdf"), 4.5, 7)
    par(cex=1.2)

    ymax = max(c(paralogsDists$min/10^3, nonParalogsDists$min/10^3))
    ymin = min(c(paralogsDists$min/10^3, nonParalogsDists$min/10^3))

    my.boxplot(list(paralogsDists$min/10^3, nonParalogsDists$min/10^3), log="y", ylim=c(ymin, 5*ymax), names=NA, border=COL2, lwd=2, ylab="Distance to nearest enhancer (kb)")

    xlabels=paste(
        c("Paralogs", "Non-paralogs"), 
        "\n n =", 
        c(nrow(paralogsDists), nrow(nonParalogsDists)), 
        "\nMean =", 
        signif(c(mean(paralogsDists$min/10^3), mean(nonParalogsDists$min/10^3)),3)
        )
    
    axis(1, at = c(1, 2), labels=xlabels, line=1, tick=FALSE)
    pval_mat = matrix(c(NA, signif(min.test$p.value, 3)), 1)
    add_pval_all_allPairs(c(1,2), ymax=2*ymax, pval_mat,offset_fac=.5)
dev.off()

#-----------------------------------------------------------------------
# Shared enhancer between paralog genes
#-----------------------------------------------------------------------

# build a data set with the frequency of occurrence numbers
mergedCounts = merge(
    percentCounts(closePairs$commonEnhancer), 
    percentCounts(rbind.fill(lapply(sampDistEhPairs, function(gP) gP[abs(gP$dist) <= MAX_DIST,]))[,"commonEnhancer"]), all=TRUE, by="count")
        
names(mergedCounts) = c("count", "paralogs", "sampled")

# plot distribution of number of enhancers
pdf(paste0(outPrefix, ".paralogPairs.common_enhancer.barplot.pdf"))

    maxCount = 4
    countTab = mergedCounts[seq(maxCount+1),]
    countTab[maxCount+1,] = countTab[maxCount+1,] +  colSums(mergedCounts[seq(maxCount+2, nrow(mergedCounts)),], na.rm=TRUE)
    height=t(countTab[,seq(2, ncol(countTab))])
    names=countTab[,1]
    names[maxCount+1] = paste0(maxCount, "+")
    par(cex=1.5)
    bp = barplot(height, names.arg=names,
        beside=TRUE, col=COL, ylim=c(0,100),
        main=paste("Shared enhancers"),
        xlab="Number of shared enhancers", ylab="%")
    text(bp, height, signif(height,2), pos=3)
    legend("topright", legend=c("Paralog pairs", "Sampled pairs"), fill=COL)
dev.off()

# create contingency table and run Fisher test
paraHasShared = closePairs$commonEnhancer >= 1
sampHasShared = lapply(sampDistEhPairsClose, function(df){df$commonEnhancer >= 1})
contab = rbind(
    para=table(paraHasShared),  
    rand=table(unlist(sampHasShared))
)
fs.test = fisher.test(contab)

# plot fraction of pairs with at least one shared enhancer
paraPercent = percentTrue(paraHasShared)
sampPercent = sapply(sampHasShared, percentTrue) 
heights = c(paraPercent, mean(sampPercent, na.rm=TRUE))
sdSamp = sd(sapply(sampHasShared, percentTrue), na.rm=TRUE)

pdf(paste0(outPrefix, ".paralogPairs.has_shared_enhancer.barplot.pdf"), 3.5, 7)

    par(cex=1.5, lwd=2)
    ymax = max(heights) + sdSamp + .2*max(heights)
    bp = my.barplot(heights, yMax=1.4*max(heights), addValues=TRUE,
        names=c("Paralog pairs", "Sampled pairs"),
        main=paste("Gene pairs with\n shared enhancers"),
        ylab="Gene pairs with shared enhancer [%]", col=COL)
    error.bar(bp,heights, c(NA,  sdSamp), lwd=2)
    
    # pvalue matrix for only two group comparision
    pval_mat = matrix(c(NA, signif(fs.test$p.value, 3)), 1)
    add_pval_all_allPairs(bp, ymax=1.1*max(heights)+sdSamp, pval_mat)

dev.off()


#-----------------------------------------------------------------------
# Expression analysis
#-----------------------------------------------------------------------

paraExpCorSummaryList=list()

# iterate over all expression data sets
for (expName in names(expDFlist)){

    expDF = expDFlist[[expName]]

    # plot distribution of correlation coefficient
    paraExpCor = allCisPairs[,paste0(expName, "_expCor")]
    sampledExpCor = sampCisPairsCombined[,paste0(expName, "_expCor")]
    
    pdf(paste0(outPrefix,".sampled_genes.Expression_correlation.", expName, ".hist.pdf"))
        par(cex=1.5, lwd=2, mfrow=c(2,1))
        breaks = seq(-1,1,length.out=51)
        xlab = paste0("Gene expression correlation over n=", ncol(expDF), " tissues [Pearson's R]")
        hist(paraExpCor, breaks, col=COL[1],
        main=paste0("Paralog genes co-expression (",expName,")"), xlab=xlab)
        hist(sampledExpCor, breaks, col=COL[2],
        main=paste0("Sampled genes co-expression (",expName,")"), xlab=xlab)    
    dev.off()
    
    # add vectors to summary plot list
    paraExpCorSummaryList = c(paraExpCorSummaryList, list(paraExpCor, sampledExpCor))

    #-------------------------------------------------------------------
    # Compare expression and shared enhancers
    #-------------------------------------------------------------------
    maxCount = 5
    paraExpCorClose <- closePairs[,paste0(expName, "_expCor")]
    commonEhGroups = getGroupsByCounts(closePairs$commonEnhancer, maxCount)
    
    # make ggplot 
    p = ggplot(closePairs, aes(x=commonEhGroups, y=paraExpCorClose)) + 
        geom_boxplot(fill=NA, aes(colour = commonEhGroups), lwd=1.5) + 
        geom_jitter(alpha=.3, size=2, aes(colour = commonEhGroups)) + 
        scale_colour_manual(values= colorRampPalette(c("gray", COL_EH))(maxCount+1), guide=FALSE) + theme_bw() + 
        theme(text = element_text(size=20)) + labs(y="Pearson correlation coefficient", x = "Number of shared enhancers")
    ggsave(p, file=paste0(outPrefix,".close.expCor_by_shared_enhancers.", expName, ".boxplot.pdf"), width=7, height=7)

    #-------------------------------------------------------------------
    # Compare expression against localization in common TAD
    #-------------------------------------------------------------------
#~     for (tadName in names(allTADs)){
    for (tadName in c("Rao_IMR90", "Dixon_IMR90", "stable_TADs")){
        inTADgroup <- factor(closePairs[, tadName], levels=c(FALSE, TRUE), labels=c("Not in same TAD", "In same TAD"))
        
        # test difference with wilcoxon test
        ws.test = wilcox.test(paraExpCorClose ~  inTADgroup)
    
        p = ggplot(closePairs, aes(x=inTADgroup, y=paraExpCorClose)) + 
            geom_boxplot(fill=NA, aes(colour = inTADgroup), lwd=1.5) + 
            geom_jitter(alpha=.3, size=2, aes(colour = inTADgroup)) + 
            scale_colour_manual(values= colorRampPalette(c("gray", COL_TAD[1]))(2), guide=FALSE) + 
            theme_bw() + theme(text = element_text(size=20), axis.text.x = element_text(angle = 45, hjust = 1)) + labs(y="Pearson correlation coefficient", x = "", title=paste0("p = ",signif(ws.test$p.value,2)))
        ggsave(p, file=paste0(outPrefix,".close.expCor_", expName, ".by_sahred_TAD_", tadName, ".boxplot.pdf"), width=3.5, height=7)
    }

    # plot for all cell types: 
    expPlotDF = data.frame(
        expCor = rep(paraExpCorClose, length(allTADs)),
        withinTAD = unlist(lapply(names(allTADs), function(tadName) factor(closePairs[, tadName], levels=c(TRUE, FALSE), labels=c("In same TAD", "Not in same TAD")))),
        TadSource = rep(factor(gsub("_", " ", names(allTADs)), gsub("_", " ", names(allTADs))) , each=length(paraExpCorClose))
    )
    p = ggplot(expPlotDF, aes(x=withinTAD, y=expCor, colour=withinTAD)) + geom_violin(aes(colour = withinTAD), lwd=1, alpha=.25) + geom_boxplot(aes(color=withinTAD), fill=NA, width=.25, lwd=1)  + 
    facet_wrap(~ TadSource) + scale_colour_manual(values= colorRampPalette(c(COL_TAD[1], "gray"))(2), guide=FALSE) + theme_bw() + theme(text = element_text(size=20), axis.text.x = element_text(angle = 45, hjust = 1)) + labs(y="Pearson correlation coefficient", x = "") +  theme(legend.position="bottom")

    # save ggplot as pdf
    ggsave(p, file=paste0(outPrefix,".close.expCor_", expName, ".by_sahred_TAD.allTADs.boxplot.pdf"), w=7, h=7)
        
}

# boxplot (violineplot) for all data sets

pdf(paste0(outPrefix,".sampled_genes.Expression_correlation.all_data.boxplot.pdf"))
    par(cex=1.3, lwd=2, mar=c(8, 4, 1, 2))
    boxplot(paraExpCorSummaryList, border=COL, ylab="Pearson correlation", names=NA, xaxt="n")
    
    # add x-axis labels
    labPos = colMeans(matrix(seq(2*nExp), nrow=2))
    axis(1, at=labPos, labels=FALSE)
    text(x=labPos, y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]),labels=names(expDFlist), srt=45, adj=1, xpd=TRUE)
    
dev.off()

# make data.frame with additional columns indicating the groups
paraExpCorSummaryDF <- data.frame(
                expCor = unlist(paraExpCorSummaryList), 
                dataset = rep(
                    rep(
                        paste0(gsub('_', ' ', names(expDFlist)), "\n n=", sapply(expDFlist, ncol), " tissues" )
                        , each=2), 
                    times = sapply(paraExpCorSummaryList,length)
                    ),
                genepairs = rep(
                    rep(c("Paralog genes", "Sampled genes"), nExp),
                    sapply(paraExpCorSummaryList,length)
                    )
                )


p = ggplot(paraExpCorSummaryDF, na.rm=TRUE, aes(x=genepairs, y = expCor, fill=genepairs)) + 
    geom_violin(adjust = .2) + 
    geom_boxplot(aes(fill=NULL), width=.2, colour="black", outlier.shape = NA) + 
    facet_grid(.~dataset) + 
    theme_bw() + labs(y="Pearson correlation coefficient", x = "", title= "Co-expression correlation over tissues") + 
    scale_x_discrete(labels=NULL)  + 
    theme(legend.position = "bottom") + guides(fill=guide_legend(title="")) +  
    theme(text = element_text(size=15)) + scale_fill_manual(values=COL) 
ggsave(paste0(outPrefix,".sampled_genes.Expression_correlation.all_data.violinplot.pdf"))


p + coord_flip() + facet_grid(dataset~.)
ggsave(paste0(outPrefix,".sampled_genes.Expression_correlation.all_data.violinplot.rotated.pdf"))

#=======================================================================
# 5.) Compare one2one orthologs of human paralogs in mouse and dog
#=======================================================================

# iterate over all species:
for (orgStr in names(orgStr2Name)){
    
    orgName = orgStr2Name[orgStr]
    
    # check who many pairs have one-to-one orthologs
    pdf(paste0(outPrefix, ".paralogPairs_One2One_orthologs.", orgName, ".barplot.pdf"), width=3.5)
        randPercent <- sapply(randPairsInCis, function(gP) percentTrue(gP[,paste0(orgStr, "_one2one")]))
        randN <- sapply(randPairsInCis, function(gP) sum(gP[,paste0(orgStr, "_one2one")]))
        height <- c(percentTrue(allCisPairs[,paste0(orgStr, "_one2one")]), mean(randPercent))
        heightN <- c(sum(allCisPairs[,paste0(orgStr, "_one2one")]), round(mean(randN)))
        
        par(lwd=2, cex=1.5)
        bp = my.barplot(height, 
            names=c("Paralog pairs\n(on same chrom)", "Random pairs\n(on same chrom)"), 
            addValues=TRUE, col=COL_RAND,
            customValues=paste0(signif(height,3),"%\n(n=",heightN, ")"), yMax=1.3*max(height),
            main=paste("One-to-one\n orthologs\n in ", orgName), 
            ylab="Pairs with one-to-one orthologs [%]")
        
        error.bar(bp,height, c(NA, sd(randPercent)))
    dev.off()
    
    # check who many of the ortholog pairs are located on the same chrom.
    pdf(paste0(outPrefix, ".paralogPairs_orthologs_onSameChrom.", orgName, ".barplot.pdf"), width=3.5)

        paraPercent = percentTrue(allCisPairs[allCisPairs[,paste0(orgStr, "_one2one")], paste0(orgStr, "_sameChrom")])
        
        randPercent = sapply(randPairsInCis, function(gP) percentTrue(gP[gP[,paste0(orgStr, "_one2one")], paste0(orgStr, "_sameChrom")]))        

        paraN = sum(allCisPairs[allCisPairs[,paste0(orgStr, "_one2one")], paste0(orgStr, "_sameChrom")])        
        randN <- sapply(randPairsInCis, function(gP) sum(gP[gP[,paste0(orgStr, "_one2one")], paste0(orgStr, "_sameChrom")]))

        height = c(paraPercent , mean(randPercent))
        heightN = c(paraN , round(mean(randN)))

        par(lwd=2, cex=1.5)
        bp = my.barplot(height, 
            names=c("Orthologs of\n paralog pairs\n", "Orthologs of\n random pairs\n"), 
            addValues=TRUE, col=c(COL_ORTHO[1], COL_RAND[2]),
            customValues=paste0(signif(height,3),"%\n(n=",heightN,")"), yMax=1.3*max(height),
            main=paste("One-to-one \n orthologs\n on same chrom\n in", orgName), 
            ylab="Pairs on same chromosomes [%]")        
        error.bar(bp,height, c(NA, sd(randPercent)))
    dev.off()
    
    # compare linear distance between orthologs of human paralogs and random genes
    orthoDistStr = paste0(orgStr, "_dist")
    orthoDist = abs(allCisPairs[abs(allCisPairs[, orthoDistStr]) <= MAX_DIST, orthoDistStr]) /10^3
    randDist = unlist(
        lapply(randPairsInCis, function(d){
            abs(d[abs(d[, orthoDistStr]) <= MAX_DIST, orthoDistStr]) / 10^3
            })
        )
    
    randDistPval <- wilcox.test(orthoDist, randDist)$p.value
    
    pdf(paste0(outPrefix, ".paralogPairs_orthologs_dist_", orgName, ".hist.pdf"))

        par(lwd=2, mfrow=c(2,1), cex=1.5, mar=c(4, 4, 3, 2))
        # ‘c(bottom, left, top, right)’ The default is ‘c(5, 4, 4, 2)
        
        hist(orthoDist, 50, col=COL_ORTHO[1],
            main="Distance between orthologs\nof paralogs", xlab="")
        hist(randDist, 50, col=COL_RAND[2],
            main=paste0("Distance between orthologs\n of random gene pairs (p=", signif(randDistPval, 3), ")"), xlab="")    
        mtext("Distance (kb)", side=1, line=2, cex=1.5)   
    dev.off()
    
    # compare distance to dist-wighted sampled pairs    

    # take only those pairs into acount that have one-to-one orthologs on same chrom within 1 MB.
    caseSubset = !is.na(closePairs[, orthoDistStr]) & abs(closePairs[, orthoDistStr]) <= MAX_DIST
    sampSubset = lapply(sampClosePairs, function(gP) !is.na(gP[, orthoDistStr]) & abs(gP[, orthoDistStr]) <= MAX_DIST)
    
    # get distances of four sets
    paraDist = abs(closePairs[caseSubset,"dist"]) /10^3
    randDist = abs(unlist(lapply(1:N_RAND, function(i) sampClosePairs[[i]][sampSubset[[i]],"dist"]))) / 10^3

    orthoDist = abs(closePairs[caseSubset, orthoDistStr]) /10^3
    randOrthoDist = abs(unlist(
            lapply(1:N_RAND, function(i) sampClosePairs[[i]][sampSubset[[i]], orthoDistStr]) 
        )) / 10^3
        
    pdf(paste0(outPrefix, ".close.orthologs_dist_sampled_", orgName, ".hist.pdf"))
        par(lwd=2, cex=1.5, mfrow=c(2,2))

        hist(paraDist, 50, col=COL[1],
            main="Distance between\n paralogs", xlab="Distance (kb)")
        hist(randDist, 50, col=COL[2],
            main="Distance between\n sampled gene pairs", xlab="Distance (kb)")    
        hist(orthoDist, 50, col=COL_ORTHO[1],
            main="Distance between\n orthologs of paralogs", xlab="Distance (kb)")
        hist(randOrthoDist, 50, col=COL_ORTHO[2],
            main="Distance between\n orthologs of sampled gene pairs", xlab="Distance (kb)")    
    dev.off()

    # test for differences
    ws.test.para = wilcox.test(paraDist, randDist)
    ws.test.ortho = wilcox.test(orthoDist, randOrthoDist)
    
    # same as boxplot
    pdf(paste0(outPrefix, ".close.orthologs_dist_sampled_", orgName, ".boxplot.pdf"))
        par(lwd=2, cex=1.5)
        my.boxplot(list(paraDist, orthoDist, randDist, randOrthoDist), 
        ylab="Distance [kb]", border=c(COL, COL_ORTHO)[c(1,3,2,4)], names=c("paralogs", "orthologs of paralogs", "sampled pairs", "orthologs of\n sampled pairs"))
    dev.off()

    # correlate linear distances of paralogs in human with mouse orthologs
    pdf(paste0(outPrefix, ".close.orthologs_dist_", orgName, ".dotplot.pdf"))

        r = cor(paraDist, orthoDist)
        #p = cor.test(paraDist, orthoDist)$p.value
        
        par(lwd=2, cex=1.5)
        plot(paraDist,  orthoDist,
            xlim=c(1, 1000), ylim=c(1, 1000), log="xy",
            main=paste("Distance in human and", orgName),
            xlab="Distance in human [kb]", ylab=paste("Distance in", orgName, "[kb]"), col=COL_ORTHO[1] ) #, col=rgb(0,0,0,.5)
        legend("topleft", paste("R =", signif(r, 3)))
         
    dev.off()

    # correlate linear distances of paralogs in human with mouse orthologs
    pdf(paste0(outPrefix, ".close.orthologs_sampled_dist_", orgName, ".dotplot.pdf"))
        
        r = cor(randDist, randOrthoDist)
        #p = cor.test(randDist, randOrthoDist)$p.value
        
        par(lwd=2, cex=1.5)
        plot(randDist,  randOrthoDist,
            xlim=c(1, 1000), ylim=c(1, 1000), log="xy",
            main=paste("Distance in human and", orgName),
            xlab="Distance in human [kb]", ylab=paste("Distance in", orgName, "[kb]"), col=COL_ORTHO[2] ) #, col=rgb(0,0,0,.5)
        legend("topleft", paste("R =", signif(r, 3)))
         
    dev.off()
    
    #-------------------------------------------------------------------
    # co-occurances in the same TAD in other organism 
    # get boolean vectors
    orthoInTAD = closePairs[closePairs[,paste0(orgStr, "_one2one")], paste0(orgStr, "_TAD")]
    sampInTAD = lapply(sampClosePairs, function(gP) gP[gP[,paste0(orgStr, "_one2one")], paste0(orgStr, "_TAD")])

    # create contingency table and run Fisher test
    contab = rbind(
        ortho=table(orthoInTAD),  
        rand=table(unlist(sampInTAD))
    )
    fs.test = fisher.test(contab)
        
    pdf(paste0(outPrefix, ".close.orthologs_inSameTAD", orgName, ".barplot.pdf"), width=3.5)

        paraPercent = percentTrue(orthoInTAD)
        randPercent = sapply(sampInTAD, percentTrue)        
        height = c(paraPercent , mean(randPercent))

        par(lwd=2, cex=1.5)
        bp = my.barplot(height, yMax=1.4*max(height),
            names=c("Paralog pairs\n", "Sampled pairs\n"), 
            addValues=TRUE, col=COL_ORTHO,
            main=paste("One-to-one \n orthologs\n in same TAD\n in", orgName), 
            ylab="Cis pairs in same TAD [%]")        
        error.bar(bp,height, c(NA, sd(randPercent)))

        # pvalue matrix for only two group comparision
        pval_mat = matrix(c(NA, signif(fs.test$p.value, 3)), 1)
        add_pval_all_allPairs(bp, ymax=1.1*max(height)+sd(randPercent), pval_mat)
        
    dev.off()
    #-------------------------------------------------------------------
    # co-occurances in the same TAD in other organism as fraction of only those that are also in same TAD in human 
    
    # iterate over all TAD types
    for (tadName in names(allTADs)){

        paraInTAD <- closePairs[, tadName]
        sampInTAD = lapply(sampClosePairs, function(gP) gP[, tadName])

        humanAndOrthoInTAD = closePairs[closePairs[,paste0(orgStr, "_one2one")] & paraInTAD, paste0(orgStr, "_TAD")]
        humanAndRandInTAD = lapply(1:N_RAND, function(i){ 
            gP = sampClosePairs[[i]]
            gP[gP[,paste0(orgStr, "_one2one")] & sampInTAD[[i]], paste0(orgStr, "_TAD")]
            })

        # create contingency table and run Fisher test
        contab = rbind(
            para=table(humanAndOrthoInTAD),  
            rand=table(unlist(humanAndRandInTAD))
        )
        fs.test = fisher.test(contab)

        pdf(paste0(outPrefix, ".close.orthologs_inSameTAD", orgName, ".from_", tadName,".barplot.pdf"), width=3.5)

            paraPercent = percentTrue(humanAndOrthoInTAD)
            randPercent = sapply(humanAndRandInTAD, percentTrue)        
            height = c(paraPercent , mean(randPercent))
    
            par(lwd=2, cex=1.5)
            bp = my.barplot(height, yMax=1.4*max(height),
                names=c("Paralog pairs\n", "Sampled pairs\n"), 
                addValues=TRUE, col=COL_ORTHO,
                main=paste("One-to-one \n orthologs\n with conserved TAD\n in", orgName), 
                ylab="Conserved shared TAD [%]")        
            error.bar(bp,height, c(NA, sd(randPercent)))
    
            # pvalue matrix for only two group comparision
            pval_mat = matrix(c(NA, signif(fs.test$p.value, 3)), 1)
            add_pval_all_allPairs(bp, ymax=1.1*max(height)+sd(randPercent), pval_mat)
            
        dev.off()
    }
    
    #-------------------------------------------------------------------
    # Compare Hi-C counts of close orthologs
      
    plotDF = data.frame(
        group = c(rep("paralogs", nrow(closePairs)), rep("sampled", nrow(sampClosePairsCombined))),
        HiCraw = c(closePairs[,paste0(orgStr, "_HiC")], sampClosePairsCombined[,paste0(orgStr, "_HiC")]),
        HiCnorm = c(closePairs[,paste0(orgStr, "_HiCnorm")], sampClosePairsCombined[,paste0(orgStr, "_HiCnorm")])
    )
    
    #------------------------------------------------------------------------
    xlabels=paste0(c("Paralogs", "Sampled"), 
        "\n n=", applyToSubset(plotDF, function(v) sum(!is.na(v)), "HiCraw", "group"), 
        "\n med=", signif(applyToSubset(plotDF, median, "HiCraw", "group", na.rm=TRUE), 3),
        "\n avg=", signif(applyToSubset(plotDF, mean, "HiCraw", "group", na.rm=TRUE), 3)
        )
    ws.test = wilcox.test(HiCnorm ~ group, data=plotDF)
    p = ggplot(plotDF, aes(x=group, y=HiCraw, colour=group)) + 
        #geom_jitter(alpha=.3, size=2, aes(colour = group), guide=FALSE) + 
        geom_boxplot(lwd=1.5) + scale_y_log10() + annotation_logticks(sides="l") +
        scale_color_manual(values=COL_ORTHO, guide=FALSE) +
        theme_bw() + theme(text = element_text(size=20)) + 
        labs(y=paste("Hi-C counts in", orgName), x="", title=paste0("p = ", signif(ws.test$p.value, 3))) + scale_x_discrete(labels=xlabels)
    ggsave(p, file=paste0(outPrefix, ".paralogPairs_orthologs_HiCraw", orgName, ".boxplot.pdf"), width=3.5, height=7)
    
    # Normalized Hi-C counts
    #------------------------------------------------------------------------
    xlabels=paste0(c("Paralogs", "Sampled"), 
        "\n n=", applyToSubset(plotDF, function(v) sum(!is.na(v)), "HiCnorm", "group"), 
        "\n med=", signif(applyToSubset(plotDF, median, "HiCnorm", "group", na.rm=TRUE), 3),
        "\n avg=", signif(applyToSubset(plotDF, mean, "HiCnorm", "group", na.rm=TRUE), 3)
        )
    ws.test = wilcox.test(HiCnorm ~ group, data=plotDF)
    p = ggplot(plotDF, aes(x=group, y=HiCnorm, colour = group)) + 
        #geom_jitter(alpha=.3, size=2, aes(colour = group), guide=FALSE) + 
        geom_boxplot(lwd=1.5) + 
        scale_color_manual(values=COL_ORTHO, guide=FALSE) +
        theme_bw() + theme(text = element_text(size=20)) + 
        labs(y=paste("log2(observed / expected) Hi-C counts in", orgName), x="", title=paste0("p = ", signif(ws.test$p.value, 3))) + scale_x_discrete(labels=xlabels)
    ggsave(p, file=paste0(outPrefix, ".paralogPairs_orthologs_HiCnorm", orgName, ".boxplot.pdf"), width=3.5, height=7)

    #-------------------------------------------------------------------
    # do the same for distal cis pairs
    #-------------------------------------------------------------------
    plotDF = data.frame(
        group = c(rep("paralogs", nrow(distalPairs)), rep("sampled", nrow(sampDistalPairsCombined))),
        HiCraw = c(distalPairs[,paste0(orgStr, "_HiC")], sampDistalPairsCombined[,paste0(orgStr, "_HiC")]),
        HiCnorm = c(distalPairs[,paste0(orgStr, "_HiCnorm")], sampDistalPairsCombined[,paste0(orgStr, "_HiCnorm")])
    )
    plotDF$HiCrawNoZero = plotDF$HiCraw
    plotDF$HiCrawNoZero[plotDF$HiCraw == 0] = NA
    plotDF$HiCnormNoZero = plotDF$HiCnorm
    plotDF$HiCnormNoZero[plotDF$HiCnorm == 0] = NA
    
    #------------------------------------------------------------------------
    # distal raw Hi-C
    xlabels=paste0(c("Paralogs", "Sampled"), 
        "\n n=", applyToSubset(plotDF, function(v) sum(!is.na(v)), "HiCraw", "group"), 
        "\n med=", signif(applyToSubset(plotDF, median, "HiCraw", "group", na.rm=TRUE), 3),
        "\n avg=", signif(applyToSubset(plotDF, mean, "HiCraw", "group", na.rm=TRUE), 3)
        )
    ws.test = wilcox.test(HiCraw ~ group, data=plotDF)    
    p = ggplot(plotDF, aes(x=group, y=HiCraw, colour = group)) + 
        #geom_jitter(alpha=.3, size=2, aes(colour = group), guide=FALSE) + 
        geom_boxplot(lwd=1.5) + scale_y_log10() + annotation_logticks(sides="l") + 
        scale_color_manual(values=COL_ORTHO, guide=FALSE) +
        theme_bw() + theme(text = element_text(size=20)) + 
        labs(y=paste("Hi-C counts in ", orgName), x="", title=paste0("p = ", signif(ws.test$p.value, 3))) + scale_x_discrete(labels=xlabels)
    ggsave(p, file=paste0(outPrefix, ".paralogPairs_orthologs.distalPairs_HiCraw.", orgName, ".boxplot.pdf"), width=3.5, height=7)
    
    #------------------------------------------------------------------------
    # Normalized Hi-C counts
    xlabels=paste0(c("Paralogs", "Sampled"), 
        "\n n=", applyToSubset(plotDF, function(v) sum(!is.na(v)), "HiCnorm", "group"), 
        "\n med=", signif(applyToSubset(plotDF, median, "HiCnorm", "group", na.rm=TRUE), 3),
        "\n avg=", signif(applyToSubset(plotDF, mean, "HiCnorm", "group", na.rm=TRUE), 3)
        )    
    ws.test = wilcox.test(HiCnorm ~ group, data=plotDF)
    p = ggplot(plotDF, aes(x=group, y=HiCnorm, colour = group)) + 
        #geom_jitter(alpha=.3, size=2, aes(colour = group), guide=FALSE) + 
        geom_boxplot(lwd=1.5) + scale_y_log10(breaks=trans_breaks("log10", function(y) 10^y), labels=trans_format("log10", math_format(10^.x))) + 
        scale_color_manual(values=COL_ORTHO, guide=FALSE) +
        theme_bw() + theme(text = element_text(size=20)) + 
        labs(y=paste("Normalized Hi-C counts in ", orgName), x="", title=paste0("p = ", signif(ws.test$p.value, 3))) + scale_x_discrete(labels=xlabels)
    ggsave(p, file=paste0(outPrefix, ".paralogPairs_orthologs.distalPairs_HiCnorm.", orgName, ".boxplot.pdf"), width=3.5, height=7)
    
    #------------------------------------------------------------------------
    # No zero raw counts
    xlabels=paste0(c("Paralogs", "Sampled"), 
        "\n n=", applyToSubset(plotDF, function(v) sum(!is.na(v)), "HiCrawNoZero", "group"), 
        "\n med=", signif(applyToSubset(plotDF, median, "HiCrawNoZero", "group", na.rm=TRUE), 3),
        "\n avg=", signif(applyToSubset(plotDF, mean, "HiCrawNoZero", "group", na.rm=TRUE), 3)
        )
    ws.test = wilcox.test(HiCrawNoZero ~ group, data=plotDF)
    p = ggplot(plotDF, aes(x=group, y=HiCrawNoZero, colour = group)) + 
        #geom_jitter(alpha=.3, size=2, aes(colour = group), guide=FALSE) + 
        geom_boxplot(lwd=1.5) + scale_y_log10() + annotation_logticks(sides="l") + 
        scale_color_manual(values=COL_ORTHO, guide=FALSE) +
        theme_bw() + theme(text = element_text(size=20)) + 
        labs(y=paste("Hi-C counts in ", orgName), x="", title=paste0("p = ", signif(ws.test$p.value, 3))) + scale_x_discrete(labels=xlabels)
    ggsave(p, file=paste0(outPrefix, ".paralogPairs_orthologs.distalPairs_HiCrawNoZero.", orgName, ".boxplot.pdf"), width=3.5, height=7)    
    #------------------------------------------------------------------------
    # Normalized Hi-C counts without Zero counts
    xlabels=paste0(c("Paralogs", "Sampled"), 
        "\n n=", applyToSubset(plotDF, function(v) sum(!is.na(v)), "HiCnormNoZero", "group"), 
        "\n med=", signif(applyToSubset(plotDF, median, "HiCnormNoZero", "group", na.rm=TRUE), 3),
        "\n avg=", signif(applyToSubset(plotDF, mean, "HiCnormNoZero", "group", na.rm=TRUE), 3)
        )    
    ws.test = wilcox.test(HiCnormNoZero ~ group, data=plotDF)
    p = ggplot(plotDF, aes(x=group, y=HiCnormNoZero, colour=group)) + 
        #geom_jitter(alpha=.3, size=2, aes(colour = group), guide=FALSE) + 
        geom_boxplot(lwd=1.5) + scale_y_log10(breaks=trans_breaks("log10", function(y) 10^y), labels=trans_format("log10", math_format(10^.x))) + 
        scale_color_manual(values=COL_ORTHO, guide=FALSE) +
        theme_bw() + theme(text = element_text(size=20)) + 
        labs(y=paste("Normalized Hi-C counts in ", orgName), x="", title=paste0("p = ", signif(ws.test$p.value, 3))) + scale_x_discrete(labels=xlabels)
    
    # save ggplot as pdf
    ggsave(p, file=paste0(outPrefix, ".paralogPairs_orthologs.distalPairs_HiCnormNoZero.", orgName, ".boxplot.pdf"), width=3.5, height=7)
    
}


#=======================================================================
# 6.) Mouse and dog paralog analysis
#=======================================================================

speciesParDFList <- lapply(orgStr2Name, function(orgName) read.table(paste0(outPrefix, ".", orgName, ".allPairs_broad.csv"), header=TRUE, sep="\t"))
names(speciesParDFList) <- orgStr2Name


speciesParDF <- do.call("rbind", speciesParDFList)
speciesParDF$species <- factor(rep(orgStr2Name, sapply(speciesParDFList, nrow)), orgStr2Name)


#------------------------------------------------------------------------
# Analyse paralogs in TAD
#------------------------------------------------------------------------

# select only close pairs
subDF <- speciesParDF[speciesParDF$distGroup == "close",]

d <- ddply(subDF, .(species, group, replicate), summarize, count=sum(TAD), percent=percentTrue(TAD), n=length(TAD))

# combine replicates by taking average and sd of counts and percentages
freqDF <- ddply(d, .(species, group), summarize, avgCount=mean(count, na.rm=TRUE), sdCount=sd(count, na.rm=TRUE), avgPercent=mean(percent, na.rm=TRUE), sdPercent=sd(percent, na.rm=TRUE))


pvalDF <- ddply(subDF, .(species), summarize, p=fisher.test(TAD, group)$p.value)

p <- ggplot(freqDF, aes(x=group, y=avgPercent)) +
    geom_errorbar(aes(ymax = avgPercent + sdPercent , ymin=avgPercent - sdPercent), position=position_dodge(width=0.9), width=.25) +
    geom_bar(aes(fill=group), stat="identity", color="black") + 
    facet_grid(.~species) + 
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="none") + scale_fill_manual(values=COL) + 
    ylab("Gene pairs in same TAD [%]") + xlab("") +
    geom_text(aes(label=paste0(signif(avgPercent, 3), "%")), vjust=1.25, size=5) +
    geom_text(aes(label=paste0("p=",signif(p,2)), x=1.5, y=1.1*max(freqDF$avgPercent, na.rm=TRUE)), size=5, data=pvalDF)
    
ggsave(p, file=paste0(outPrefix, ".close.TAD_by_group_and_species.barplot.pdf"), width=3.5)    

#-----------------------------------------------------------------------
# Analyse HiC ontact frequencies
message("Start to analyse Hi-C contacts...")


HiCcolumnsSpecies = c(c("HiC", "HiCnorm"), paste0(c("HiC", "HiCnorm"), "NoZero"))
HiClabelsSpecies = c(
HiC="Hi-C contacts",
HiCnorm="Normalized Hi-C contacts",
HiCNoZero="Hi-C contacts",
HiCnormNoZero="Normalized Hi-C contacts"
)


for (HiCcol in HiCcolumnsSpecies ) {
    
    message(paste("INFO: Plotting:", HiCcol))
    HiClab <- HiClabelsSpecies[HiCcol]
    
    #------------------------------------------------------------------------
    # Hi-C of all pairs
    #------------------------------------------------------------------------
    subDF <- speciesParDF
    
    summaryDF <- ddply(subDF, .(species, group), summarize, n=sum(!is.na(get(HiCcol))), med=median(get(HiCcol), na.rm=TRUE), avg=mean(get(HiCcol), na.rm=TRUE))
    
    pvalDF <- ddply(subDF, .(species), summarize, p=wilcox.test(get(HiCcol) ~ group)$p.value)
    pvalDF$group="paralog"
        
    p = ggplot(subDF, aes(x=group, y=get(HiCcol))) +
        geom_boxplot(aes(color=group), lwd=1.5) + 
        facet_grid(.~species) + 
        scale_y_log10(breaks=trans_breaks("log10", function(y) 10^y), labels=trans_format("log10", math_format(10^.x))) +
        scale_color_manual(values=COL, guide_legend(title = "")) +
        theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="none") + 
        labs(y=HiClab, x="") +
        geom_text(aes(label=paste0("p=",signif(p,2)), x=1.5, y=1.1*max(subDF[,HiCcol], na.rm=TRUE)), size=5, data=pvalDF) +
#~         geom_text(aes(label=paste0("n=",n, "\nmed=", signif(med,3), "\navg=", signif(avg,3)),  y=min(subDF[,HiCcol], na.rm=TRUE)), vjust = -.1, data=summaryDF)
        geom_text(aes(label=paste0("n=",n),  y=min(subDF[,HiCcol], na.rm=TRUE)), vjust = -.1, data=summaryDF)
    
    ggsave(p, file=paste0(outPrefix, ".all_pairs.", HiCcol, ".by_group_and_species.ggboxplot.pdf"), width=3.5)

    #------------------------------------------------------------------------
    # Hi-C of close pairs
    #------------------------------------------------------------------------
    
    # select only close pairs
    subDF <- speciesParDF[speciesParDF$distGroup == "close",]
    
    summaryDF <- ddply(subDF, .(species, group), summarize, n=sum(!is.na(get(HiCcol))), med=median(get(HiCcol), na.rm=TRUE), avg=mean(get(HiCcol), na.rm=TRUE))
    
    pvalDF <- ddply(subDF, .(species), summarize, p=wilcox.test(get(HiCcol) ~ group)$p.value)
    pvalDF$group="paralog"

    p = ggplot(subDF, aes(x=group, y=get(HiCcol))) +
        geom_boxplot(aes(color=group), lwd=1.5) + 
        facet_grid(.~species) + 
        scale_y_log10(breaks=trans_breaks("log10", function(y) 10^y), labels=trans_format("log10", math_format(10^.x))) +
        scale_color_manual(values=COL, guide_legend(title = "")) +
        theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="none") + 
        labs(y=HiClab, x="") +
        geom_text(aes(label=paste0("p=",signif(p,2)), x=1.5, y=1.1*max(subDF[,HiCcol], na.rm=TRUE)), size=5, data=pvalDF) +
#~         geom_text(aes(label=paste0("n=",n, "\nmed=", signif(med,3), "\navg=", signif(avg,3)),  y=min(subDF[,HiCcol], na.rm=TRUE)), vjust = -.1, data=summaryDF)
        geom_text(aes(label=paste0("n=",n),  y=min(subDF[,HiCcol], na.rm=TRUE)), vjust = -.1, data=summaryDF)
    
    ggsave(p, file=paste0(outPrefix, ".close.", HiCcol, ".by_group_and_species.ggboxplot.pdf"), width=3.5)
    
    #------------------------------------------------------------------------
    # Hi-C of distal pairs
    #------------------------------------------------------------------------
    
    # select only close pairs
    subDF <- speciesParDF[speciesParDF$distGroup == "distal",]
    
    summaryDF <- ddply(subDF, .(species, group), summarize, n=sum(!is.na(get(HiCcol))), med=median(get(HiCcol), na.rm=TRUE), avg=mean(get(HiCcol), na.rm=TRUE))
    
    pvalDF <- ddply(subDF, .(species), summarize, p=wilcox.test(get(HiCcol) ~ group)$p.value)
    pvalDF$group="paralog"

    p = ggplot(subDF, aes(x=group, y=get(HiCcol))) +
        geom_boxplot(aes(color=group), lwd=1.5) + 
        facet_grid(.~species) + 
        scale_y_log10(breaks=trans_breaks("log10", function(y) 10^y), labels=trans_format("log10", math_format(10^.x))) +
        scale_color_manual(values=COL, guide_legend(title = "")) +
        theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="none") + 
        labs(y=HiClab, x="") +
        geom_text(aes(label=paste0("p=",signif(p,2)), x=1.5, y=1.1*max(subDF[,HiCcol], na.rm=TRUE)), size=5, data=pvalDF) +
#~         geom_text(aes(label=paste0("n=",n, "\nmed=", signif(med,3), "\navg=", signif(avg,3)),  y=min(subDF[,HiCcol], na.rm=TRUE)), vjust = -.1, data=summaryDF)
        geom_text(aes(label=paste0("n=",n),  y=min(subDF[,HiCcol], na.rm=TRUE)), vjust = .1, data=summaryDF)
    
    ggsave(p, file=paste0(outPrefix, ".distal.", HiCcol, ".by_group_and_species.ggboxplot.pdf"), width=3.5)

}


#~ # save temporary image file
#~ save.image(paste0(WORKIMAGE_FILE, ".after_species_analysis.Rdata"))
#=======================================================================
stop("INFO: Stopped script HERE!")
#=======================================================================


#=======================================================================
# Compare sahred enhancer positions
#=======================================================================

nEhPos <- 3
nParaPairs <- nrow(closePairs)
nRandPairs <- nrow(sampDistEhPairsCloseCombined)
ehColNames <- c("eh_up", "eh_cent", "eh_down")
ehPosFactor <- factor(c("upstream", "center", "downstream"), levels=c("upstream", "center", "downstream"), labels=c("upstream", "center", "downstream"))

ehPosDF = data.frame(
        "group"=c(rep("paralogs", nParaPairs*nEhPos), rep("sampled", nRandPairs*nEhPos)),
        "enahncer_position"= unlist(list(rep(ehPosFactor, each=nParaPairs), rep(ehPosFactor, each=nRandPairs))),
        "eh_counts"=c(unlist(closePairs[,ehColNames]), unlist(sampDistEhPairsCloseCombined[,ehColNames])),
        "sameIMR90_TAD"=factor(c(rep(closePairs[,"Rao_IMR90"], nEhPos), rep(sampDistEhPairsCloseCombined[,"Rao_IMR90"], nEhPos)), levels=c(TRUE, FALSE), labels=c("same TAD (Rao IMR90)", "not same TAD"))
)

pdf(paste0(outPrefix,".enhancer_position_same_strand.paralogs_sampled.1-5_enhancers.barplot.pdf"))
    ggplot(ehPosDF[!is.na(ehPosDF$eh_count) & ehPosDF$eh_count <= 6 & ehPosDF$eh_count >=1,]) + 
        geom_bar(aes(factor(eh_counts), fill=enahncer_position), width=.6, position="dodge") +
        theme_bw() + facet_grid(group~sameIMR90_TAD, scales="free_y") + 
        labs(y="Counts", x="Number of enhancers", title= "Position of shared enhancer") + scale_fill_manual(values=COL_EH_POS) 
dev.off()

pdf(paste0(outPrefix,".enhancer_position_same_strand.paralogs_sampled.number_of_enhancers.barplot.pdf"))
    ggplot(ehPosDF[!is.na(ehPosDF$eh_count) & ehPosDF$eh_count >=1,]) + 
        geom_bar(aes(enahncer_position, fill=factor(eh_counts), weight=eh_counts)) +
        theme_bw() + facet_grid(group~sameIMR90_TAD, scales="free_y", margins=TRUE) + 
        labs(y="Summed counts", x="Enhancer position", title= "Shared enhancer position relative to genes pairs") + 
        scale_fill_manual(values=colorRampPalette(brewer.pal(9, "Oranges"))(length(unique(ehPosDF$eh_count))), guide = guide_legend("Enhancer\ncounts")) + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf(paste0(outPrefix,".enhancer_position_same_strand.paralogs_sampled.all_summed.barplot.pdf"))
    ggplot(ehPosDF) +  #, 
        geom_bar(aes(enahncer_position, fill=enahncer_position, weight=eh_counts)) +
        theme_bw() + facet_grid(group~sameIMR90_TAD, scales="free_y") + 
        labs(y="Combined counts", title= "Position of shared enhancer") + scale_fill_manual(values=COL_EH_POS) 
        
dev.off()

#=======================================================================
# Compare jung vs. old pairs
#=======================================================================

nSpecies=length(orgStr2Name) 
nTAD <- length(allTADs)
nParaPairs <- nrow(closePairs)
nRandPairs <- sum(sapply(sampClosePairs, nrow))


# check if "one2one ortholog" is a good approximation for duplication age
seqSimDF <- data.frame(
    "species"= rep(factor(orgStr2Name, orgStr2Name), each=nParaPairs),
    "one2one" = c(closePairs[,"mmusculus_one2one"], closePairs[,"cfamiliaris_one2one"]),
    "commonOrtholog" = factor(c(closePairs[,"mmusculus_commonOrtholg"], closePairs[,"cfamiliaris_commonOrtholg"]), levels=c(TRUE, FALSE)),
    "Ds"= rep(closePairs[,"hsapiens_paralog_ds"],nSpecies),
    "Dn"= rep(closePairs[,"hsapiens_paralog_dn"],nSpecies)
)

pdf(paste0(outPrefix,"old_vs_jung.Ds_vs_one2one.violin_plot.pdf"), w=3.5, h=7)
    ggplot(seqSimDF, aes(x=one2one, y=Ds, fill=one2one)) + geom_violin() + scale_y_log10() + 
        theme_bw() + theme(legend.position = "bottom") + scale_fill_manual(values=COL_AGE) +
        geom_boxplot(aes(fill=NULL), width=.2, colour="black", outlier.shape = NA) +
        facet_grid(.~species) + 
        labs(y="Rate of synonymouse mutations Ds") + guides(fill=guide_legend(title="One-to-one orthologs"))
dev.off()

pdf(paste0(outPrefix,"old_vs_jung.Dn_vs_one2one.violin_plot.pdf"), w=3.5, h=7)
    ggplot(seqSimDF, aes(x=one2one, y=Dn, fill=one2one)) + geom_violin() + scale_y_log10() + 
        theme_bw() + theme(legend.position = "bottom") + scale_fill_manual(values=COL_AGE) +
        geom_boxplot(aes(fill=NULL), width=.2, colour="black", outlier.shape = NA) +
        facet_grid(.~species) + 
        labs(y="Rate of non-synonymouse mutations Ds") + guides(fill=guide_legend(title="One-to-one orthologs"))
dev.off()

# now the same for commonOrtholog
pdf(paste0(outPrefix,"old_vs_jung.Ds_vs_commonOrtholog.violin_plot.pdf"), w=3.5, h=7)
    ggplot(seqSimDF, aes(x=commonOrtholog, y=Ds, fill=commonOrtholog)) + geom_violin() + scale_y_log10() + 
        theme_bw() + theme(legend.position = "bottom") + scale_fill_manual(values=COL_AGE) +
        geom_boxplot(aes(fill=NULL), width=.2, colour="black", outlier.shape = NA) +
        facet_grid(.~species) + 
        labs(y="Rate of synonymouse mutations Ds") + guides(fill=guide_legend(title="Common Ortholog"))
dev.off()

pdf(paste0(outPrefix,"old_vs_jung.Dn_vs_commonOrtholog.violin_plot.pdf"), w=3.5, h=7)
    ggplot(seqSimDF, aes(x=commonOrtholog, y=Dn, fill=commonOrtholog)) + geom_violin() + scale_y_log10() + 
        theme_bw() + theme(legend.position = "bottom") + scale_fill_manual(values=COL_AGE) +
        geom_boxplot(aes(fill=NULL), width=.2, colour="black", outlier.shape = NA) +
        facet_grid(.~species) + 
        labs(y="Rate of non-synonymouse mutations Ds") + guides(fill=guide_legend(title="Common Ortholog"))
dev.off()


# define different markers for age
ageMarks <- c("mmusculus_NotOne2one", "cfamiliaris_NotOne2one", "mmusculus_commonOrtholg", "cfamiliaris_commonOrtholg")

for (AGE_MARK in ageMarks){
    
    print(AGE_MARK)
    
    # make data.frame with additional columns indicating the groups
    paraExpCorDF <- data.frame(
        expCor = unlist(paraExpCorSummaryList), 
        dataset = rep(
            rep(
                paste0(gsub('_', ' ', names(expDFlist)), "\n n=", sapply(expDFlist, ncol), " tissues" )
                , each=2), 
            times = sapply(paraExpCorSummaryList,length)
            ),
        genepairs = rep(
            rep(c("Paralog", "Sampled"), nExp),
            sapply(paraExpCorSummaryList,length)
            ),
        dist = rep(c(allCisPairs[,"dist"], sampCisPairsCombined[,"dist"]), nExp),
        DupAge = factor(rep(c(allCisPairs[,AGE_MARK], sampCisPairsCombined[,AGE_MARK]), nExp), c(TRUE, FALSE), c("Young", "Old"))
    )

    # make data.frame with additional columns indicating the groups
    paraEhDF <- data.frame(
        commonEnhancer = c(closePairs[,"commonEnhancer"], sampDistEhPairsCloseCombined[,"commonEnhancer"]),
        genepairs = rep(
            c("Paralog", "Sampled"),
            c(nrow(closePairs), nrow(sampDistEhPairsCloseCombined))
            ),
        DupAge = factor(c(closePairs[,AGE_MARK], sampDistEhPairsCloseCombined[,AGE_MARK]), c(TRUE, FALSE), c("Young", "Old"))
    )
    paraEhDF$HasCommonEh <- factor(paraEhDF$commonEnhancer >= 1, c(TRUE, FALSE), c("Common enhancer", "No common enhancer"))        
    
    #-------------------------------------------------------------------
    # Expression correlation
    #-------------------------------------------------------------------
    g <- ggplot(paraExpCorDF, na.rm=TRUE, aes(x=genepairs, y = expCor, fill=DupAge)) + theme_bw() + labs(y="Pearson correlation coefficient", title= "Co-expression correlation over tissues")  + facet_grid(.~dataset) + theme(legend.position = "bottom") + theme(text = element_text(size=15)) + scale_color_manual(values=COL) + scale_fill_manual(values=COL_AGE) + geom_boxplot(lwd=1)    
    ggsave(paste0(outPrefix,"old_vs_jung.", AGE_MARK, ".Expression_correlation.all_data.boxplot.pdf"))
    
    g <- ggplot(paraExpCorDF, na.rm=TRUE, aes(x=genepairs, y = expCor, fill=DupAge)) + theme_bw() + labs(y="Pearson correlation coefficient", title= "Co-expression correlation over tissues")  + facet_grid(.~dataset) + theme(legend.position = "bottom") + theme(text = element_text(size=15)) + scale_color_manual(values=COL) + scale_fill_manual(values=COL_AGE) + geom_violin(adjust = .2)
    ggsave(paste0(outPrefix,"old_vs_jung.", AGE_MARK, ".Expression_correlation.all_data.violinplot.pdf"))

    #-------------------------------------------------------------------
    # sahred Enhancer
    #-------------------------------------------------------------------
#~     g <- ggplot(paraEhDF, aes(commonEnhancer, group=DupAge, fill=DupAge)) + 
#~         geom_bar(aes(y = ..density..), width=.5, position="dodge") + xlim(0,10) +
#~         theme_bw() + 
#~         labs(y="Density", title= "Number of shared enhancer")  + theme(legend.position = "bottom") + scale_color_manual(values=COL) + scale_fill_manual(values=COL_AGE) + 
#~         facet_grid(genepairs~.)
#~     ggsave(paste0(outPrefix,"old_vs_jung.", AGE_MARK, ".close.sharedEnhancer.pdf"))
#~ 
#~     g <- ggplot(paraEhDF, aes(HasCommonEh, group=DupAge, fill=DupAge)) + 
#~         geom_bar(aes(y = ..density..), width=.5, position="dodge")  +
#~         theme_bw() + 
#~         labs(y="Density", title= "Number of shared enhancer")  + theme(legend.position = "bottom") + scale_color_manual(values=COL) + scale_fill_manual(values=COL_AGE) + 
#~         facet_grid(genepairs~.)
#~     ggsave(paste0(outPrefix,"old_vs_jung.", AGE_MARK, ".close.sharedEnhancer.pdf"))

    #-------------------------------------------------------------------
    # Distance
    #-------------------------------------------------------------------
    g <- ggplot(paraExpCorDF[1:(nrow(paraExpCorDF)/nExp),], aes(abs(dist)/10^3, group=DupAge, colour=DupAge, fill=DupAge)) + 
        geom_density(alpha = 0.1, lwd=1.2) +
        theme_bw() + 
        labs(y="Density", x="Distance [kb]")  + theme(legend.position = "bottom") + scale_color_manual(values=COL_AGE) + scale_fill_manual(values=COL_AGE) +
        facet_grid(genepairs~.)
    ggsave(paste0(outPrefix,"old_vs_jung.", AGE_MARK, ".dist.pdf"))
    
    g <- ggplot(paraExpCorDF[1:(nrow(paraExpCorDF)/nExp),], aes(abs(dist), group=DupAge, colour=DupAge, fill=DupAge)) + 
        scale_x_log10() +
    #~     geom_histogram() +
        geom_density(alpha = 0.1, lwd=1.2) +
        theme_bw() + 
        labs(y="Density", x="Distance [bp]")  + theme(legend.position = "bottom") + scale_color_manual(values=COL_AGE) + scale_fill_manual(values=COL_AGE) +
        facet_grid(genepairs~.)
    ggsave(paste0(outPrefix,"old_vs_jung.", AGE_MARK, ".dist.log10.pdf"))
    
    #-------------------------------------------------------------------
    # same TAD
    #-------------------------------------------------------------------
    tadDF <- data.frame(
        "sameTAD" = c(
                unlist(rbind(closePairs[,names(allTADs)])),
                unlist(rbind(sampClosePairsCombined[,names(allTADs)]))
                    ),
        "genepairs" = c(rep("Paralog", nParaPairs*nTAD), rep("Sampled", nRandPairs*nTAD)),    
        "cellType" = c(rep(names(allTADs), each=nParaPairs), rep(names(allTADs), each=nRandPairs)),    
        DupAge = factor(c(rep(closePairs[,AGE_MARK], nTAD), rep(sampClosePairsCombined[,AGE_MARK], nTAD)), c(TRUE, FALSE), c("Young", "Old"))
    )

    # get the fraction for each combination
    fracTadDF <- ddply(tadDF, .(genepairs, cellType, DupAge), summarize, frac=percentTrue(sameTAD))
    
    g <- ggplot(fracTadDF, aes(y=frac, x=cellType, fill=DupAge)) + 
        geom_bar(stat="identity", position="dodge") +
        theme_bw() + 
        labs(y="Percent of pairs in same TAD", x="")  + theme(legend.position = "bottom") +  scale_fill_manual(values=COL_AGE) +
        geom_text(aes(label=signif(frac,3)), position=position_dodge(width=1), vjust=-0.25, size=3) + 
        facet_grid(genepairs~.) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(paste0(outPrefix,"old_vs_jung.", AGE_MARK, ".sameTAD.pdf"))

    #-----------------------------------------------------------------------
    # Hi-C and Capture-C of distal pairs
    #-----------------------------------------------------------------------
    nParaDist <- nrow(distalPairs)
    nRandDist <- nrow(sampDistalPairsCombined)
    
    hicDF = data.frame(
            group = c(rep("paralogs", nParaDist), rep("sampled", nRandDist)),
            DupAge = factor(c(distalPairs[,AGE_MARK], sampDistalPairsCombined[,AGE_MARK]), c(TRUE, FALSE), c("Young", "Old")),
            HiC = c(distalPairs[,"HiC"], sampDistalPairsCombined[,"HiC"]),
            HiCobsExp = c(distalPairs[,"HiCobsExp"], sampDistalPairsCombined[,"HiCobsExp"]),
            captureC_raw = c(distalPairs[,"captureC_raw"], sampDistalPairsCombined[,"captureC_raw"]),
            captureC_ObsExp = c(distalPairs[,"captureC_ObsExp"], sampDistalPairsCombined[,"captureC_ObsExp"]),
            dist=abs(c(distalPairs[,"dist"], sampDistalPairsCombined[,"dist"]))/10^3,
            distBin=as.factor(breaks[.bincode(abs(c(distalPairs[,"dist"], sampDistalPairsCombined[,"dist"]))/10^3, breaks)])
    )
    hicDF$HiCNoZero = hicDF$HiC
    hicDF$HiCNoZero[hicDF$HiC == 0] = NA
    hicDF$HiCobsExpNoZero = hicDF$HiCobsExp
    hicDF$HiCobsExpNoZero[hicDF$HiCobsExp == 0] = NA

    g <- ggplot(hicDF, aes(x=group, y=HiC, colour = DupAge))  +
        geom_boxplot(lwd=1.5) + scale_y_log10() + 
        scale_color_manual(values=COL_AGE, name="") +
        theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + 
        guides(fill=guide_legend(title="")) +
        labs(y="Hi-C counts", x="")
    ggsave(paste0(outPrefix,"old_vs_jung.", AGE_MARK, ".Hi-C.pdf"), w=3.5, h=7)

    g <- ggplot(hicDF, aes(x=group, y=HiCobsExp, colour = DupAge))  +
        geom_boxplot(lwd=1.5) + scale_y_log10() + 
        scale_color_manual(values=COL_AGE, name="") +
        theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + 
        guides(fill=guide_legend(title="")) +
        labs(y="Normalized Hi-C counts", x="")
    ggsave(paste0(outPrefix,"old_vs_jung.", AGE_MARK, ".Hi-C_norm.pdf"), w=3.5, h=7)

    g <- ggplot(hicDF, aes(x=group, y=captureC_raw, colour = DupAge))  +
        geom_boxplot(lwd=1.5) + scale_y_log10() + 
        scale_color_manual(values=COL_AGE, name="") +
        theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + 
        guides(fill=guide_legend(title="")) +
        labs(y="Capture-C counts", x="")
    ggsave(paste0(outPrefix,"old_vs_jung.", AGE_MARK, ".CaptureC_raw.pdf"), w=3.5, h=7)
    
    g <- ggplot(hicDF, aes(x=group, y=captureC_ObsExp, colour = DupAge))  +
        geom_boxplot(lwd=1.5) + scale_y_log10() + 
        scale_color_manual(values=COL_AGE, name="") +
        theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + 
        guides(fill=guide_legend(title="")) +
        labs(y="Capture-C (Obs/Exp)", x="")
    ggsave(paste0(outPrefix,"old_vs_jung.", AGE_MARK, ".captureC_ObsExp.pdf"), w=3.5, h=7)

}

#=======================================================================
# Best candidate picking...
#=======================================================================
#~ candidates = with(closePairs, commonEnhancer>3 & sameStrand & abs(ENCODE_cell_lines_expCor) > .6 ) & closePairsGR$stable_TADs 
#~ table(candidates) 

#closePairs[candidates,]

# write the data of Fetuin-A and Fetuin-B (ENSG00000145192, ENSG00000090512)
#grep("ENSG00000145192", closePairs[,2])

#~ write.table(t(as.matrix(closePairs[515,])), file=paste0(outPrefix, ".Fetuin_AB.annotation.txt"), sep="\t", quote=FALSE, col.names=FALSE)


#=======================================================================
# PRC1 paralogs and candidate gene families
#=======================================================================

candFams = c("CBX", "PCGF")

for (FAM in candFams){

    famGR <- tssGR[grep(FAM, tssGR$hgnc_symbol)]
    
    setPairs <- data.frame(t(combn(names(famGR),2 )))
    setPairs <- addPairDist(setPairs, tssGR)
    # add HGNC symbols
    setPairs = addHGNC(setPairs, tssGR)
    # add same Strand info:
    setPairs = addSameStrand(setPairs, tssGR)
    # add number of common enhancers
    setPairs = addCommonEnhancer(setPairs, gene2ehID)
    
    setPairs = addPairExp(setPairs, expDFlist[["ENCODE_cell_lines"]], expCol="IMR_90", label="exp_IMR90")
    
    # add pairwise correlations of gene expression over all tissues
    for (expName in names(expDFlist)) {
        
        message(paste("INFO: annotate pairs with expression correlation form:", expName))
        expDF = expDFlist[[expName]]
        
        setPairs = addCor(setPairs, expDF, colName=paste0(expName, "_expCor"))
    }
    
    # make GRanges objects for cis paralog pairs and random paris on same chromosome
    setPairsGR = getPairAsGR(setPairs, tssGR)
    for(tadName in names(allTADs)){
        message(paste("INFO: Compute overlap with TADs from:", tadName))
        # co-occurance within the same domain
        setPairsGR = addWithinSubject(setPairsGR, allTADs[[tadName]], tadName)
    }
    
    # assign annotation in GRanges object to gene pair data.frames
    setPairs[,names(allTADs)] <- data.frame( mcols(setPairsGR)[, names(allTADs)] )
    
    for(tadName in names(allTADs)){
        message(paste("INFO: Compute common subset of overallping TADs from:", tadName))
        setPairs = addSubTADmode(setPairs, allTADs[[tadName]], tssGR, paste0(tadName, "_subTAD"))
    }
    
    candPairs <- setPairs[!is.na(setPairs$dist) & abs(setPairs$dist) <= 10^6,]
    
    write.table(candPairs, col.names=TRUE, row.names=FALSE, file=paste0(outPrefix, ".candPairs.", FAM, ".annotation.txt"), sep="\t", quote=FALSE)

}


#=======================================================================
# save workspace image
#=======================================================================
#~ save.image(WORKIMAGE_FILE)
#INFO: FANTOM enhancer map: 65953 of 66942 enhancer-associated genes could be mapped uniquelly to ENSG ID
