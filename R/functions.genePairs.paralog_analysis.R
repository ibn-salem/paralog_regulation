########################################################################
#
#   This script implements pipeline like analysis functions and is used
#   by 'paralog_regulation.R'
#   It depends on functions in 'functions.genePairs.R'
# 
########################################################################

#-----------------------------------------------------------------------
# runs the basic paralog analysis for given input data. 
# TODO: define COL2 and COL, MAX_DIST and some others...
#-----------------------------------------------------------------------
runBasicParalogAnalysis <- function(outPrefix, paralogPairs, pairScoreCol, tssGR, TAD, tissueName, HiClist, HiClistNorm, pairScoreFac=-1){

    
    #-------------------------------------------------------------------
    # Filtering:
    #-------------------------------------------------------------------
    message("Start to filter paralog pairs...")

    # remove double entries of the form A-B and B-A
    paralogPairsUniqP = uniquePair(paralogPairs)
    
    # get for each gene only one unique pair, the one with highest similarity
    # this is computed by an maximum weight matching
    paralogPairsWithDN = paralogPairs[!is.na(paralogPairs[,pairScoreCol]),]
    paralogPairsUniqG = uniquePairPerGeneBySim(paralogPairsWithDN, pairScoreFac*paralogPairsWithDN[,pairScoreCol])
    
    # get only a unique pair order (one of A-B, B-A) form the unique pairs
    paralogPairsUniqGuniqP = uniquePair(paralogPairsUniqG)
    
    # subset of paralog pairs that are located on the same chromosome
    allCisPairs = getCisPairs(paralogPairsUniqGuniqP, tssGR)
    
    
    #-----------------------------------------------------------------------
    # annotate gene pairs
    #-----------------------------------------------------------------------
    # Annotate:  add linear distance between TSS
    allCisPairs = addPairDist(allCisPairs, tssGR)

    # add same TAD annotation
    tadName = paste0(tissueName, "_TAD")
    
    # make GRanges objects for cis paralog pairs and random paris on same chromosome
    allCisPairsGR = getPairAsGR(allCisPairs, tssGR)
    
    allCisPairsGR = addWithinSubject(allCisPairsGR, TAD, tadName)
    # assign annotation in GRanges object to gene pair data.frames
    allCisPairs[,tadName] <- data.frame( mcols(allCisPairsGR)[, tadName] )
    
    # Adds Hi-C contact frequencies to a gene pair data set
    allCisPairs = addHiCfreq(allCisPairs, tssGR, HiClist)
    allCisPairs = addHiCfreq(allCisPairs, tssGR, HiClistNorm, label="HiCnorm")

    #-----------------------------------------------------------------------
    # Separate pairs by distance
    #-----------------------------------------------------------------------

    # get close cis pairs
    closePairs = allCisPairs[abs(allCisPairs$dist) <= MAX_DIST & abs(allCisPairs$dist) > 0,] 
    
    # get distal pairs
    distalPairs = allCisPairs[abs(allCisPairs$dist) > DISTAL_MIN_DIST & abs(allCisPairs$dist) <= DISTAL_MAX_DIST,] 

    # pairs with zero distance (artefact of annotation?)
    zereDistPairs = allCisPairs[abs(allCisPairs$dist) == 0,]


    #-----------------------------------------------------------------------
    # paralog pair filtering numbers
    #-----------------------------------------------------------------------
    nPairs = c(
        "paralogPairs"=nrow(paralogPairs), 
        "paralogPairsUniqP"=nrow(paralogPairsUniqP), 
        "paralogPairsWithDN"=nrow(paralogPairsWithDN), 
        "paralogPairsUniqG"=nrow(paralogPairsUniqG), 
        "paralogPairsUniqGuniqP"=nrow(paralogPairsUniqGuniqP),
        "allCisPairs"=nrow(allCisPairs),
        "zereDistPairs"=nrow(zereDistPairs),
        "closePairs"=nrow(closePairs),
        "distalPairs"=nrow(distalPairs)
        )
    write.table(nPairs, file=paste0(outPrefix, ".paralog_pairs_filtering.txt"),
        sep="\t", quote=FALSE, col.names=FALSE)


    #-----------------------------------------------------------------------
    # Sample cis and dist pairs according to linear distance in paralog gene pairs
    #-----------------------------------------------------------------------
    message("Start to sample random gene pairs...")

    # First sample randomly from all gnes
    randPairs = bplapply(1:N_RAND, function(x){getRandomPairs(nrow(paralogPairsUniqGuniqP), names(tssGR))})
    Sys.sleep(1) # hack to fix problems with bplapply on MOGON
        
    # filter for random pairs in Cis
    randPairsInCis = lapply(randPairs, getCisPairs, tssGR)
    randPairsInCis = lapply(randPairsInCis, addPairDist, tssGR)
    
    # get all possible gene pairs within MAX_DIST bp
    allCloseGenePairs = getAllGenePairs(tssGR, maxDist=MAX_DIST, minDist=1)
    
    # get sample weights according distance
    closeWeights = getSampleWeightsByDist(allCloseGenePairs, closePairs, adjust=DENSITY_BW_ADJUST)
    
    # sample gene pairs
    sampClosePairs = bplapply(1:N_RAND, function(x){ 
        sampleFromAllPairsByWeight(n=nrow(closePairs), hitDF=allCloseGenePairs, tssGR, weight=closeWeights)
        })
    Sys.sleep(1) # hack to fix problems with bplapply on MOGON
    
    message("INFO: Start to collect all distal gene pairs...")
    # Now sample from all possible gene pairs within DISTAL_MIN_DIST - DISTAL_MAX_DIST bp
    allDistalGenePairs = getAllGenePairs(tssGR, maxDist=DISTAL_MAX_DIST, minDist=DISTAL_MIN_DIST)
    message("INFO: Finish to collect all distal gene pairs.")
    
    # get sample weights according to distance
    distalWeight = getSampleWeightsByDist(allDistalGenePairs, distalPairs, adjust=DENSITY_BW_ADJUST)
    
    # sample according to distance
    sampDistalPairs = bplapply(1:N_RAND, function(x){ 
        sampleFromAllPairsByWeight(n=nrow(distalPairs), hitDF=allDistalGenePairs, tssGR, weight=distalWeight)
        })
    Sys.sleep(1) # hack to fix problems with bplapply on MOGON

    #-----------------------------------------------------------------------
    # annotate sampled gene pairs
    #-----------------------------------------------------------------------

    # make GRanges objects for cis paralog pairs and random paris on same chromosome
    sampClosePairsGR = bplapply(sampClosePairs, getPairAsGR, tssGR)
    Sys.sleep(1) # hack to fix problems with bplapply on MOGON

    sampClosePairsGR = bplapply(sampClosePairsGR, addWithinSubject, TAD, tadName)
    Sys.sleep(1) # hack to fix problems with bplapply on MOGON

    # assign annotation in GRanges object to gene pair data.frames
    for (i in seq(N_RAND)){
        sampClosePairs[[i]][,tadName] <- data.frame( mcols(sampClosePairsGR[[i]])[,tadName] )
    }
    
    # combine all sampling replicates to one data frame
    randPairsCombined <- do.call("rbind", randPairs)
    randPairsInCisCombined <- do.call("rbind", randPairsInCis)
    sampClosePairsCombined <- do.call("rbind", sampClosePairs)
    sampDistalPairsCombined <- do.call("rbind", sampDistalPairs)

    # add Hi-C contact frequencies
    sampDistalPairsCombined = addHiCfreq(sampDistalPairsCombined, tssGR, HiClist)    
    sampDistalPairsCombined = addHiCfreq(sampDistalPairsCombined, tssGR, HiClistNorm, label="HiCnorm")

    #-----------------------------------------------------------------------
    # plot number of genes in each group
    #-----------------------------------------------------------------------
    message("Start with basic plots...")
    paralogs = unique(unlist(paralogPairs[,1:2]))
    nonParalogs = setdiff(names(tssGR), paralogs)
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
        plot(groupSizes, as.vector(countSizes), xlim=c(0,40), type="o", pch=19, xlab="Number of paralogs", ylab="Number of genes", main="Paralog group size distribution", col=COL[1])
    dev.off()

    #-----------------------------------------------------------------------
    # fraction of paralog pairs on same chromosom
    #-----------------------------------------------------------------------
    pdf(paste0(outPrefix, ".paralogPairs_on_same_chrom.barplot.pdf"), width=3.5)
        paraPercent = 100 * nrow(allCisPairs) / nrow(paralogPairsUniqGuniqP)
        randPercent = sapply(randPairsInCis, nrow) * 100 / nrow(randPairs[[1]])
        height=c(paraPercent, mean(randPercent))
        par(cex=1.5, lwd=2)
        bp = my.barplot(height, 
            names=paste0(c("Paralogs\n (n=", "Random pairs\n (n="), c(nrow(paralogPairsUniqGuniqP), paste0(N_RAND, "x",nrow(paralogPairsUniqGuniqP))), c(")", ")")), 
            addValues=TRUE, col=COL_RAND,
            main=paste(tissueName, "gene pairs\n on same\n chromosome"), 
            ylab="%")
        error.bar(bp,height, c(NA, sd(randPercent)))
    dev.off()
    
    #-----------------------------------------------------------------------
    # Linear distance of all cis pairs
    #-----------------------------------------------------------------------

    pdf(paste0(outPrefix, ".sampling_cis_and_distal.pdf"), w=7, h=3.5)
        par(cex=1, lwd=1.5, mfrow=c(2,4))
    
        # sampled close by dist
        hist(abs(closePairs[,"dist"] / 10^3), 20, col=COL[1],
        main="Close cis paralogs", xlab="Distance (kb)")
        hist(abs(sampClosePairsCombined[,"dist"] / 10^3), 20, col=COL[2],
        main="Close cis sampled genes", xlab="Distance (kb)")   
        qqplot(abs(closePairs[,"dist"] / 10^3), abs(sampClosePairsCombined[,"dist"] / 10^3), xlab="Paralogs", ylab="Sampled gene pairs", main="QQ-Plot")
        abline(0,1, col="red")    
        qqplot(abs(closePairs[,"dist"] / 10^3), abs(sampClosePairsCombined[,"dist"] / 10^3), log="xy", xlab="Paralogs", ylab="Sampled gene pairs", main="QQ-Plot\n(log scale)")
        abline(0,1, col="red")    
        
        # sampled distal pairs
        hist(abs(distalPairs[,"dist"] / 10^3), 20, col=COL[1],
        main="Distal cis paralogs", xlab="Distance (kb)")
        hist(abs(sampDistalPairsCombined[,"dist"] / 10^3), 20, col=COL[2],
        main="Sampled by dist", xlab="Distance (kb)")   
        qqplot(abs(distalPairs[,"dist"] / 10^3), abs(sampDistalPairsCombined[,"dist"] / 10^3), xlab="Paralogs", ylab="Sampled gene pairs",main="QQ-Plot")
        abline(0,1, col="red")    
        qqplot(abs(distalPairs[,"dist"] / 10^3), abs(sampDistalPairsCombined[,"dist"] / 10^3), log="xy", xlab="Paralogs", ylab="Sampled by dist and enhancer", main="QQ-Plot\n(log scale)")
        abline(0,1, col="red")    
        
    dev.off()

    # Show bias of equaly sampled genes to paralog genes
    
    # get paralog distances <= MAX_DIST
    paraDist = closePairs[,"dist"]
    # get cispairs and distance form uniformaly random genes
    randDist <- randPairsInCisCombined[abs(randPairsInCisCombined[,"dist"])<= MAX_DIST, "dist"]
    # distance form sampled pairs
    sampledDist = sampClosePairsCombined[,"dist"]
        
    pdf(paste0(outPrefix, ".random_genes_distance.hist.pdf"))
        par(cex=1.5, lwd=2, mfrow=c(2,1))
        hist(abs(paraDist)/10^3, 50, col=COL_RAND[1],
        main="Distance between paralog genes", xlab="Distance (kb)")
        hist(abs(randDist)/10^3, 50, col=COL_RAND[2],
        main="Distance between random genes", xlab="Distance (kb)")    
    dev.off()
    
    pdf(paste0(outPrefix, ".samped_random_para_dist.hist.pdf"))
        par(lwd=2, mfrow=c(3,1))
        par(cex=1.5,  mar=c(3, 4.1, 1.5, 2.1))
        
        hist(abs(paraDist)/10^3, 50, col=COL_RAND[1],
        main="Paralog gene pairs", xlab="")
        hist(abs(randDist)/10^3, 50, col=COL_RAND[2],
        main="Random gene pairs", xlab="")    
        hist(abs(sampledDist)/10^3, 50, col=COL[2],
        main="Sampled gene pairs", xlab="")
        mtext("Distance (kb)", side=1, line=2)    
    dev.off()
    
    #-----------------------------------------------------------------------
    # Analyse TADs
    message("Start to analyse TADs...")

    # Plot fraction of paralog pairs within TAD
    # get boolean vectors
    paraInTAD = closePairs[, tadName]
    sampInTAD = lapply(sampClosePairs, function(gP) gP[, tadName])

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
    sdRand = sd(sampPercent, na.rm=TRUE)

    pdf(paste0(outPrefix, ".paralogPairs.within_", tadName, ".barplot.pdf"), width=3.5)

        par(cex=1.5, lwd=2)
        yMax = 1.3*max(heights) + sdRand
        
        bp = my.barplot(heights, addValues=TRUE, yMax=yMax,
            names=c("Paralogs", "Sampled pairs"),
            main=paste("Pairs in", gsub("_", " ", tadName)), ylab="Gene pairs with TSSs in same TAD [%]", col=COL)

        error.bar(bp,heights, c(NA,  sdRand), lwd=2)
        
        # pvalue matrix for only two group comparision
        pval_mat = matrix(c(NA, signif(fs.test$p.value, 3)), 1)
        add_pval_all_allPairs(bp, ymax=1.1*max(heights)+sdRand, pval_mat)

    dev.off()

    #-----------------------------------------------------------------------
    # Analyse HiC ontact frequencies
    message("Start to analyse Hi-C contacts...")
    
    # Make data.frame for plotting with ggplot2
    breaks = seq(DISTAL_MIN_DIST, DISTAL_MAX_DIST, length.out=10) / 10^3
    
    plotDF = data.frame(
            group = c(rep("paralogs", nrow(distalPairs)), rep("sampled", nrow(sampDistalPairsCombined))),
            HiCraw = c(distalPairs[,"HiCfreq"], sampDistalPairsCombined[,"HiCfreq"]),
            HiCnorm = c(distalPairs[,"HiCnorm"], sampDistalPairsCombined[,"HiCnorm"]),
            dist=abs(c(distalPairs[,"dist"], sampDistalPairsCombined[,"dist"]))/10^3,
            distBin=as.factor(breaks[.bincode(abs(c(distalPairs[,"dist"], sampDistalPairsCombined[,"dist"]))/10^3, breaks)])
        )
    plotDF$HiCrawNoZero = plotDF$HiCraw
    plotDF$HiCrawNoZero[plotDF$HiCraw == 0] = NA
    plotDF$HiCnormNoZero = plotDF$HiCnorm
    plotDF$HiCnormNoZero[plotDF$HiCnorm == 0] = NA

    #------------------------------------------------------------------------
    # Hi-C raw
    xlabels=paste0(c("Paralogs", "Sampled"), 
        "\n n=", applyToSubset(plotDF, function(v) sum(!is.na(v)), "HiCraw", "group"), 
        "\n med=", signif(applyToSubset(plotDF, median, "HiCraw", "group", na.rm=TRUE), 3),
        "\n avg=", signif(applyToSubset(plotDF, mean, "HiCraw", "group", na.rm=TRUE), 3)
        )
    ws.test = wilcox.test(plotDF[,"HiCraw"] ~ plotDF[,"group"])
    
    p = ggplot(plotDF, aes(x=group, y=HiCraw, color=group)) +
        geom_boxplot(lwd=1.5) + scale_y_log10() + annotation_logticks(sides="l") +
        scale_color_manual(values=COL, guide_legend(title = "")) +
        theme_bw() + theme(text = element_text(size=20), legend.position="none") + 
        guides(fill=guide_legend(title="")) +
        labs(y="Hi-C counts", x="", title=paste0("p = ", signif(ws.test$p.value, 3))) +
        scale_x_discrete(labels=xlabels )
    ggsave(p, file=paste0(outPrefix, ".distal_pairs.Hi-C_raw_contacts.ggboxplot.pdf"), width=3.5)

    #------------------------------------------------------------------------
    # Hi-C normed
    xlabels=paste0(c("Paralogs", "Sampled"), 
        "\n n=", applyToSubset(plotDF, function(v) sum(!is.na(v)), "HiCnorm", "group"), 
        "\n med=", signif(applyToSubset(plotDF, median, "HiCnorm", "group", na.rm=TRUE), 3),
        "\n avg=", signif(applyToSubset(plotDF, mean, "HiCnorm", "group", na.rm=TRUE), 3)
        )
    ws.test = wilcox.test(plotDF[,"HiCnorm"] ~ plotDF[,"group"])
    
    p = ggplot(plotDF, aes(x=group, y=HiCnorm, color=group))  +
        geom_boxplot(lwd=1.5) + scale_y_log10(breaks=trans_breaks("log10", function(y) 10^y), labels=trans_format("log10", math_format(10^.x))) +# annotation_logticks(sides="l") + 
        scale_color_manual(values=COL, guide_legend(title = "")) +
        theme_bw() + theme(text = element_text(size=20), legend.position="none") + 
        guides(fill=guide_legend(title="")) +
        labs(y="log2(observed / expected) Hi-C counts", x="", title=paste0("p = ", signif(ws.test$p.value, 3))) +
        scale_x_discrete(labels=xlabels )
    ggsave(p, file=paste0(outPrefix, ".distal_pairs.Hi-C_normalized_contacts.ggboxplot.pdf"), width=3.5)
    #------------------------------------------------------------------------
    # Hi-C raw without zero Hi-C counts
    xlabels=paste0(c("Paralogs", "Sampled"), 
        "\n n=", applyToSubset(plotDF, function(v) sum(!is.na(v)), "HiCrawNoZero", "group"), 
        "\n med=", signif(applyToSubset(plotDF, median, "HiCrawNoZero", "group", na.rm=TRUE), 3),
        "\n avg=", signif(applyToSubset(plotDF, mean, "HiCrawNoZero", "group", na.rm=TRUE), 3)
        )
    ws.test = wilcox.test(plotDF[,"HiCrawNoZero"] ~ plotDF[,"group"])
    
    p = ggplot(plotDF, aes(x=group, y=HiCrawNoZero, color=group))  +
        geom_boxplot(lwd=1.5) + scale_y_log10() + annotation_logticks(sides="l") +
        scale_color_manual(values=COL, guide_legend(title = "")) +
        theme_bw() + theme(text = element_text(size=20), legend.position="none") + 
        guides(fill=guide_legend(title="")) +
        labs(y="Hi-C counts", x="", title=paste0("p = ", signif(ws.test$p.value, 3))) +
        scale_x_discrete(labels=xlabels)
    ggsave(p, file=paste0(outPrefix, ".distal_pairs.Hi-C_raw_contacts_noZero.ggboxplot.pdf"), width=3.5)
    #------------------------------------------------------------------------
    # Hi-C normed without zero Hi-C counts
    xlabels=paste0(c("Paralogs", "Sampled"), 
        "\n n=", applyToSubset(plotDF, function(v) sum(!is.na(v)), "HiCnormNoZero", "group"), 
        "\n med=", signif(applyToSubset(plotDF, median, "HiCnormNoZero", "group", na.rm=TRUE), 3),
        "\n avg=", signif(applyToSubset(plotDF, mean, "HiCnormNoZero", "group", na.rm=TRUE), 3)
        )
    ws.test = wilcox.test(plotDF[,"HiCnormNoZero"] ~ plotDF[,"group"])
    
    p = ggplot(plotDF, aes(x=group, y=HiCnormNoZero, color=group)) +
        geom_boxplot(lwd=1.5) + scale_y_log10(breaks=trans_breaks("log10", function(y) 10^y), labels=trans_format("log10", math_format(10^.x))) + 
        scale_color_manual(values=COL, guide_legend(title = "")) +
        theme_bw() + theme(text = element_text(size=20), legend.position="none") + 
        guides(fill=guide_legend(title="")) +
        labs(y="log2(observed / expected) Hi-C counts", x="", title=paste0("p = ", signif(ws.test$p.value, 3))) +
        scale_x_discrete(labels=xlabels )
    ggsave(p, file=paste0(outPrefix, ".distal_pairs.Hi-C_normalized_contacts_noZero.ggboxplot.pdf"), width=3.5)
    
}
#~ outPrefix=paste0(outPrefix, ".Mouse")
#~ paralogPairs=paralogPairsMouse
#~ pairScoreCol="mmusculus_paralog_ds"
#~ tssGR=tssGRmouse
#~ TAD=speciesTADs[["mmusculus"]]
#~ tissueName="Mouse"
#~ HiClist=speciesHiC[["mmusculus"]][[1]]
#~ HiClistNorm=speciesHiC[["mmusculus"]][[2]]
#~ pairScoreFac=-1


   
