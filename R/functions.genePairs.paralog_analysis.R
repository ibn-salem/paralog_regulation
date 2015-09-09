########################################################################
#
#   This script implements pipeline like analysis functions and is used
#   by 'paralog_regulation.R'
#   It depends on functions in 'functions.genePairs.R'
# 
########################################################################

#-----------------------------------------------------------------------
# filters paralog paris and returns a list with the variables
#-----------------------------------------------------------------------
filterParalogPairs <- function(paralogPairs, tssGR, similarity="hsapiens_paralog_perc_id"){

    # remove double entries of the form A-B and B-A
    paralogPairsUniqP = uniquePair(paralogPairs)
    
    # get for each gene only one unique pair, the one with highest similarity
    # this is computed by an maximum weight matching
    paralogPairsUniqG = uniquePairPerGeneBySim(paralogPairs, paralogPairs[,similarity])
    
    # get only a unique pair order (one of A-B, B-A) form the unique pairs
    paralogPairsUniqGuniqP = uniquePair(paralogPairsUniqG)
    
    # subset of paralog pairs that are located on the same chromosome
    allCisPairs = getCisPairs(paralogPairsUniqGuniqP, tssGR)
    
    return(list(
        "paralogPairs"=paralogPairs,
        "paralogPairsUniqP"=paralogPairsUniqP,
        "paralogPairsUniqG"=paralogPairsUniqG,
        "paralogPairsUniqGuniqP"=paralogPairsUniqGuniqP,
        "allCisPairs"=allCisPairs
    ))

}

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
#~     paralogPairsUniqG = uniquePairPerGeneBySim(paralogPairs, paralogPairs[,pairScoreCol])
    paralogPairsWithDN = paralogPairs[!is.na(paralogPairs[,pairScoreCol]),]
    paralogPairsUniqG = uniquePairPerGeneBySim(paralogPairsWithDN, pairScoreFac*paralogPairsWithDN[,pairScoreCol])
    
    # get only a unique pair order (one of A-B, B-A) form the unique pairs
    paralogPairsUniqGuniqP = uniquePair(paralogPairsUniqG)
    
    # subset of paralog pairs that are located on the same chromosome
    allCisPairs = getCisPairs(paralogPairsUniqGuniqP, tssGR)
    
    
    
    # Annotate:  add linear distance between TSS
    allCisPairs = addPairDist(allCisPairs, tssGR)
    
    
    # get close cis pairs
    cisPairs = allCisPairs[abs(allCisPairs$dist) <= MAX_DIST,] 
    
    # get distal pairs
    distalCisPairs = allCisPairs[abs(allCisPairs$dist) > DISTAL_MIN_DIST & abs(allCisPairs$dist) <= DISTAL_MAX_DIST,] 
    #-----------------------------------------------------------------------
    # paralog pair filtering numbers
    #-----------------------------------------------------------------------
    nPairs = c(
        "paralogPairs"=nrow(paralogPairs), 
        "paralogPairsUniqP"=nrow(paralogPairsUniqP), 
        "paralogPairsUniqG"=nrow(paralogPairsUniqG), 
        "paralogPairsUniqGuniqP"=nrow(paralogPairsUniqGuniqP),
        "allCisPairs"=nrow(allCisPairs),
        "cisPairs"=nrow(cisPairs),
        "distalCisPairs"=nrow(distalCisPairs)
        )
    write.table(nPairs, file=paste0(outPrefix, ".paralog_pairs_filtering.txt"),
        sep="\t", quote=FALSE, col.names=FALSE)
    
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
    # Sample cis and dist pairs according to linear distance in paralog gene pairs
    #-----------------------------------------------------------------------
    message("Start to sample random gene pairs...")

    # First sample randomly from all gnes
    randPairs = replicate(N_RAND, getRandomPairs(nrow(paralogPairsUniqGuniqP), names(tssGR)), simplify=FALSE)
        
    # filter for random pairs in Cis
    randPairsInCis = lapply(randPairs, getCisPairs, tssGR)
    randPairsInCis = lapply(randPairsInCis, addPairDist, tssGR)
    
    # get all possible gene pairs within MAX_DIST bp
    allGenePairs = getAllGenePairs(tssGR, maxDist=MAX_DIST)
    
    # get sample weights according distance
    cisWeights = getSampleWeightsByDist(allGenePairs, cisPairs, adjust=DENSITY_BW_ADJUST)
    
    # sample gene pairs
    randCisPairs = replicate(N_RAND, 
        sampleFromAllPairsByWeight(n=nrow(cisPairs), hitDF=allGenePairs, tssGR, weight=cisWeights)
        , simplify=FALSE)
    
    # Now sample from all possible gene pairs within DISTAL_MIN_DIST - DISTAL_MAX_DIST bp
    allDistalGenePairs = getAllGenePairs(tssGR, maxDist=DISTAL_MAX_DIST, minDist=DISTAL_MIN_DIST)
    
    # get sample weights according to distance
    distalWeight = getSampleWeightsByDist(allDistalGenePairs, distalCisPairs, adjust=DENSITY_BW_ADJUST)
    
    # sample according to distance
    randDistalCisPairs = replicate(N_RAND, 
        sampleFromAllPairsByWeight(n=nrow(distalCisPairs), hitDF=allDistalGenePairs, tssGR, weight=distalWeight)
        , simplify=FALSE)
    
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
    pdf(paste0(outPrefix, ".linear_distance.random.hist.pdf"))
        paraDist = abs(cisPairs$dist)
        randDist = abs(do.call("rbind", randPairsInCis)[,"dist"])
        randDist = randDist[randDist <= MAX_DIST]
        
        par(cex=1.5, lwd=2, mfrow=c(2,1))
        hist(paraDist/10^3, 50, xlim=c(0, 1000), col=COL_RAND[1], xlab="Distance (kb)", main="Paralog pairs")
        hist(randDist/10^3, 50, xlim=c(0, 1000), col=COL_RAND[2], xlab="Distance (kb)", main="Random Pairs")
    dev.off()

    pdf(paste0(outPrefix, ".linear_distance.sampled.hist.pdf"))
        par(cex=1.5, lwd=2, mfrow=c(2,1))
        hist(abs(cisPairs$dist)/10^3, 50, col=COL[1], xlab="Distance (kb)", main="Paralog pairs")
        hist(abs(do.call("rbind", randCisPairs)[,"dist"])/10^3, 50, col=COL[2], xlab="Distance (kb)", main="Sampled Pairs")
    dev.off()

    #-----------------------------------------------------------------------
    # Analyse TADs
    message("Start to analyse TADs...")

    tadName = paste0(tissueName, "_TAD")
    
    # make GRanges objects for cis paralog pairs and random paris on same chromosome
    cisParaGR = getPairAsGR(cisPairs, tssGR)
    cisRandGR = lapply(randCisPairs, getPairAsGR, tssGR)
    
    # co-occurance within the same domain
    cisParaGR = addWithinSubject(cisParaGR, TAD, tadName)
    cisRandGR = lapply(cisRandGR, addWithinSubject, TAD, tadName)
    
    
    # Plot fraction of paralog pairs within TAD
    for ( D in DIST_TH){
        
        # get boolean vectors
        paraInTAD = mcols(cisParaGR)[width(cisParaGR)<=D, tadName]
        randInTAD = lapply(cisRandGR, function(gr) mcols(gr)[width(gr)<=D, tadName])
    
        # create contingency table and run Fisher test
        contab = rbind(
            para=table(paraInTAD),  
            rand=table(unlist(randInTAD))
        )
        fs.test = fisher.test(contab)
        
        # get percent values for barplot
        paraPercent = percentTrue(paraInTAD)
        randPercent = sapply(randInTAD, percentTrue) 
        heights = c(paraPercent, mean(randPercent, na.rm=TRUE))
        sdRand = sd(randPercent, na.rm=TRUE)
    
        pdf(paste0(outPrefix, ".paralogPairs_cis_within_", D/10^3, "_kb_within_", tadName, ".barplot.pdf"), width=3.5)
    
            par(cex=1.5, lwd=2)
            yMax = 1.3*max(heights) + sdRand
            
            bp = my.barplot(heights, addValues=TRUE, yMax=yMax,
                names=c("Paralogs", "Sampled pairs"),
                main=paste("Pairs (within", D/10^3, "kb)\n in", tadName), ylab="Gene pairs with TSSs in same TAD [%]", col=COL)
    
            error.bar(bp,heights, c(NA,  sdRand), lwd=2)
            
            # pvalue matrix for only two group comparision
            pval_mat = matrix(c(NA, signif(fs.test$p.value, 3)), 1)
            add_pval_all_allPairs(bp, ymax=1.1*max(heights)+sdRand, pval_mat)
    
        dev.off()
    
    }

    #-----------------------------------------------------------------------
    # Analyse HiC ontact frequencies
    message("Start to analyse Hi-C contacts...")
    
    # combine all sampling replicates to one data frame
    randDistalCisPairsCombined = do.call("rbind", randDistalCisPairs)
    
    
    # Adds Hi-C contact frequencies to a gene pari data set
    distalCisPairs = addHiCfreq(distalCisPairs, tssGR, HiClist)
    randDistalCisPairsCombined = addHiCfreq(randDistalCisPairsCombined, tssGR, HiClist)
    #~ randCisPairs = lapply(randCisPairs, addHiCfreq, tssGR, HiClist)
    
    distalCisPairs = addHiCfreq(distalCisPairs, tssGR, HiClistNorm, label="HiCnorm")
    randDistalCisPairsCombined = addHiCfreq(randDistalCisPairsCombined, tssGR, HiClistNorm, label="HiCnorm")

    # save or load downloaded data 
    #save(distalCisPairs, randDistalCisPairsCombined, file=paste0(outDataPrefix, ".Hi-C.", CELL, ".", HIC_RESOLUTION, ".distal_and_sampled_pairs.RData"))
    # Make data.frame for plotting with ggplot2

    breaks = seq(DISTAL_MIN_DIST, DISTAL_MAX_DIST, length.out=10) / 10^3
    
    plotDF = data.frame(
            group = c(rep("paralogs", nrow(distalCisPairs)), rep("sampled", nrow(randDistalCisPairsCombined))),
            HiCraw = c(distalCisPairs[,"HiCfreq"], randDistalCisPairsCombined[,"HiCfreq"]),
            HiCnorm = c(distalCisPairs[,"HiCnorm"], randDistalCisPairsCombined[,"HiCnorm"]),
            dist=abs(c(distalCisPairs[,"dist"], randDistalCisPairsCombined[,"dist"]))/10^3,
            distBin=as.factor(breaks[.bincode(abs(c(distalCisPairs[,"dist"], randDistalCisPairsCombined[,"dist"]))/10^3, breaks)])
        )
    plotDF$HiCrawNoZero = plotDF$HiCraw
    plotDF$HiCrawNoZero[plotDF$HiCraw == 0] = NA
    plotDF$HiCnormNoZero = plotDF$HiCnorm
    plotDF$HiCnormNoZero[plotDF$HiCnorm == 0] = NA

    #------------------------------------------------------------------------
    xlabels=paste0(c("Paralogs", "Sampled"), 
        "\n n=", applyToSubset(plotDF, function(v) sum(!is.na(v)), "HiCraw", "group"), 
        "\n med=", signif(applyToSubset(plotDF, median, "HiCraw", "group", na.rm=TRUE), 3),
        "\n avg=", signif(applyToSubset(plotDF, mean, "HiCraw", "group", na.rm=TRUE), 3)
        )
    ws.test = wilcox.test(plotDF[,"HiCraw"] ~ plotDF[,"group"])
    
    p = ggplot(plotDF, aes(x=group, y=HiCraw, color=group)) + scale_y_log10() +
        geom_boxplot(lwd=1.5) + 
        scale_color_manual(values=COL, guide_legend(title = "")) +
        theme_bw() + theme(text = element_text(size=20), legend.position="none") + 
        guides(fill=guide_legend(title="")) +
        labs(y="Hi-C counts", x="", title=paste0("P = ", signif(ws.test$p.value, 3))) +
        scale_x_discrete(labels=xlabels )
    
    # save ggplot as pdf
    ggsave(p, file=paste0(outPrefix, ".distal_pairs.Hi-C_raw_contacts.ggboxplot.pdf"), width=3.5)
    #------------------------------------------------------------------------
    xlabels=paste0(c("Paralogs", "Sampled"), 
        "\n n=", applyToSubset(plotDF, function(v) sum(!is.na(v)), "HiCnorm", "group"), 
        "\n med=", signif(applyToSubset(plotDF, median, "HiCnorm", "group", na.rm=TRUE), 3),
        "\n avg=", signif(applyToSubset(plotDF, mean, "HiCnorm", "group", na.rm=TRUE), 3)
        )
    ws.test = wilcox.test(plotDF[,"HiCnorm"] ~ plotDF[,"group"])
    
    p = ggplot(plotDF, aes(x=group, y=HiCnorm, color=group)) + scale_y_log10() +
        geom_boxplot(lwd=1.5) + 
        scale_color_manual(values=COL, guide_legend(title = "")) +
        theme_bw() + theme(text = element_text(size=20), legend.position="none") + 
        guides(fill=guide_legend(title="")) +
        labs(y="log2(observed / expected) Hi-C counts", x="", title=paste0("P = ", signif(ws.test$p.value, 3))) +
        scale_x_discrete(labels=xlabels )
    
    # save ggplot as pdf
    ggsave(p, file=paste0(outPrefix, ".distal_pairs.Hi-C_normalized_contacts.ggboxplot.pdf"), width=3.5)
    #------------------------------------------------------------------------
    # do the same without zero Hi-C counts
    xlabels=paste0(c("Paralogs", "Sampled"), 
        "\n n=", applyToSubset(plotDF, function(v) sum(!is.na(v)), "HiCrawNoZero", "group"), 
        "\n med=", signif(applyToSubset(plotDF, median, "HiCrawNoZero", "group", na.rm=TRUE), 3),
        "\n avg=", signif(applyToSubset(plotDF, mean, "HiCrawNoZero", "group", na.rm=TRUE), 3)
        )
    ws.test = wilcox.test(plotDF[,"HiCrawNoZero"] ~ plotDF[,"group"])
    
    p = ggplot(plotDF, aes(x=group, y=HiCrawNoZero, color=group)) + scale_y_log10() +
        geom_boxplot(lwd=1.5) + 
        scale_color_manual(values=COL, guide_legend(title = "")) +
        theme_bw() + theme(text = element_text(size=20), legend.position="none") + 
        guides(fill=guide_legend(title="")) +
        labs(y="Hi-C counts", x="", title=paste0("P = ", signif(ws.test$p.value, 3))) +
        scale_x_discrete(labels=xlabels)
    
    # save ggplot as pdf
    ggsave(p, file=paste0(outPrefix, ".distal_pairs.Hi-C_raw_contacts_noZero.ggboxplot.pdf"), width=3.5)
    #------------------------------------------------------------------------
    xlabels=paste0(c("Paralogs", "Sampled"), 
        "\n n=", applyToSubset(plotDF, function(v) sum(!is.na(v)), "HiCnormNoZero", "group"), 
        "\n med=", signif(applyToSubset(plotDF, median, "HiCnormNoZero", "group", na.rm=TRUE), 3),
        "\n avg=", signif(applyToSubset(plotDF, mean, "HiCnormNoZero", "group", na.rm=TRUE), 3)
        )
    ws.test = wilcox.test(plotDF[,"HiCnormNoZero"] ~ plotDF[,"group"])
    
    p = ggplot(plotDF, aes(x=group, y=HiCnormNoZero, color=group)) + scale_y_log10() +
        geom_boxplot(lwd=1.5) + 
        scale_color_manual(values=COL, guide_legend(title = "")) +
        theme_bw() + theme(text = element_text(size=20), legend.position="none") + 
        guides(fill=guide_legend(title="")) +
        labs(y="log2(observed / expected) Hi-C counts", x="", title=paste0("P = ", signif(ws.test$p.value, 3))) +
        scale_x_discrete(labels=xlabels )
    
    # save ggplot as pdf
    ggsave(p, file=paste0(outPrefix, ".distal_pairs.Hi-C_normalized_contacts_noZero.ggboxplot.pdf"), width=3.5)
    
#~     # compare contact frequencies 
#~     paraHiCfreq = distalCisPairs$HiCfreq
#~     randHiCfreq = randDistalCisPairsCombined$HiCfreq
#~     
#~     # Wilcoxon-rank-sum test
#~     ws.test = wilcox.test(paraHiCfreq, randHiCfreq)
#~     
#~     pdf(paste0(outPrefix, ".distal_pairs.Hi-C_raw_contacts.boxplot.pdf"))
#~         par(cex=1.5, mgp=c(3,1,0))
#~         boxplot(paraHiCfreq+1, randHiCfreq+1, 
#~             log="y", border=COL, lwd=2, 
#~             main=paste0("Willcoxon p-value: ", signif(ws.test$p.value, 2)))
#~         #points(c(1,2), c(mean(paraHiCfreq+1, na.rm=TRUE), mean(randHiCfreq+1, na.rm=TRUE)), pch=18, col=COL, cex=2)
#~         xlabels=paste(c("Paralog genes", "Sampled genes"), "\n n = ", c(length(paraHiCfreq), length(randHiCfreq)),
#~             "\nmedian: ",
#~             signif(c(median(paraHiCfreq, na.rm=TRUE), median(randHiCfreq, na.rm=TRUE)), 3),
#~             "\nmean: ",
#~             signif(c(mean(paraHiCfreq, na.rm=TRUE), mean(randHiCfreq, na.rm=TRUE)), 3))
#~         axis(1, at = c(1, 2), labels=xlabels, line=3, tick=FALSE)
#~         mtext("Hi-C contacts \n between TSS of distal gene pairs", side=2, line=2, cex=1.5)
#~     dev.off()
#~     
#~     paraHiCnorm = distalCisPairs$HiCnorm
#~     randHiCnorm = randDistalCisPairsCombined$HiCnorm
#~         
#~     # Wilcoxon-rank-sum test
#~     ws.test = wilcox.test(paraHiCnorm, randHiCnorm)
#~     #ttest = t.test(paraHiCnorm, randHiCnorm)
#~     
#~     pdf(paste0(outPrefix, ".distal_pairs.Hi-C_normalized_contacts.boxplot.pdf"))
#~     
#~         par(mgp=c(3,1,0),cex=1.5)
#~         boxplot(paraHiCnorm, randHiCnorm,
#~             border=COL, lwd=2,
#~             ylab="log2(observed / expected) Hi-C contacts", 
#~             main=paste0("Willcoxon p-value: ", signif(ws.test$p.value, 2)))
#~                 #,"\nT-test p-value: ", signif(ttest$p.value, 2)))
#~         xlabels=paste(
#~             c("Paralog genes", "Sampled genes"), 
#~             "\n n = ", 
#~             c(length(paraHiCnorm), length(randHiCnorm)),
#~             "\nmedian: ",
#~             signif(c(median(paraHiCnorm, na.rm=TRUE), median(randHiCnorm, na.rm=TRUE)), 3),
#~             "\nmean: ",
#~             signif(c(mean(paraHiCnorm, na.rm=TRUE), mean(randHiCnorm, na.rm=TRUE)), 3))
#~         axis(1, at = c(1, 2), labels=xlabels, line=3, tick=FALSE)
#~     
#~     dev.off()
#~     
    #boxplot(c(list(cisPairs$HiCfreq), lapply(randCisPairs, function(df) df$HiCfreq)))

    
}
#~ outPrefix=paste0(outPrefix, ".Mouse")
#~ paralogPairs=paralogPairsMouse
#~ pairScoreCol="mmusculus_paralog_dn"
#~ tssGR=tssGRmouse
#~ TAD=speciesTADs[["mmusculus"]]
#~ tissueName="Mouse"
#~ HiClist=speciesHiC[["mmusculus"]][[1]]
#~ HiClistNorm=speciesHiC[["mmusculus"]][[2]]
#~ pairScoreFac=-1


   
