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
runBasicParalogAnalysis <- function(outPrefix, paralogPairs, pairScoreCol, tssGR, genesGR, TAD, tissueName, HiClist, HiClistNorm){

    if (!LOAD_PAIRS) {
    
    
        #-------------------------------------------------------------------
        # Filtering:
        #-------------------------------------------------------------------
        message("Start to filter paralog pairs...")
    
        # remove double entries of the form A-B and B-A
        paralogPairsUniqP = uniquePair(paralogPairs)
    
        # filter out overlapping gene pairs
        nonOVL <- nonOverlappingGenePairs(paralogPairsUniqP, genesGR)
        paralogPairsUniqPnonOVL <- paralogPairsUniqP[nonOVL,]
    
        # write all paralog pairs to output file:
        write.table(paralogPairsUniqPnonOVL, file=paste0(outPrefix, ".paralog_pairs.paralogPairsUniqPnonOVL.txt"),
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    
        # get for each gene only one unique pair, the one with highest similarity
        # this is computed by an maximum weight matching
        paralogPairsUniqPnonOVLWithDS = paralogPairsUniqPnonOVL[!is.na(paralogPairsUniqPnonOVL[,pairScoreCol]),]
        
        pairSelectFact = ifelse(SELECT_OLD_PAIRS, 1, -1)
        paralogPairsUniqPnonOVLUniqG = uniquePairPerGeneBySim(paralogPairsUniqPnonOVLWithDS, pairSelectFact*paralogPairsUniqPnonOVLWithDS[,pairScoreCol])
        
        # get only a unique pair order (one of A-B, B-A) form the unique pairs
        paralogPairsUniq = uniquePair(paralogPairsUniqPnonOVLUniqG)
        
        # subset of paralog pairs that are located on the same chromosome
        allCisPairs = getCisPairs(paralogPairsUniq, tssGR)
        
        #-----------------------------------------------------------------------
        # annotate gene pairs
        #-----------------------------------------------------------------------
        # add same Strand info to all paralogs
        paralogPairsUniq = addSameStrand(paralogPairsUniq, tssGR)
        
        # add linear distance between TSS
        paralogPairsUniq = addPairDist(paralogPairsUniq, tssGR)
        allCisPairs = addPairDist(allCisPairs, tssGR)
        
        # add HGNC symbols
        allCisPairs = addHGNC(allCisPairs, tssGR)
        
        # add same Strand info:
        allCisPairs = addSameStrand(allCisPairs, tssGR)
        
        # make GRanges objects for cis paralog pairs and random paris on same chromosome
        allCisPairsGR = getPairAsGR(allCisPairs, tssGR)
        
        allCisPairsGR = addWithinSubject(allCisPairsGR, TAD, "TAD")
        # assign annotation in GRanges object to gene pair data.frames
        allCisPairs[,"TAD"] <- data.frame( mcols(allCisPairsGR)[, "TAD"] )
        
        # Adds Hi-C contact frequencies to a gene pair data set
        allCisPairs = addHiCfreq(allCisPairs, tssGR, HiClist)
        allCisPairs = addHiCfreq(allCisPairs, tssGR, HiClistNorm, label="HiCnorm")
    
        # save allCisPairs with all annotations
        write.table(allCisPairs, file=paste0(outPrefix, ".paralog_pairs.allCisPairs.annotated.txt"),
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
        
        #-----------------------------------------------------------------------
        # Separate pairs by distance
        #-----------------------------------------------------------------------
    
        # get close cis pairs
        closePairs = allCisPairs[abs(allCisPairs$dist) <= MAX_DIST,] 
        
        # get distal pairs
        distalPairs = allCisPairs[abs(allCisPairs$dist) > MAX_DIST,] 
    
        #-----------------------------------------------------------------------
        # paralog pair filtering numbers
        #-----------------------------------------------------------------------
        nPairs = c(
            "paralogPairs"=nrow(paralogPairs), 
            "paralogPairsUniqP"=nrow(paralogPairsUniqP), 
            "paralogPairsUniqPnonOVL"=nrow(paralogPairsUniqPnonOVL), 
            "paralogPairsUniqPnonOVLWithDS"=nrow(paralogPairsUniqPnonOVLWithDS), 
            "paralogPairsUniqPnonOVLUniqG"=nrow(paralogPairsUniqPnonOVLUniqG), 
            "paralogPairsUniq"=nrow(paralogPairsUniq),
            "allCisPairs"=nrow(allCisPairs),
            "closePairs"=nrow(closePairs),
            "distalPairs"=nrow(distalPairs)
            )
    
        write.table(nPairs, file=paste0(outPrefix, ".paralog_pairs_filtering.txt"),
            sep="\t", quote=FALSE, col.names=FALSE)
    
        message("INFO: Start sampling of gene pairs...")
        #===================================================================
        # Sample gene pairs
        #===================================================================
        
        #-----------------------------------------------------------------------
        # Sample pairs with equal probability from all genes
        #-----------------------------------------------------------------------
        randPairs <- lapply(1:N_RAND, function(x){getRandomPairs(nrow(paralogPairsUniq), names(tssGR))})
        Sys.sleep(3) # hack to fix problems with bplapply on MOGON
        
        # add same strand information
        randPairs <- lapply(randPairs, addSameStrand, tssGR)
        randPairs <- lapply(randPairs, addPairDist, tssGR)
    
        # filter for random pairs in Cis
        randPairsInCis <- lapply(randPairs, getCisPairs, tssGR)
        
        #-----------------------------------------------------------------------
        # Sample all cis pairs by distance using 90 bins 
        #-----------------------------------------------------------------------
        # Now sample from all possible gene pairs within 1 - DISTAL_MAX_DIST bp
        allCisGenePairsOVL <- getAllGenePairs(tssGR, maxDist=DISTAL_MAX_DIST, minDist=1)
        
        # filter out overlapping pairs
        allNonOVL <- nonOverlappingGenePairs(allCisGenePairsOVL, genesGR, useIDs=TRUE)
        allCisGenePairs <- allCisGenePairsOVL[allNonOVL,]
    
        # get sampling weights for distal pairs according to distance
        cisDistBreaks <- seq(log10(1), log10(DISTAL_MAX_DIST), length.out=91)
        cisDistProp <- weightsByBin(log10(abs(allCisPairs$dist)), log10(abs(allCisGenePairs$dist)), breaks=cisDistBreaks)
        
        # sample according to distance
        sampCisPairs <- bplapply(1:N_RAND, function(x){ 
            sampleFromAllPairsByWeight(n=nrow(allCisPairs), hitDF=allCisGenePairs, tssGR, weight=cisDistProp)
            })
        Sys.sleep(3) # hack to fix problems with bplapply on MOGON
    
        #===================================================================
        # Annotate sampled gene pairs
        #===================================================================
    
        # make GRanges objects for cis paralog pairs and random paris on same chromosome
        sampCisPairsGR = bplapply(sampCisPairs, getPairAsGR, tssGR)
        Sys.sleep(3) # hack to fix problems with bplapply on MOGON
    
        sampCisPairsGR = bplapply(sampCisPairsGR, addWithinSubject, TAD, "TAD")
        Sys.sleep(3) # hack to fix problems with bplapply on MOGON
    
        # assign annotation in GRanges object to gene pair data.frames
        for (i in seq(N_RAND)){
            sampCisPairs[[i]][,"TAD"] <- data.frame( mcols(sampCisPairsGR[[i]])[,"TAD"] )
        }
        
        # combine all sampling replicates to one data frame
        randPairsCombined <- do.call("rbind", randPairs)
        randPairsInCisCombined <- do.call("rbind", randPairsInCis)
        sampCisPairsCombined <- do.call("rbind", sampCisPairs)
    
        # add Hi-C contact frequencies
        sampCisPairsCombined = addHiCfreq(sampCisPairsCombined, tssGR, HiClist)
        sampCisPairsCombined = addHiCfreq(sampCisPairsCombined, tssGR, HiClistNorm, label="HiCnorm")
        
        #-----------------------------------------------------------------------
        # save a work image after sampling and annotation.
        #-----------------------------------------------------------------------
        save(
        paralogPairsUniq, 
        allCisPairs, 
        closePairs, 
        distalPairs,
        randPairs,  randPairsCombined,
        randPairsInCis, randPairsInCisCombined,
        sampCisPairs, sampCisPairsCombined,
        file=paste0(WORKIMAGE_FILE, ".annotated_gene_pairs.Rdata"))
    
    }else{
        message(paste("INFO: Start loading annotated gene pairs form this file:",  paste0(WORKIMAGE_FILE, ".annotated_gene_pairs.Rdata")))
        load(paste0(WORKIMAGE_FILE, ".annotated_gene_pairs.Rdata"))
        message("INFO: Finished loading.")
    }
    
    #===================================================================
    # Run basic analysis
    #===================================================================

    #-----------------------------------------------------------------------
    # plot number of genes in each group
    #-----------------------------------------------------------------------
    message("Start with basic plots...")
    paralogs = unique(unlist(paralogPairs[,1:2]))
    nonParalogs = setdiff(names(tssGR), paralogs)
    nrGenes = sapply(list(paralogs, nonParalogs), length)    
    
    pdf(paste0(outPrefix, ".number_of_genes.pdf"), width=3.5)    
        par(cex=1.5, mgp=c(3,1,0), mar = c(7, 4, 4, 2) + 0.1)
        barP = barplot(nrGenes, beside=TRUE, col=COL2, names.arg=NA, ylab="Number of genes", ylim=c(0, 1.25*max(nrGenes)), main=tissueName) 
    
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
            main=paste(tissueName, "gene pairs\n on same\n chromosome"), 
            ylab="Pairs on same chromosome [%]")
        error.bar(bp,height, c(NA, sd(randPercent)))

        # p-value matrix for only two group comparison
        pval_mat = matrix(c(NA, signif(pval, 3)), 1)
        add_pval_all_allPairs(bp, ymax=1.1*max(height)+ sd(randPercent), pval_mat, min_pval=10^-16)


    dev.off()
    
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
    
    #-----------------------------------------------------------------------
    # Show bias of equaly sampled genes to paralog genes
    #-----------------------------------------------------------------------
    
    # get paralog distances <= MAX_DIST
    paraDist = abs(closePairs[,"dist"])/10^3
    # get cispairs and distance form uniformaly random genes
    randDist <- abs(randPairsInCisCombined[abs(randPairsInCisCombined[,"dist"])<= MAX_DIST, "dist"])/10^3
    # distance form sampled pairs
    sampledDist = abs(sampCisPairsCombined[abs(sampCisPairsCombined[,"dist"]) <= MAX_DIST,"dist"])/10^3
    
    randDistPval <- wilcox.test(paraDist, randDist)$p.value
    sampleDistPval <- wilcox.test(paraDist, sampledDist)$p.value
                
    pdf(paste0(outPrefix, ".random_genes_distance.hist.pdf"))
        par(cex=1.5, lwd=2, mfrow=c(2,1))
        hist(paraDist, 50, col=COL_RAND[1],
        main="Distance between paralog genes", xlab="Distance (kb)")
        hist(randDist, 50, col=COL_RAND[2],
        main="Distance between random genes", xlab="Distance (kb)")    
    dev.off()
    
    pdf(paste0(outPrefix, ".samped_random_para_dist.hist.pdf"))
#~         par(lwd=2, mfrow=c(3,1))
#~         par(cex=1.5,  mar=c(3, 4.1, 1.5, 2.1))
#~         
#~         hist(paraDist, 50, col=COL_RAND[1],
#~         main="Paralog gene pairs", xlab="")
#~         hist(randDist, 50, col=COL_RAND[2],
#~         main="Random gene pairs", xlab="")    
#~         hist(sampledDist, 50, col=COL[2],
#~         main="Sampled gene pairs", xlab="")
#~         mtext("Distance (kb)", side=1, line=2)    

        par(lwd=2, mfrow=c(3,1))
        par(cex=1.5,  mar=c(3, 4.1, 1.5, 2.1))
        
        hp <- hist(paraDist, 50, col=COL_RAND[1],
        main=paste("Paralog gene pairs in", tissueName), xlab="")
        hr <- hist(randDist, 50, col=COL_RAND[2],
        main=paste0("Random gene pairs (p=", signif(randDistPval, 2), ")"), xlab="")
        hs <- hist(sampledDist, 50, col=COL[2],
        main=paste0("Sampled gene pairs (p=", signif(sampleDistPval, 2), ")"), xlab="")
        mtext("Distance (kb)", side=1, line=2, cex=1.5)    

    dev.off()
    
    #-----------------------------------------------------------------------
    # Analyse TADs
    message("Start to analyse TADs...")

    # Plot fraction of paralog pairs within TAD
    # get boolean vectors
    paraInTAD = closePairs[, "TAD"]
    sampInTAD = lapply(sampCisPairs, function(gP) gP[abs(gP$dist) <= MAX_DIST, "TAD"])

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

    pdf(paste0(outPrefix, ".paralogPairs.within_TAD.barplot.pdf"), width=3.5)

        par(cex=1.5, lwd=2)
        yMax = 1.3*max(heights) + sdRand
        
        bp = my.barplot(heights, addValues=TRUE, yMax=yMax,
            names=c("Paralogs", "Sampled"),
            main=paste("Pairs in", tissueName, "TAD"), ylab="Gene pairs in same TAD [%]", col=COL)

        error.bar(bp,heights, c(NA,  sdRand), lwd=2)
        
        # pvalue matrix for only two group comparision
        pval_mat = matrix(c(NA, signif(fs.test$p.value, 3)), 1)
        add_pval_all_allPairs(bp, ymax=1.1*max(heights)+sdRand, pval_mat)

    dev.off()

    #-----------------------------------------------------------------------
    # build data.frame with all annotations and all pairs
    #-----------------------------------------------------------------------
    breaksCis = 10^(0:9) / 10^3
    breaksTAD = c(0, 10, 1000, Inf)

    nPairs <- nrow(allCisPairs)
    allCisPairs$g1 <- allCisPairs[,1]
    allCisPairs$g2 <- allCisPairs[,2]
    
    allPairDFs <- list(allCisPairs, sampCisPairsCombined)
    names(allPairDFs) <- c("paralogs", "sampled by distance")
    
    # get common columns in all data sets
    cC <- Reduce(intersect, lapply(allPairDFs, names))
    
    # combine all pairs by using only the common columns
    aP <- do.call("rbind", lapply(allPairDFs, function(df) df[,cC]))
    
    # add noZero HI-C 
    aP <- addNoZero(aP, c("HiC", "HiCnorm"))
    
    
    # add further columns
    aP[,"group"] <- factor(rep(c("paralog", "sampled"), c(nrow(allPairDFs[[1]]), sum(sapply(allPairDFs[2:length(allPairDFs)], nrow)))), levels=c("paralog", "sampled"))
    aP[,"sampType"] <- factor(rep(names(allPairDFs), times=sapply(allPairDFs, nrow)), levels=names(allPairDFs))
    aP[,"replicate"] <- c(rep(1, nPairs), replicate(length(allPairDFs)-1, rep(1:N_RAND, each=nPairs)))
    aP[,"dist"] <- abs(aP[,"dist"]) / 10^3
    aP[,"distGroup"] <- factor(ifelse(aP[,"dist"] <= (MAX_DIST/10^3), "close", "distal"), levels=c("close", "distal"))
    aP[,"distBin"] <- as.factor(breaksCis[.bincode(aP[,"dist"], breaksCis)])
    aP[,"distTadBin"] <- factor(breaksTAD[.bincode(aP[,"dist"], breaksTAD)], levels=breaksTAD[1:3], labels=c("<10kb", "10-1000kb", ">1000kb"))
    
    # save bright allDF data set with column for each source:
    write.table(aP, file=paste0(outPrefix, ".allPairs_broad.csv"),
        sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    
    #-----------------------------------------------------------------------
    # Analyse HiC ontact frequencies
    message("Start to analyse Hi-C contacts...")
    
    for (HiCcol in c("HiC", "HiCnorm", "HiCNoZero", "HiCnormNoZero") ) {
        
        #------------------------------------------------------------------------
        # Hi-C of all pairs
        #------------------------------------------------------------------------
                
        xlabels=paste0(c("Paralogs", "Sampled"), 
            "\n n=", applyToSubset(aP, function(v) sum(!is.na(v)), HiCcol, "group"), 
            "\n med=", signif(applyToSubset(aP, median, HiCcol, "group", na.rm=TRUE), 3),
            "\n avg=", signif(applyToSubset(aP, mean, HiCcol, "group", na.rm=TRUE), 3)
            )
        ws.test = wilcox.test(aP[,HiCcol] ~ aP[,"group"])
    
        p = ggplot(aP, aes(x=group, y=get(HiCcol), color=group)) +
            geom_boxplot(lwd=1.5) + scale_y_log10() + annotation_logticks(sides="l") +
            scale_color_manual(values=COL, guide_legend(title = "")) +
            theme_bw() + theme(text = element_text(size=20), legend.position="none") + 
            guides(fill=guide_legend(title="")) +
            labs(y=HiCcol, x="", title=paste0("p = ", signif(ws.test$p.value, 3))) +
            scale_x_discrete(labels=xlabels )
        ggsave(p, file=paste0(outPrefix, ".all_pairs.", HiCcol, ".by_group.ggboxplot.pdf"), width=3.5)

        #------------------------------------------------------------------------
        # Hi-C of close pairs
        #------------------------------------------------------------------------
        
        # select only close pairs
        subDF <- aP[aP$distGroup == "close",]
        
        xlabels=paste0(c("Paralogs", "Sampled"), 
            "\n n=", applyToSubset(subDF, function(v) sum(!is.na(v)), HiCcol, "group"), 
            "\n med=", signif(applyToSubset(subDF, median, HiCcol, "group", na.rm=TRUE), 3),
            "\n avg=", signif(applyToSubset(subDF, mean, HiCcol, "group", na.rm=TRUE), 3)
            )
        ws.test = wilcox.test(subDF[,HiCcol] ~ subDF[,"group"])
    
        p = ggplot(subDF, aes(x=group, y=get(HiCcol), color=group)) +
            geom_boxplot(lwd=1.5) + scale_y_log10() + annotation_logticks(sides="l") +
            scale_color_manual(values=COL, guide_legend(title = "")) +
            theme_bw() + theme(text = element_text(size=20), legend.position="none") + 
            guides(fill=guide_legend(title="")) +
            labs(y=HiCcol, x="", title=paste0("p = ", signif(ws.test$p.value, 3))) +
            scale_x_discrete(labels=xlabels )
        ggsave(p, file=paste0(outPrefix, ".close.", HiCcol, ".by_group.ggboxplot.pdf"), width=3.5)

        #------------------------------------------------------------------------
        # Hi-C of distal pairs
        #------------------------------------------------------------------------
        
        # select only close pairs
        subDF <- aP[aP$distGroup == "distal",]
        
        xlabels=paste0(c("Paralogs", "Sampled"), 
            "\n n=", applyToSubset(subDF, function(v) sum(!is.na(v)), HiCcol, "group"), 
            "\n med=", signif(applyToSubset(subDF, median, HiCcol, "group", na.rm=TRUE), 3),
            "\n avg=", signif(applyToSubset(subDF, mean, HiCcol, "group", na.rm=TRUE), 3)
            )
        ws.test = wilcox.test(subDF[,HiCcol] ~ subDF[,"group"])
    
        p = ggplot(subDF, aes(x=group, y=get(HiCcol), color=group)) +
            geom_boxplot(lwd=1.5) + scale_y_log10() + annotation_logticks(sides="l") +
            scale_color_manual(values=COL, guide_legend(title = "")) +
            theme_bw() + theme(text = element_text(size=20), legend.position="none") + 
            guides(fill=guide_legend(title="")) +
            labs(y=HiCcol, x="", title=paste0("p = ", signif(ws.test$p.value, 3))) +
            scale_x_discrete(labels=xlabels )
        ggsave(p, file=paste0(outPrefix, ".distal.", HiCcol, ".by_group.ggboxplot.pdf"), width=3.5)

    }
    #------------------------------------------------------------------------
    # plot Hi-C contacts versus genomic distance
    
    p1 <- dotplotWithDensityLogXY(revDF(aP[aP$distGroup == "close",]), "dist", "HiC", "group", COL)
    p2 <- dotplotWithDensityLogXY(revDF(aP[aP$distGroup == "close",]), "dist", "HiCnorm", "group", COL)
    
    p3 <- dotplotWithDensityLogXY(revDF(aP[aP$distGroup == "distal",]), "dist", "HiC", "group", COL)
    p4 <- dotplotWithDensityLogXY(revDF(aP[aP$distGroup == "distal",]), "dist", "HiCnorm", "group", COL)

    pdf(paste0(outPrefix, ".close_and_distal_pairs.Hi-C_vs_dist_all.ggboxplot.pdf"), w=12, h=12)
        do.call(grid.arrange, list(p1,p2,p3,p4)) 
    dev.off()

    p1 <- dotplotWithDensityLogXY(revDF(aP), "dist", "HiC", "group", COL)
    p2 <- dotplotWithDensityLogXY(revDF(aP), "dist", "HiCnorm", "group", COL)
    pdf(paste0(outPrefix, ".all_pairs.Hi-C_vs_dist_all.ggboxplot.pdf"), w=6, h=12)
        do.call(grid.arrange, list(p1,p2)) 
    dev.off()

    
}
#~ outPrefix=paste0(outPrefix, ".Mouse")
#~ paralogPairs=paralogPairsMouse
#~ pairScoreCol="mmusculus_paralog_ds"
#~ tssGR=speciesTssGR[["mmusculus"]]
#~ genesGR=speciesGenesGR[["mmusculus"]]
#~ TAD=speciesTADs[["mmusculus"]]
#~ tissueName="Mouse"
#~ HiClist=speciesHiC[["mmusculus"]][[1]]
#~ HiClistNorm=speciesHiC[["mmusculus"]][[2]]
#~ 


   
