########################################################################
#
# A script to analyse the co-regulation by distal enhancers of paralog 
# genes and functional related genes.
# It looks for association (compared to randomized paralog assignment) 
# for a paralog pair (if paralog group is larger, remove randomly genes) to
# - co-localization on linear genome (distance between them)
# - common enhancers associations in correlation based maps
# - co-occurances in same interaction domains
# - more contact to each other in Hi-C map or more contacts to common enhancers.
########################################################################

#--------------------------#
# TODOs and known issues:  #
#--------------------------#
# DONE:
# X list of paralog pairs contains duplicates (each direction)
# X Some paralog pairs might not be contained in both direction!!!
# X randomization with all genes
# X annotate the GR objects of pairs (is in same domain TREU/FALSE)
# X faster mapping of genes to enhancer IDs (write custom function)
# - randomized ehancer-promoter map
# X random gene pairs with same distance distribution as paralogs
# X take only one unique pair per gene
# x exclude real paralog pairs in randomly sampled control set
# X select unique pair per gene by maximal sequence similarity
# x in sampling process, consider density smooth parameters (might be the reason for slightly to less very close (nearly zero) distance sampled pairs)
#   Try to take wight not from absolut, but relativ distances
# After discussion with Miguel (16.04.15):
# x Find a metric to get exclusive expression (threshold on low exp genes, than correlation?)
# X Use maximal information coefficient (MIC) 
# X Search for TF binding motifs in promoter and enhancer of paralogs 
# X check for co-expression of all genes in same TAD (not paralogs)
# X used observed/expected matrix provided by vectors from Rao et al. 2014
# X Use list of all TADs from Rao14 and Dixon12
# X redesin analysis: 
# X     - make general function for gene-pairs (instead of paralogs)
# X     - Do all annotation on all possible gene pairs before sampling
# => RECHECK: MIC score with toy example and dot plot (R and MIC should correlate)
# - preselect non-coexistence gene pairs and cluster tissues with them
# X To check which distal pair cut-off makes most sense bin pairs by size and make boxplots using ggplot2
# RECHECK: The maximum matching code might overwrite non-symetric similarities in case of A-B B-A pairs
# X check for same strand correlation with distance
# X Use capture Hi-C data (from latest Peter Fraser paper) to quantify contacts between paralogs
#       - recheck parsing and sparse matrix object (manually check function parseCaptureHiC in parseHiC.R)
# # take duplication age (dS) into account
# # for faster pipeline Run: outsource parsing of (and scanning) of TF motifs
# - Recheck sampling of gene pairs by distance.
#       - try to sample separatly for dist and enhancer
#       - plot distal pair distance with sampled pairs with log-log qqplot
#       - Maybe Sample according to distance of allCisPairs and divide close and distal afterwards

# FEDBACK from RECOMB-CG 2015
# - take duplication age into account, separate young and old pairs
# - Are there paralogs with different functions? (non-coexistence expression?)
# - For TAD vs. rearrangement correlation take length of intergenes into account.
# - Direct correlation with syntenic blocks (use Magsimus or similar tools)


# DISCUSSION: (after discussion with Miguel on 16.06.15
# X leave out expression data for the first publication
# X Question of interest: Paralogs fitting to TAD structure of genome?
# x RECHECK carefully mouse and dog Hi-C data (and comparision to sampled genes)
# x change scale of distance box plot of orthologs
# x Synteny breaks of between mouse and human around TAD boundaries.
# X Use >= 1MB as distal pair cut-off
# X linear distance correlation for orthologs of sampled genes (should be less correlation?)
# X Put numbers to Hi-C boxplot in ortholog analysis

# FURTHER INTERESTING STUFF:
# X GTEx Consortium RNA-seq data for expression analysis
# - use other functional gene pairs (KEGG, GO (level?), PPI)
# X Colocalization of paralog pairs in other organisms (% species with shared chrom)
# - Use mouse Hi-C data from the Rao et al. 2014 paper

# EXAMPLE:
# - Check HoxA and HoxD locus in detail, as well as, IGf2/H19 locus (Kurukut et al. PNAS 2006)
# X Use the PRC1 complex as example: 
# X     - check for CBX2,4,8 on chr17 and CBX6,7 on chr22
# X     - PHC1,2,3 are on diff. chrom
# chr17: CBX1,2,8,4

# ADDITIONAL ANALYSIS:
# - repeat all analysis with all pairs, only two pairs
# X CHECK: ratio of synonymous mutations used for pair choosing?
# - check robustness to 1MB distance cutoff
# - build sampled background separately for enhancer, and TAD/Hi-C 
# - Check stable TADs, check definition. Why not significant?
# X Significance test on expression correlation
# - check linear distance conservation with correlation p-value
# - Compute fraction of one-to-one orthologs within the same TAD from 
#   only those human paralogs that are in the same TAD
# - include size of mouse and dog TADs in the size boxplot of all TADs
# - for expression analysis replace boxplot with density plot (similar to Fortin2015)
# - Enhancer positioning pattern around pairs of paralogs (and within TADs)
# X Evolutionary breakpoint of TADs
# - check source of slight enrichment of negative expression correlation (use random pairs from different chromosomes)
# - expression correlation of distal pairs

# ISSUES TO BE FIXED BEFORE FINAL NUMBERS:
# - TAD bed file and ranges +/- 1 error (Rao org data use 0-based including coords)
# - make all gene pair data frames as.character()
# - check mapping of enhancers to genes with tssGR object and id remapping
#   DONE: 65953 of 66942 genes could be mapped uniquelly to ENSG IDs.

# TO FINALIZE PIPELINE FOR SUBMISSION:
# X redesign code to run on MOGON server
# X use orghologMouse instead of orthologAll data set
# X use seed() command for reproducible randomizations
# - remove the following parts from the analysis (if not included in manuscript):
#   X TF motif analysis
#   X conserved TADs over all cell types
#   x MIC score of expression correlation
# - put sameTAD annotations in the closePairs not in the GR and combine replicated sampled pairs early

require(biomaRt)        # to retrieve human paralogs from Ensembl
require(stringr)        # for some string functionality
require(RColorBrewer)   # for nice colors
require(colorspace)     # for some more colors
require(GenomicRanges)  # for genomic intervals and overlap calculation
require(rtracklayer)    # for import.bed
require(plyr)           # count() function
require(data.table)     # for data.table object
require(gridExtra)      # for dotplot with denisty at axis
require(gplots)         # heatmap.2 function
require(ggplot2)        # for nice plots
require(scales)         # for proper logarithmic scales in ggplot
require(BiocParallel)   # for parallel computing

#-----------------------------------------------------------------------
# Load parameters from external script
#-----------------------------------------------------------------------
# read parameter script location from command-line argument
#~ args <- commandArgs(trailingOnly = TRUE)
#~ PARAM_SCRIPT=args[1]

PARAM_SCRIPT="R/paralog_regulation.param.v15.R"
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


#-----------------------------------------------------------------------
# load some custom functions
#-----------------------------------------------------------------------
source("R/functions.plot.R")
source("R/functions.regMap.R")
source("R/functions.GRanges.R")
source("R/functions.genePairs.R")
source("R/functions.genePairs.randomization.R")
source("R/functions.genePairs.paralog_analysis.R")
source("R/functions.Hi-C.R")
source("R/parseHiC.R")
#~ 
#~ # load ensemble data sets of genes and paralog pairs
#~ source("R/data.ensembl.R")     # load ensambl data
#~ source("R/data.expression.R")  # load expression data from EBI expression atlas
#~ source("R/data.captureHiC.R")  # load capture Hi-C data between promoters from Mifsud et al. 2015
#~ source("R/data.gene_age.R")  # load duplication age for paralogs as computed by Pablo Mier Munoz

#=======================================================================
# 1.) Load data from exported tables
#=======================================================================
#TODO

#=======================================================================
# 4.) Run analysis
#=======================================================================

#-----------------------------------------------------------------------
# verify same distribution of linked enhancer in paralog genes and randomly sampled genes
#-----------------------------------------------------------------------
pdf(paste0(outPrefix, ".sampling_cis_dist.pdf"))
    par(cex=1, lwd=1.5, mfrow=c(2,1))
    # verify same distribution of distances
    paraDist = closePairs$dist / 10^3
    sampDist = sampClosePairsCombined$dist / 10^3
    hist(abs(paraDist), 50, col=COL[1],
    main="Paralog pairs", xlab="Distance (kb)")
    hist(abs(sampDist), 50, col=COL[2],
    main="Sampled pairs", xlab="Distance (kb)")    
dev.off()

# get distributions of linked enhancers and distances
paraLinkedEnhancer = tssGR[c(closePairs[,1], closePairs[,2])]$linked_enhancer
sampEhCloseLinkedEnhancer = tssGR[c(sampEhClosePairsCombined[,1], sampEhClosePairsCombined[,2])]$linked_enhancer

pdf(paste0(outPrefix, ".sampling_cis_and_distal.pdf"), w=9, h=10.5)
    par(cex=1, lwd=1.5, mfrow=c(4,4))

    # sampled close by dist
    hist(abs(closePairs[,"dist"] / 10^3), 20, col=COL[1],
    main="Close cis paralogs", xlab="Distance (kb)")
    hist(abs(sampClosePairsCombined[,"dist"] / 10^3), 20, col=COL[2],
    main="Close cis sampled genes", xlab="Distance (kb)")   
    qqplot(abs(closePairs[,"dist"] / 10^3), abs(sampClosePairsCombined[,"dist"] / 10^3), xlab="Paralogs", ylab="Sampled gene pairs", main="QQ-Plot")
    abline(0,1, col="red")    
    qqplot(abs(closePairs[,"dist"] / 10^3), abs(sampClosePairsCombined[,"dist"] / 10^3), log="xy", xlab="Paralogs", ylab="Sampled gene pairs", main="QQ-Plot\n(log scale)")
    abline(0,1, col="red")    

    # sampled close by dist and enhancer
    hist(abs(closePairs[,"dist"] / 10^3), 20, col=COL[1],
    main="Close cis paralogs", xlab="Distance (kb)")
    hist(abs(sampEhClosePairsCombined[,"dist"] / 10^3), 20, col=COL[2],
    main="Sampled by dist\n and enhancer", xlab="Distance (kb)")   
    qqplot(abs(closePairs[,"dist"] / 10^3), abs(sampEhClosePairsCombined[,"dist"] / 10^3), xlab="Paralogs", ylab="Sampled gene pairs", main="QQ-Plot")
    abline(0,1, col="red")    
    qqplot(abs(closePairs[,"dist"] / 10^3), abs(sampEhClosePairsCombined[,"dist"] / 10^3), log="xy", xlab="Paralogs", ylab="Sampled by dist and enhancer", main="QQ-Plot\n(log scale)")
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

    # number of linked enhancers
    hist(paraLinkedEnhancer[paraLinkedEnhancer<=50], 50, col=COL[1],
    main="Close cis paralogs", xlab="Enhancers")
    hist(sampEhCloseLinkedEnhancer[sampEhCloseLinkedEnhancer<=50], 50, col=COL[2],
    main="Sampled by dist\n and enhancer", xlab="Enhancers")    

    qqplot(paraLinkedEnhancer, sampEhCloseLinkedEnhancer, xlab="Enhancers in paralog genes", ylab="Enhancers in sampled genes",  main="QQ-Plot")
    abline(0,1, col="red")
    qqplot(paraLinkedEnhancer, sampEhCloseLinkedEnhancer, xlab="Enhancers in paralog genes", log="xy", ylab="Enhancers in sampled genes",  main="QQ-Plot\n(log scale)")
    abline(0,1, col="red")

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
# fraction of paralog pairs on same chromosom
#-----------------------------------------------------------------------

# Plot the percent of paralog pairs on the same chromosome
pdf(paste0(outPrefix, ".paralogPairs_on_same_chrom.barplot.pdf"), width=3.5)

    paraPercent = 100 * nrow(allCisPairs) / nrow(paralogPairsUniq)
    randPercent = sapply(randPairsInCis, nrow) * 100 / nrow(randPairs[[1]])
    height=c(paraPercent, mean(randPercent))

    par(cex=1.5, lwd=2)
    bp = my.barplot(height, 
        names=paste0(c("Paralogs\n (n=", "Random pairs\n (n="), c(nrow(paralogPairsUniq), paste0(N_RAND, "x",nrow(paralogPairsUniq))), c(")", ")")), 
        addValues=TRUE, col=COL_RAND,
        main="Gene pairs on\n same chromosome", 
        ylab="%")
    
    error.bar(bp,height, c(NA, sd(randPercent)))
        
dev.off()

#---------------------------------------------------------------
# fraction of paralog pairs with same strand on same chrom
#---------------------------------------------------------------
paraPercent = sapply(list(paralogPairsUniq, distalPairs, closePairs), function(gP) percentTrue(gP[,"sameStrand"]))
randPercent = lapply(list(randPairs, sampDistalPairs,sampClosePairs), function(gpl) {
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
# Paralog pairs across chromosomes
#-----------------------------------------------------------------------
# get number of pairs across each pair of chromosomes
interChromPairsMatrix = interChromPairMatrix(paralogPairsUniq, tssGR)
diag(interChromPairsMatrix) = NA

pdf(paste0(outPrefix, ".paralogPairs_interchrom_counts.heatmap.pdf"))
    my.heatmap.2(interChromPairsMatrix, Rowv=FALSE, Colv=FALSE, 
        revC=FALSE, dendrogram="none", 
        notecex=1, notecol="black", trace="none",
        margins=c(2, 2), key.xlab="Paralog pairs", 
        main="Human paralog gene pairs\n across chromosomes", col=colorRampPalette(brewer.pal(9,"Blues")), na.color="darkgray",
        tracecol="darkgreen", keysize=1.5,key.title="Color Key",
        ylab="Chromosomes", xlab="Chromosomes", ClabSide=3, RlabSide=2
    )
dev.off()

# take enrichment over random gene pairs on these chromosomes
# make random pairwise chrom matrix
randChromPairMatrix = bplapply(randPairs, interChromPairMatrix, tssGR)
Sys.sleep(3) # hack to fix problems with bplapply on MOGON

allRandChromPairMatrix =  Reduce("+", randChromPairMatrix)
diag(allRandChromPairMatrix) = NA

# calculate the log fold enrichment over random pairs
logFoldMatrix = log2(
        ((interChromPairsMatrix + 1) / sum(interChromPairsMatrix+1, na.rm=TRUE)) / 
        ((allRandChromPairMatrix +1) / sum(allRandChromPairMatrix+1, na.rm=TRUE))
    )

pdf(paste0(outPrefix, ".paralogPairs_interchrom_log_obs_vs_rand.heatmap.pdf"))
    my.heatmap.2(logFoldMatrix, Rowv=FALSE, Colv=FALSE, dendrogram="none", 
        notecex=1, notecol="black", trace="none",
        margins=c(2, 2), key.xlab="log_2(obs/exp)", 
        main="Enrichment of paralog pairs\n over random pairs", col=colorRampPalette(rev(brewer.pal(11,"RdBu"))), na.color="darkgray", 
        tracecol="darkgreen", keysize=1.5,key.title="Color Key",
        ylab="Chromosomes", xlab="Chromosomes", ClabSide=3, RlabSide=2
    )
dev.off()

pdf(paste0(outPrefix, ".paralogPairs_interchrom_log_obs_vs_rand.heatmap_dendrogram.pdf"))
    my.heatmap.2(logFoldMatrix, revC=TRUE,
        notecex=1, notecol="black", trace="none",
        margins=c(4, 4), key.xlab="log_2(obs/exp)", 
        main="Enrichment of paralog pairs\n over random pairs", col=colorRampPalette(rev(brewer.pal(11,"RdBu"))), na.color="darkgray", 
        tracecol="darkgreen", keysize=1.5,key.title="Color Key",
        ylab="Chromosomes", xlab="Chromosomes"
    )
dev.off()

#-----------------------------------------------------------------------
# 3) linear distance between paralogs
#-----------------------------------------------------------------------

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

# plot linear distance distribution for all distance thresholds:
allSampCloseDist = sampClosePairsCombined[,"dist"]

for (D in DIST_TH) {
    pdf(paste0(outPrefix, ".paralogPairs_cis_within_", D/10^3, "kb.hist.pdf"))
        par(cex=1.3, lwd=2, mfrow=c(2,1))
        hist(abs(closePairs[abs(closePairs$dist) <= D, "dist"])/10^3, 50, col=COL[1], main=paste("Linear distance distribution between paralogs\n on same chromosome within", D/10^3, "kb"), xlab="Distance (kb)")

        hist(abs(allSampCloseDist[abs(allSampCloseDist) <= D])/10^3, 50, col=COL[2], main=paste("Linear distance distribution between sampled gene pairs\n on same chromosome within", D/10^3, "kb"), xlab="Distance (kb)")
    dev.off()
}

pdf(paste0(outPrefix, ".paralogPairs_on_same_chrom_distance.boxplot.pdf"))
    par(cex=1.3, lwd=2)
    dists = lapply(c(list(abs(closePairs$dist)), lapply(sampClosePairs, function(df) abs(df[,"dist"]))), function(x){x/10^6}) 
    boxplot(dists, border=c(COL[1], rep(COL[2], N_RAND)),
        names=c("Paralog\n Pairs", 1:N_RAND),
        ylab="Distance (Mb)",main="Linear distance distribution\n between gene pairs on same chromosome"
    )
    mtext("Random pairs", side=1, line=2, at=mean(1:N_RAND +1), cex=1.3)
dev.off()

pdf(paste0(outPrefix, ".paralogPairs_on_same_chrom_distance.log_boxplot.pdf"))
    par(cex=1.3, lwd=2)
    dists = lapply(c(list(log10(abs(closePairs$dist)+1)), lapply(sampClosePairs, function(df) log10(abs(df[,"dist"])+1))), function(x){x/10^6}) 
    boxplot(dists, border=c(COL[1], rep(COL[2], N_RAND)),
        names=c("Paralog\n Pairs", 1:N_RAND),
        ylab="Distance (log_10 Mb)",main="Linear distance distribution\n between gene pairs on same chromosome"
    )
    mtext("Random pairs", side=1, line=2, at=mean(1:N_RAND +1), cex=1.3)
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
        ylab="TAD size (log_10 bp)", col=COL_DOMAIN)
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

    pdf(paste0(outPrefix, ".paralogPairs.within_", LABEL, ".barplot.pdf"), width=3.5)

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
    
    #---------------------------------------------------------------
    # compare same strand
    #---------------------------------------------------------------
    # create contingency table and run Fisher test
    paraStrand <- closePairs[,"sameStrand"]
    sampStrand <- lapply(sampClosePairs, function(gP) gP[, "sameStrand"])
    
    paraStrandInTAD <- paraStrand[paraInTAD]
    #sampStrandInTAD <- unlist(sampStrand)[unlist(sampInTAD)]
    sampStrandInTAD <- lapply(sampClosePairs, function(gP) gP[gP[,LABEL], "sameStrand"])

    contab <- rbind(
        para=table(paraStrand),  
        rand=table(unlist(sampStrand))
    )
    fs.test <- fisher.test(contab)

    contab.inTAD <- rbind(
        para=table(paraStrandInTAD),  
        rand=table(unlist(sampStrandInTAD))
    ) 
    fs.test.inTAD <- fisher.test(contab.inTAD)

    # get percent values for barplot
    paraStrandPercent = percentTrue(paraStrand)
    sampStrandPercent = sapply(sampStrand, percentTrue) 
    paraStrandInTADPercent = percentTrue(paraStrandInTAD)
    sampStrandInTADPercent = sapply(sampStrandInTAD, percentTrue)
    
    heights = matrix(c(paraStrandPercent, mean(sampStrandPercent, na.rm=TRUE), paraStrandInTADPercent, mean(sampStrandInTADPercent, na.rm=TRUE)), 2)
    sdSamp = c(sd(sampStrandPercent, na.rm=TRUE), sd(sampStrandInTADPercent, na.rm=TRUE))
    
    pdf(paste0(outPrefix, ".paralogPairs.SameStrand_within_", LABEL, ".barplot.pdf")) #, width=3.5

        par(cex=1.5, lwd=2)
        yMax = 1.6*max(heights) + max(sdSamp, na.rm=TRUE)
        
        bp = my.barplot(heights, addValues=TRUE, yMax=yMax,
            names=c("Close Pairs", "Close pairs in same TAD"),
            main=paste("Pairs with same strand\n in", LABEL, "TAD"), ylab="Gene pairs with same strand [%]", col=COL, srt=0, adj=NULL)

        error.bar(bp[2,], heights[2,], sdSamp, lwd=2)
        
        add_pval_two_pairs(bp, values=1.1*heights, pvalues=c(fs.test$p.value, fs.test.inTAD$p.value))
        
        legend("top", c("Paralogs", "Sampled"), fill=COL)

    dev.off()
    
}

# Plot for TADs from all cell-types together
pdf(paste0(outPrefix, ".paralogPairs.within_allTADs.barplot.pdf"), width=14)

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
    barplot(height, names.arg=names,
        beside=TRUE, col=COL2, main="Number of enhancers\n linked to single genes",
        xlab="Number of enhancers", ylab="%")
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
    percentCounts(rbind.fill(sampEhClosePairs)[,"commonEnhancer"]), all=TRUE, by="count")    
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
sampHasShared = lapply(sampEhClosePairs, function(df){df$commonEnhancer >= 1})
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

#save.image(paste0(WORKIMAGE_FILE, ".before_Hi-C_analysis.Rdata"))
#-----------------------------------------------------------------------
# Hi-C analysis with data from Rao et al. 2014
#-----------------------------------------------------------------------
breaksExp <- c(0, 1, 10, Inf)
breaksExpCor <- c(-1, -.5, 0, .5, 1)

#~ breaksClose = seq(1, MAX_DIST, length.out=9) / 10^3
breaksClose = 10^seq(2,6,1) / 10^3

#~ avgExp <- apply(closePairs[,c("g1_exp_IMR90", "g1_exp_IMR90")], 1, mean)

plotDFcloseHiC = data.frame(
        group = rep(c("paralogs", "sampled"), c(nrow(closePairs), nrow(sampClosePairsCombined))),
        HiCRaw = c(closePairs[,"HiCRaw"], sampClosePairsCombined[,"HiCRaw"]),
        HiC = c(closePairs[,"HiC"], sampClosePairsCombined[,"HiC"]),
        HiCobsExp = c(closePairs[,"HiCobsExp"], sampClosePairsCombined[,"HiCobsExp"]),
        captureC_raw = c(closePairs[,"captureC_raw"], sampClosePairsCombined[,"captureC_raw"]),
        captureC_ObsExp = c(closePairs[,"captureC_ObsExp"], sampClosePairsCombined[,"captureC_ObsExp"]),
        inTAD = factor(c(closePairs[,"Rao_IMR90"], sampClosePairsCombined[,"Rao_IMR90"]), levels=c(TRUE, FALSE), labels=c("same TAD", "diff TAD")),
        sameStrand = factor(c(closePairs[,"sameStrand"], sampClosePairsCombined[,"sameStrand"]), levels=c(TRUE, FALSE), labels=c("same strand", "opposite strand")),
        dist=abs(c(closePairs[,"dist"], sampClosePairsCombined[,"dist"]))/10^3,
        distBin=as.factor(breaksClose[.bincode(abs(c(closePairs[,"dist"], sampClosePairsCombined[,"dist"]))/10^3, breaksClose)]),
        avgExp=c(apply(closePairs[,c("g1_exp_IMR90", "g2_exp_IMR90")], 1, mean), apply(sampClosePairsCombined[,c("g1_exp_IMR90", "g2_exp_IMR90")], 1, mean)),
        expBin=factor(breaksExp[.bincode(c(apply(closePairs[,c("g1_exp_IMR90", "g2_exp_IMR90")], 1, mean), apply(sampClosePairsCombined[,c("g1_exp_IMR90", "g2_exp_IMR90")], 1, mean))+1, breaksExp)], labels=c("no", "low", "high"))
    )
plotDFcloseHiC <- addNoZero(plotDFcloseHiC)


# Make data.frame for plotting with ggplot2
breaksDistal = 10^seq(6,9,.25) / 10^3 

plotDFdistalHiC = data.frame(
        group = c(rep("paralogs", nrow(distalPairs)), rep("sampled", nrow(sampDistalPairsCombined))),
        HiCRaw = c(distalPairs[,"HiCRaw"], sampDistalPairsCombined[,"HiCRaw"]),
        HiC = c(distalPairs[,"HiC"], sampDistalPairsCombined[,"HiC"]),
        HiCobsExp = c(distalPairs[,"HiCobsExp"], sampDistalPairsCombined[,"HiCobsExp"]),
        captureC_raw = c(distalPairs[,"captureC_raw"], sampDistalPairsCombined[,"captureC_raw"]),
        captureC_ObsExp = c(distalPairs[,"captureC_ObsExp"], sampDistalPairsCombined[,"captureC_ObsExp"]),
        dist=abs(c(distalPairs[,"dist"], sampDistalPairsCombined[,"dist"]))/10^3,
        distBin=as.factor(breaksDistal[.bincode(abs(c(distalPairs[,"dist"], sampDistalPairsCombined[,"dist"]))/10^3, breaksDistal)]),
        avgExp=c(apply(distalPairs[,c("g1_exp_IMR90", "g2_exp_IMR90")], 1, mean), apply(sampDistalPairsCombined[,c("g1_exp_IMR90", "g2_exp_IMR90")], 1, mean)),
        expBin=factor(breaksExp[.bincode(c(apply(distalPairs[,c("g1_exp_IMR90", "g2_exp_IMR90")], 1, mean), apply(sampDistalPairsCombined[,c("g1_exp_IMR90", "g2_exp_IMR90")], 1, mean))+1, breaksExp)], labels=c("no", "low", "high"))
    )
plotDFdistalHiC <- addNoZero(plotDFdistalHiC)

#-----------------------------------------------------------------------
# Make data.frame for plotting with ggplot2
breaksAll = 10^(0:9) / 10^3

cC <- Reduce(intersect, lapply(list(closePairs, distalPairs, sampClosePairsCombined, sampDistalPairsCombined), names))
paraPairs <- rbind(closePairs[,cC], distalPairs[,cC])
sampPairs <- rbind(sampClosePairsCombined[,cC], sampDistalPairsCombined[,cC])
aP <- rbind(paraPairs, sampPairs)

plotDFallHiC = data.frame(
        group = c(rep("paralogs", nrow(paraPairs)), rep("sampled", nrow(sampPairs))),
        HiCRaw = aP[,"HiCRaw"],
        HiC = aP[,"HiC"],
        HiCobsExp = aP[,"HiCobsExp"],
        captureC_raw = aP[,"captureC_raw"],
        captureC_ObsExp = aP[,"captureC_ObsExp"],
        dist=abs(aP[,"dist"])/10^3,
        distBin=as.factor(breaksAll[.bincode(abs(aP[,"dist"])/10^3, breaksAll)]),
        avgExp=c(apply(paraPairs[,c("g1_exp_IMR90", "g2_exp_IMR90")], 1, mean), apply(sampPairs[,c("g1_exp_IMR90", "g2_exp_IMR90")], 1, mean)),
        expBin=factor(breaksExp[.bincode(c(apply(paraPairs[,c("g1_exp_IMR90", "g2_exp_IMR90")], 1, mean), apply(sampPairs[,c("g1_exp_IMR90", "g2_exp_IMR90")], 1, mean))+1, breaksExp)], labels=c("no", "low", "high"))
    )
plotDFallHiC <- addNoZero(plotDFallHiC)

#-----------------------------------------------------------------------
# Make data.frame for plotting with ggplot2 for all cis pairs
breaksCis = 10^(0:9) / 10^3

cC <- intersect(names(allCisPairs), names(sampCisPairsCombined))
aP <- rbind(allCisPairs[,cC], sampCisPairsCombined[,cC])

plotDFcisHiC = data.frame(
        group = c(rep("paralogs", nrow(allCisPairs)), rep("sampled", nrow(sampCisPairsCombined))),
        HiCRaw = aP[,"HiCRaw"],
        HiC = aP[,"HiC"],
        HiCobsExp = aP[,"HiCobsExp"],
        captureC_raw = aP[,"captureC_raw"],
        captureC_ObsExp = aP[,"captureC_ObsExp"],
        inTAD = factor(aP[,"Rao_IMR90"], levels=c(TRUE, FALSE), labels=c("same TAD", "diff TAD")),
        dist=abs(aP[,"dist"])/10^3,
        distBin=as.factor(breaksCis[.bincode(abs(aP[,"dist"])/10^3, breaksCis)]),
        avgExp=apply(aP[,c("g1_exp_IMR90", "g2_exp_IMR90")], 1, mean),
        expBin=factor(breaksExp[.bincode(apply(aP[,c("g1_exp_IMR90", "g2_exp_IMR90")], 1, mean)+1, breaksExp)], labels=c("no", "low", "high"))
    )
plotDFcisHiC <- addNoZero(plotDFcisHiC)

#-----------------------------------------------------------------------
# Make data.frame for plotting with ggplot2 for all cis pairs
breaksCis = 10^(0:9) / 10^3
breaksTAD = c(0, 10, 1000, Inf)

cC <- intersect(names(allCisPairs), names(sampCis100PairsCombined))
aP <- rbind(allCisPairs[,cC], sampCis100PairsCombined[,cC])

plotDFcis100HiC = data.frame(
        group = c(rep("paralogs", nrow(allCisPairs)), rep("sampled", nrow(sampCis100PairsCombined))),
        replicate = c(rep(1, nrow(allCisPairs)), rep(1:N_RAND, each=nrow(allCisPairs))), 
        HiCRaw = aP[,"HiCRaw"],
        HiC = aP[,"HiC"],
        HiCobsExp = aP[,"HiCobsExp"],
        captureC_raw = aP[,"captureC_raw"],
        captureC_ObsExp = aP[,"captureC_ObsExp"],
        inTAD = factor(aP[,"Rao_IMR90"], levels=c(TRUE, FALSE), labels=c("same TAD", "not same TAD")),
        subTAD = aP[,"Rao_IMR90_subTAD"],
        dupAgeGroup = factor(aP[,"mmusculus_commonOrtholg"], c(TRUE, FALSE), c("Young", "Old")),
        age = aP[,"age"],
        HIPPIE = aP[,"HIPPIE"],
        PPI = factor(aP[,"HIPPIE"] >= HIPPIE_MEDIUM, levels=c(TRUE, FALSE), labels=c("PPI", "no PPI")),
        dist=abs(aP[,"dist"])/10^3,
        distBin=as.factor(breaksCis[.bincode(abs(aP[,"dist"])/10^3, breaksCis)]),
        distTadBin=factor(breaksTAD[.bincode(abs(aP[,"dist"])/10^3, breaksTAD)], levels=breaksTAD[1:3], labels=c("<10kb", "10-1000kb", ">1000kb")),
        avgExp=apply(aP[,c("g1_exp_IMR90", "g2_exp_IMR90")], 1, mean),
        expBin=factor(breaksExp[.bincode(apply(aP[,c("g1_exp_IMR90", "g2_exp_IMR90")], 1, mean)+1, breaksExp)], labels=c("no", "low", "high")),
        expCor=aP[,"GTEx_expCor"],
        expCorBin=factor(breaksExpCor[.bincode(aP[,"GTEx_expCor"], breaksExpCor)], labels=c("high neg", "low neg", "low pos", "high pos"))
    )
plotDFcis100HiC <- addNoZero(plotDFcis100HiC)

plotDFcis100HiC[,"inTADclose"] <- as.character(plotDFcis100HiC$distTadBin)
plotDFcis100HiC[plotDFcis100HiC$distTadBin == "10-1000kb" ,"inTADclose"] <- as.character(plotDFcis100HiC[plotDFcis100HiC$distTadBin == "10-1000kb" ,"inTAD"])
plotDFcis100HiC$inTADclose <- factor(plotDFcis100HiC$inTADclose, levels=c("<10kb", "same TAD", "not same TAD", ">1000kb"))

# add three category sub TAD structure:
sub3TAD <- as.character(plotDFcis100HiC[,"subTAD"])
sub3TAD[sub3TAD == "no TAD" | sub3TAD == "diff TAD"] <- "not same TAD"
sub3TAD[sub3TAD == "diff sub TAD"] <- "same TAD"
plotDFcis100HiC[,"sub3TAD"] <- factor(sub3TAD, levels=c("not same TAD", "same TAD", "same sub TAD"))


#-----------------------------------------------------------------------
# Make Enhancer data.frame for plotting with ggplot2 for all cis pairs
cC <- intersect(names(closePairs), names(sampEhClosePairsCombined))
aP <- rbind(closePairs[,cC], sampEhClosePairsCombined[,cC])

plotDFeh = data.frame(
        group = c(rep("paralogs", nrow(closePairs)), rep("sampled", nrow(sampEhClosePairsCombined))),
        replicate = c(rep(1, nrow(closePairs)), rep(1:N_RAND, each=nrow(closePairs))), 
        inTAD = factor(aP[,"Rao_IMR90"], levels=c(TRUE, FALSE), labels=c("same TAD", "not same TAD")),
        subTAD = aP[,"Rao_IMR90_subTAD"],
        eh = aP[,"commonEnhancer"] > 0,
        commonEnhancer = aP[,"commonEnhancer"],
        dupAgeGroup = factor(aP[,"mmusculus_commonOrtholg"], c(TRUE, FALSE), c("Young", "Old")),
        age = aP[,"age"],
        HIPPIE = aP[,"HIPPIE"],
        PPI = factor(aP[,"HIPPIE"] >= HIPPIE_MEDIUM, levels=c(TRUE, FALSE), labels=c("PPI", "no PPI")),
        dist=abs(aP[,"dist"])/10^3,
        distBin=as.factor(breaksCis[.bincode(abs(aP[,"dist"])/10^3, breaksCis)]),
        distTadBin=factor(breaksTAD[.bincode(abs(aP[,"dist"])/10^3, breaksTAD)], levels=breaksTAD[1:3], labels=c("<10kb", "10-1000kb", ">1000kb")),
        avgExp=apply(aP[,c("g1_exp_IMR90", "g2_exp_IMR90")], 1, mean),
        expBin=factor(breaksExp[.bincode(apply(aP[,c("g1_exp_IMR90", "g2_exp_IMR90")], 1, mean)+1, breaksExp)], labels=c("no", "low", "high")),
        expCor=aP[,"GTEx_expCor"],
        expCorBin=factor(breaksExpCor[.bincode(aP[,"GTEx_expCor"], breaksExpCor)], labels=c("high neg", "low neg", "low pos", "high pos"))
    )

# add three category sub TAD structure:
sub3TAD <- as.character(plotDFeh[,"subTAD"])
sub3TAD[sub3TAD == "no TAD" | sub3TAD == "diff TAD"] <- "not same TAD"
sub3TAD[sub3TAD == "diff sub TAD"] <- "same TAD"
plotDFeh[,"sub3TAD"] <- factor(sub3TAD, levels=c("not same TAD", "same TAD", "same sub TAD"))

#-----------------------------------------------------------------------
# barplot percent of pairs by group, inTAD, distTadBin

# count gene pairs by combinations
dc <- ddply(plotDFcis100HiC, .(group, replicate, distTadBin, inTAD), summarize, count=length(inTAD), .drop=FALSE)
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
ggsave(p, file=paste0(outPrefix, ".cis100_pairs.gene_pairs.by_group_in_TAD_distTadBin.barplot.pdf"), w=7, h=7)

p <- ggplot(freqDF, aes(x=inTAD, y=avgPercent, fill=group)) +
    geom_bar(stat="identity", position="dodge", colour="black") + 
    facet_grid(.~distTadBin) + 
    geom_errorbar(aes(ymax = avgPercent + sdPercent , ymin=avgPercent - sdPercent), position=position_dodge(width=0.9), width=.25) +
    geom_text(aes(label=round(avgPercent), y= avgPercent), position=position_dodge(width=1), vjust=1.25, size=5) +
    scale_fill_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + labs(y="Gene pairs [n]", x="")
ggsave(p, file=paste0(outPrefix, ".cis100_pairs.percent_gene_pairs.by_group_in_TAD_distTadBin.barplot.pdf"), w=7, h=7)

#-----------------------------------------------------------------------
# percent in TAD by group and TADsource for [10-1000kb]

n <- nrow(plotDFcis100HiC)
# build data frame with all tadSources
plotDFcis100HiCtadSource <- plotDFcis100HiC[rep(1:n, length(allTADs)),]
plotDFcis100HiCtadSource[,"inTAD"] <- factor(sapply(names(allTADs), function(tad) aP[,tad]), levels=c(TRUE, FALSE), labels=c("same TAD", "not same TAD"))
plotDFcis100HiCtadSource[,"TADsource"] <- rep(names(allTADs), each=n)

subDF <- plotDFcis100HiCtadSource[plotDFcis100HiCtadSource$distTadBin =="10-1000kb",]

# count gene pairs by combinations
dc <- ddply(subDF, .(group, replicate, TADsource, inTAD), summarize, count=length(inTAD), .drop=FALSE)
# flag get existing replicate/group combination
exitingComb <- dc$replicate <= 1 | dc$group == "sampled"
# sum up counts for sub-combination to get percentage values
summedCounts <- rep(ddply(dc, .(group, replicate, TADsource), summarize, s=sum(count))$s, each=length(levels(dc$inTAD)))
# add percentage values and counts by marking non-existing combinations with NA
freqRepDF <- cbind(dc, cbind(summedCounts, existingCount=ifelse(exitingComb, dc$count, NA), percent=dc$count * 100 / summedCounts))

# combine replicates by taking average and sd of counts and percentages
freqDF <- ddply(freqRepDF, .(group, TADsource,inTAD), summarize, avgCount=mean(existingCount, na.rm=TRUE), sdCount=sd(existingCount, na.rm=TRUE), avgPercent=mean(percent, na.rm=TRUE), sdPercent=sd(percent, na.rm=TRUE))

# calculate p-values
pvalDF <- ddply(subDF, .(TADsource), summarize, p=fisher.test(group, inTAD)$p.value)

# combine p-values with plotting coordinates (max of values per group)
pvalDF <- merge(pvalDF, ddply(freqDF[freqDF$inTAD=="same TAD",], .(TADsource), summarize, avgCount=max(avgCount)+50, avgPercent=max(avgPercent)+5))
pvalDF[,"group"] <- "paralogs"

p <- ggplot(freqDF[freqDF$inTAD=="same TAD",], aes(x=TADsource, y=avgCount, fill=group)) +
    geom_bar(stat="identity", position="dodge", colour="black") + 
    geom_errorbar(aes(ymax = avgCount + sdCount , ymin=avgCount - sdCount), position=position_dodge(width=0.9), width=.25) +
    geom_text(aes(label=round(avgCount), y= avgCount), position=position_dodge(width=1), vjust=1.25, size=5) +
#~     geom_text(aes(label=paste0(round(avgCount), "\n(",round(avgPercent), "%)"), y= avgCount), position=position_dodge(width=1), vjust=1.25, size=5) +
    geom_text(aes(label=paste0("p=", signif(p,2))), data=pvalDF, size=5) + 
    scale_fill_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + labs(y="Gene pairs in same TAD", x="")
ggsave(p, file=paste0(outPrefix, ".cis100_pairs.gene_pairs.by_group_in_TAD_TADsource.barplot.pdf"), w=14, h=7)

p <- ggplot(freqDF[freqDF$inTAD=="same TAD",], aes(x=TADsource, y=avgPercent, fill=group)) +
    geom_bar(stat="identity", position="dodge", colour="black") + 
    geom_errorbar(aes(ymax = avgPercent + sdPercent , ymin=avgPercent - sdPercent), position=position_dodge(width=0.9), width=.25) +
    geom_text(aes(label=round(avgPercent), y= avgPercent), position=position_dodge(width=1), vjust=1.25, size=5) +
    geom_text(aes(label=paste0("p=", signif(p,2))), data=pvalDF, size=5) + 
    scale_fill_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + labs(y="Gene pairs in same TAD [%]", x="")
ggsave(p, file=paste0(outPrefix, ".cis100_pairs.percent_gene_pairs.by_group_in_TAD_TADsource.barplot.pdf"), w=14, h=7)


#------------------------------------------------------------------------
# HiC by group and inTADclose

# iterate over all contact measurements (in HiCcolumns):
#~ for (HiCcol in HiCcolumns){
for (HiCcol in c("HiC", "captureC_raw", "HiCNoZero", "captureC_rawNoZero")){
    
    message(paste("INFO: Plotting:", HiCcol))

    nDF <- ddply(plotDFcis100HiC, .(inTADclose, group), summarize, n=sum(!is.na(get(HiCcol))))
    
    # calculate p-values
    pvalDF <- ddply(plotDFcis100HiC, .(inTADclose), summarize, p=wilcox.test(as.formula(paste(HiCcol, "~ group")))$p.value)
#~     avgDF <- ddply(plotDFcis100HiC, .(inTADclose, group), summarize, mean=mean(HiCcol, na.rm=TRUE), median=median(HiCcol, na.rm=TRUE))
    p = ggplot(plotDFcis100HiC, aes_string(x="group", y=HiCcol)) + 
        geom_boxplot(aes(colour = group), lwd=1.5) + scale_y_log10() +
        facet_grid(.~inTADclose) + 
        scale_color_manual(values=COL, guide_legend(title = "")) +
        theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="none") + 
        guides(fill=guide_legend(title="")) +
        labs(x="") + 
        geom_text(aes(label=paste0("n=",n),  y=0.25), data=nDF) +
        geom_text(aes(label=paste0("p=", signif(p,2)), x=1.5, y=250), data=pvalDF, size=5)   
    ggsave(p, file=paste0(outPrefix, ".cis100_pairs.", HiCcol, ".by_inTADclose.boxplot.pdf"), w=7, h=7)
}

#-----------------------------------------------------------------------
# number of gene pairs by group and subTAD combination

subDF <- plotDFcis100HiC[plotDFcis100HiC$distTadBin == "10-1000kb",]

dc <- ddply(subDF, .(group, subTAD, replicate), summarize, count=length(subTAD))
summedCounts <- rep(ddply(dc, .(group, replicate), summarize, s=sum(count))$s, each=length(levels(dc$subTAD)))
freqRepDF <- cbind(dc, cbind(summedCounts, percent=dc$count * 100 / summedCounts))
freqDF <- ddply(freqRepDF, .(group, subTAD), summarize, avgCount=mean(count, na.rm=TRUE), sdCount=sd(count, na.rm=TRUE), avgPercent=mean(percent), sd=sd(percent))

# calculate p-values
pvals <- sapply(levels(subDF$subTAD), function(sT) fisher.test(table(subDF$subTAD == sT, subDF$group))$p.value) 
pvalDF <- data.frame(
    subTAD=levels(subDF$subTAD),
    p = pvals,
    ypos = ddply(freqDF, .(subTAD), summarize, y=max(avgCount))$y + 50,
    group=NA,
    dupAgeGroup=NA,
    age=NA
)


p <- ggplot(freqDF, aes(x=group, y=avgCount, fill=group)) +
    geom_bar(stat="identity", position="dodge", colour="black") + 
    geom_errorbar(aes(ymax = avgCount + sdCount , ymin=avgCount - sdCount), position=
    position_dodge(width=0.9), width=.25) +
    facet_grid(.~subTAD) +
    scale_fill_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_blank(), legend.position="bottom") + labs(y="Gene pairs", x="") + theme(strip.text.x = element_text(size=4, angle=90)) +
    geom_text(aes(label=avgCount, y= avgCount), position=position_dodge(width=1), vjust=1.25, size=5) +
    geom_text(aes(label=paste0("p=",signif(p,2)), y=ypos, x=1.5), size=5, data=pvalDF)


g <- addPictureLabels(p, subTADfigPaths)
pdf(paste0(outPrefix, ".cis100_pairs.10_1000kb.subTAD.byGroup.barplot.pdf"))
    grid.draw(g)
dev.off()

# number of gene pairs by dupAgeGroup, group and subTAD combination
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
    facet_grid(.~subTAD) +
    scale_fill_manual(values=c(COL_AGE, COL[2]), guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + labs(y="Gene pairs", x="") + theme(strip.text.x = element_text(size=4, angle=90)) +
    geom_text(aes(label=avgCount, y= avgCount), vjust=1.25, size=5) +
    geom_text(aes(label=paste0("p=",signif(p,2)), y=ypos, x=1.5), size=5, data=pvalDF)

g <- addPictureLabels(p, subTADfigPaths)

pdf(paste0(outPrefix, ".cis100_pairs.10_1000kb.subTAD.by_dupAgeGroup_and_group.barplot.pdf"))
    grid.draw(g)
dev.off()

#---------------------
# number of gene pairs by age, group and subTAD combination

# make dupAge for sampled pairs as NA
subDF[subDF$group == "sampled", "age"] <- NA

dcAge <- ddply(subDF, .(group, age, replicate), function(d) {data.frame(
    table(d$subTAD)
    )})
names(dcAge)[c(4,5)] <- c("subTAD", "count")
freqAgeDF <- ddply(dcAge, .(group, age, subTAD), summarize, avgCount=mean(count, na.rm=TRUE), sdCount=sd(count, na.rm=TRUE))
# convert to factor that do not exclude NAs for plotting
freqAgeDF$age <- factor(freqAgeDF$age, exclude = NULL)

p <- ggplot(freqAgeDF, aes(x=group, y=avgCount, fill=age)) +
    geom_bar(stat="identity", colour="black") + 
    geom_errorbar(aes(ymax = avgCount + sdCount , ymin=avgCount - sdCount), width=.25) +
    facet_grid(.~subTAD) +
    scale_fill_manual(values=COL_AGE_LEVELS, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1)) + labs(y="Gene pairs", x="") + theme(strip.text.x = element_text(size=4, angle=90)) +
#~     geom_text(aes(label=avgCount, y= avgCount), vjust=1.25, size=5) +
    geom_text(aes(label=paste0("p=",signif(p,2)), y=ypos, x=1.5), size=5, data=pvalDF)

g <- addPictureLabels(p, subTADfigPaths)

pdf(paste0(outPrefix, ".cis100_pairs.10_1000kb.subTAD.by_age_and_group.barplot.pdf"))
    grid.draw(g)
dev.off()

#-----------------------------------------------------------------------
# number of gene pairs with shared enhancer by group and subTAD combination
#-----------------------------------------------------------------------

# get subset for the size 10-1000kb
subDF <- plotDFeh[plotDFeh$distTadBin == "10-1000kb",]

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
pdf(paste0(outPrefix, ".cis100_pairs.10_1000kb.eh.by_subTAD_and_group.barplot.pdf"))
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
pdf(paste0(outPrefix, ".cis100_pairs.10_1000kb.ehPercent.by_subTAD_and_group.barplot.pdf"))
    grid.draw(g)
dev.off()


#---------------------
# subTAD by group
dc <- ddply(plotDFcis100HiC, .(group, subTAD, replicate), summarize, count=length(subTAD))
summedCounts <- rep(ddply(dc, .(group, replicate), summarize, s=sum(count))$s, each=length(levels(dc$subTAD)))
freqRepDF <- cbind(dc, cbind(summedCounts, percent=dc$count * 100 / summedCounts))
freqDF <- ddply(freqRepDF, .(group, subTAD), summarize, avgPercent=mean(percent), sd=sd(percent))

p <- ggplot(freqDF, aes(x=subTAD, y=avgPercent, fill=group)) +
    geom_bar(stat="identity", position="dodge", colour="black") + 
    geom_errorbar(aes(ymax = avgPercent + sd , ymin=avgPercent - sd), position=position_dodge(width=0.9), width=.25) +
    geom_text(aes(label=signif(avgPercent, 2), y= avgPercent), position=position_dodge(width=1), vjust=1.25, size=5) +
    scale_fill_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + labs(y="Gene pairs [%]", x="")
ggsave(p, file=paste0(outPrefix, ".cis100_pairs.subTAD.byGroup.barplot.pdf"), w=3.5, h=7)


#-----------------------------------------------------------------------
# HiC by group and subTAD combination
subDF <- plotDFcis100HiC[plotDFcis100HiC$distTadBin == "10-1000kb",]


# iterate over all contact measurements (in HiCcolumns):
#~ for (HiCcol in HiCcolumns){
for (HiCcol in c("HiC", "captureC_raw", "HiCNoZero", "captureC_rawNoZero")){
    
    message(paste("INFO: Plotting:", HiCcol))

    nDF <- ddply(subDF, .(subTAD, group), summarize, n=sum(!is.na(get(HiCcol))))
    
    # calculate p-values
    pvalDF <- ddply(subDF, .(subTAD), summarize, p=wilcox.test(as.formula(paste(HiCcol, "~ group")))$p.value)
#~     avgDF <- ddply(plotDFcis100HiC, .(inTADclose, group), summarize, mean=mean(HiCcol, na.rm=TRUE), median=median(HiCcol, na.rm=TRUE))
    p = ggplot(subDF, aes_string(x="group", y=HiCcol)) + 
        geom_boxplot(aes(colour = group), lwd=1.5) + scale_y_log10() +
        facet_grid(.~subTAD) + 
        scale_color_manual(values=COL, guide_legend(title = "")) +
        theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="none") + 
        guides(fill=guide_legend(title="")) +
        labs(x="") + 
        geom_text(aes(label=paste0("n=",n),  y=0.25), data=nDF) +
        geom_text(aes(label=paste0("p=", signif(p,2)), x=1.5, y=250), data=pvalDF, size=5)   
    ggsave(p, file=paste0(outPrefix, ".cis100_pairs.10_1000kb.", HiCcol, ".by_subTAD.boxplot.pdf"), w=7, h=7)
}

#-----------------------------------------------------------------------
# SUBTAD WITH 3 GROUPS
#-----------------------------------------------------------------------


#-----------------------------------------------------------------------
# number of gene pairs by group and sub3TAD combination

subDF <- plotDFcis100HiC[plotDFcis100HiC$distTadBin == "10-1000kb",]

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
    dupAgeGroup=NA,
    age=NA
)


p <- ggplot(freqDF, aes(x=group, y=avgCount, fill=group)) +
    geom_bar(stat="identity", position="dodge", colour="black") + 
    geom_errorbar(aes(ymax = avgCount + sdCount , ymin=avgCount - sdCount), position=
    position_dodge(width=0.9), width=.25) +
    facet_grid(.~sub3TAD) +
    scale_fill_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_blank(), legend.position="bottom") + labs(y="Gene pairs", x="") + theme(strip.text.x = element_text(size=4, angle=90)) +
    geom_text(aes(label=avgCount, y= avgCount), position=position_dodge(width=1), vjust=1.25, size=5) +
    geom_text(aes(label=paste0("p=",signif(p,2)), y=ypos, x=1.5), size=5, data=pvalDF)


g <- addPictureLabels(p, sub3TADfigPaths)
pdf(paste0(outPrefix, ".cis100_pairs.10_1000kb.sub3TAD.byGroup.barplot.pdf"))
    grid.draw(g)
dev.off()

# number of gene pairs by dupAgeGroup, group and sub3TAD combination
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
    geom_bar(stat="identity", colour="black") + 
    geom_errorbar(aes(ymax = avgCount + sdCount , ymin=avgCount - sdCount), width=.25) +
    facet_grid(.~sub3TAD) +
    scale_fill_manual(values=c(COL_AGE, COL[2]), guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + labs(y="Gene pairs", x="") + theme(strip.text.x = element_text(size=4, angle=90)) +
    geom_text(aes(label=avgCount, y= avgCount), vjust=1.25, size=5) +
    geom_text(aes(label=paste0("p=",signif(p,2)), y=ypos, x=1.5), size=5, data=pvalDF)

g <- addPictureLabels(p, sub3TADfigPaths)

pdf(paste0(outPrefix, ".cis100_pairs.10_1000kb.sub3TAD.by_dupAgeGroup_and_group.barplot.pdf"))
    grid.draw(g)
dev.off()

#---------------------
# number of gene pairs by age, group and sub3TAD combination

# make dupAge for sampled pairs as NA
subDF[subDF$group == "sampled", "age"] <- NA

dcAge <- ddply(subDF, .(group, age, replicate), function(d) {data.frame(
    table(d$sub3TAD)
    )})
names(dcAge)[c(4,5)] <- c("sub3TAD", "count")
freqAgeDF <- ddply(dcAge, .(group, age, sub3TAD), summarize, avgCount=mean(count, na.rm=TRUE), sdCount=sd(count, na.rm=TRUE))
# convert to factor that do not exclude NAs for plotting
freqAgeDF$age <- factor(freqAgeDF$age, exclude = NULL)

p <- ggplot(freqAgeDF, aes(x=group, y=avgCount, fill=age)) +
    geom_bar(stat="identity", colour="black") + 
    geom_errorbar(aes(ymax = avgCount + sdCount , ymin=avgCount - sdCount), width=.25) +
    facet_grid(.~sub3TAD) +
    scale_fill_manual(values=COL_AGE_LEVELS, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1)) + labs(y="Gene pairs", x="") + theme(strip.text.x = element_text(size=4, angle=90)) +
#~     geom_text(aes(label=avgCount, y= avgCount), vjust=1.25, size=5) +
    geom_text(aes(label=paste0("p=",signif(p,2)), y=ypos, x=1.5), size=5, data=pvalDF)

g <- addPictureLabels(p, sub3TADfigPaths)

pdf(paste0(outPrefix, ".cis100_pairs.10_1000kb.sub3TAD.by_age_and_group.barplot.pdf"))
    grid.draw(g)
dev.off()

#-----------------------------------------------------------------------
# number of gene pairs with shared enhancer by group and sub3TAD combination
#-----------------------------------------------------------------------

# get subset for the size 10-1000kb
subDF <- plotDFeh[plotDFeh$distTadBin == "10-1000kb",]

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
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_blank(), legend.position="bottom") + labs(y="Gene pairs with shared enhancer", x="") + theme(strip.text.x = element_text(size=4, angle=90)) +
    geom_text(aes(label=avgCount, y= avgCount), vjust=1.25, size=5) +
    geom_text(aes(label=paste0("p=",signif(p,2)), y=ypos, x=1.5), size=5, data=pvalDF)

g <- addPictureLabels(p, sub3TADfigPaths)
pdf(paste0(outPrefix, ".cis100_pairs.10_1000kb.eh.by_sub3TAD_and_group.barplot.pdf"))
    grid.draw(g)
dev.off()

# same with percent values
p <- ggplot(freqDF, aes(x=group, y=avgPercent, fill=group)) +
    geom_bar(stat="identity", colour="black") + 
    geom_errorbar(aes(ymax = avgPercent + sdPercent , ymin=avgPercent - sdPercent), width=.25) +
    facet_grid(.~sub3TAD) +
    scale_fill_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_blank(), legend.position="bottom") + labs(y="Percent of pairs with shared enhancer  [%]", x="") + theme(strip.text.x = element_text(size=4, angle=90)) +
    geom_text(aes(label=signif(avgPercent,3), y= avgPercent), vjust=1.25, size=5) +
    geom_text(aes(label=paste0("p=",signif(p,2)), y=yposPercent, x=1.5), size=5, data=pvalDF)

g <- addPictureLabels(p, sub3TADfigPaths)
pdf(paste0(outPrefix, ".cis100_pairs.10_1000kb.ehPercent.by_sub3TAD_and_group.barplot.pdf"))
    grid.draw(g)
dev.off()

#---------------------
# sub3TAD by group
dc <- ddply(plotDFcis100HiC, .(group, sub3TAD, replicate), summarize, count=length(sub3TAD))
summedCounts <- rep(ddply(dc, .(group, replicate), summarize, s=sum(count))$s, each=length(levels(dc$sub3TAD)))
freqRepDF <- cbind(dc, cbind(summedCounts, percent=dc$count * 100 / summedCounts))
freqDF <- ddply(freqRepDF, .(group, sub3TAD), summarize, avgPercent=mean(percent), sd=sd(percent))

p <- ggplot(freqDF, aes(x=sub3TAD, y=avgPercent, fill=group)) +
    geom_bar(stat="identity", position="dodge", colour="black") + 
    geom_errorbar(aes(ymax = avgPercent + sd , ymin=avgPercent - sd), position=position_dodge(width=0.9), width=.25) +
    geom_text(aes(label=signif(avgPercent, 2), y= avgPercent), position=position_dodge(width=1), vjust=1.25, size=5) +
    scale_fill_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + labs(y="Gene pairs [%]", x="")
ggsave(p, file=paste0(outPrefix, ".cis100_pairs.sub3TAD.byGroup.barplot.pdf"), w=3.5, h=7)


#-----------------------------------------------------------------------
# HiC by group and sub3TAD combination
subDF <- plotDFcis100HiC[plotDFcis100HiC$distTadBin == "10-1000kb",]


# iterate over all contact measurements (in HiCcolumns):
#~ for (HiCcol in HiCcolumns){
for (HiCcol in c("HiC", "captureC_raw", "HiCNoZero", "captureC_rawNoZero")){
    
    message(paste("INFO: Plotting:", HiCcol))

    nDF <- ddply(subDF, .(sub3TAD, group), summarize, n=sum(!is.na(get(HiCcol))))
    
    # calculate p-values
    pvalDF <- ddply(subDF, .(sub3TAD), summarize, p=wilcox.test(as.formula(paste(HiCcol, "~ group")))$p.value)
#~     avgDF <- ddply(plotDFcis100HiC, .(inTADclose, group), summarize, mean=mean(HiCcol, na.rm=TRUE), median=median(HiCcol, na.rm=TRUE))
    p = ggplot(subDF, aes_string(x="group", y=HiCcol)) + 
        geom_boxplot(aes(colour = group), lwd=1.5) + scale_y_log10() +
        facet_grid(.~sub3TAD) + 
        scale_color_manual(values=COL, guide_legend(title = "")) +
        theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="none") + 
        guides(fill=guide_legend(title="")) +
        labs(x="") + 
        geom_text(aes(label=paste0("n=",n),  y=0.25), data=nDF) +
        geom_text(aes(label=paste0("p=", signif(p,2)), x=1.5, y=250), data=pvalDF, size=5)   
    ggsave(p, file=paste0(outPrefix, ".cis100_pairs.10_1000kb.", HiCcol, ".by_sub3TAD.boxplot.pdf"), w=7, h=7)
}




#-------------------------------
# Denisty of average expression of pairs in IMR90
p <- ggplot(plotDFcloseHiC, aes(avgExp+1, ..density.., fill = group, color=group)) +
   geom_density(alpha=.5) + scale_x_log10() + 
   scale_color_manual(values=COL, guide_legend(title = "")) + scale_fill_manual(values=COL, guide_legend(title = "")) +
   theme_bw() + theme(text = element_text(size=20), legend.justification=c(1,1), legend.position=c(1,1))
ggsave(p, file=paste0(outPrefix, ".close_pairs.avgExp_IMR90.by_group.density.pdf"), w=3.5, h=3.5)    

p <- ggplot(plotDFcis100HiC, aes(avgExp+1, ..density.., fill = group, color=group)) +
   geom_density(alpha=.5) + scale_x_log10() + 
   scale_color_manual(values=COL, guide_legend(title = "")) + scale_fill_manual(values=COL, guide_legend(title = "")) +
   theme_bw() + theme(text = element_text(size=20), legend.justification=c(1,1), legend.position=c(1,1))
ggsave(p, file=paste0(outPrefix, ".cis100_pairs.avgExp_IMR90.by_group.density.pdf"), w=3.5, h=3.5)    

p <- dotplotWithDensityLogXY(revDF(plotDFcloseHiC[plotDFcis100HiC$dist > 1,]), "avgExp+1", "HiC", "group", COL, fit=TRUE, ALPHA=.2)
pdf(paste0(outPrefix, ".close_pairs.HiC.vs_avgExp.dotplot.pdf"))
    grid::grid.draw(p)
dev.off()

p <- dotplotWithDensityLogXY(revDF(plotDFcis100HiC[plotDFcis100HiC$dist > 1,]), "avgExp+1", "HiC", "group", COL, fit=TRUE, ALPHA=.2)
pdf(paste0(outPrefix, ".cis100_pairs.HiC.vs_avgExp.dotplot.pdf"))
    grid::grid.draw(p)
dev.off()


# expression correlation by TAD and TAD distance bins
nDF <- data.frame(table(plotDFcis100HiC[, c("group", "distTadBin", "inTAD")], useNA="ifany"))
p = ggplot(plotDFcis100HiC, aes(x=inTAD, y=expCor^2)) + 
    geom_boxplot(aes(colour = group), lwd=1.5) + 
    facet_grid(.~distTadBin) + 
    scale_color_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + 
    guides(fill=guide_legend(title="")) +
    labs(y="Expression Correlation [R^2]", x="Linear distance bin [kb]") + 
    geom_text(aes(label=paste0("n=",Freq), vjust=(as.numeric(group)-1)*2-1,  y=0, angle=90), data=nDF)
ggsave(p, file=paste0(outPrefix, ".cis100_pairs.expCor.by_TADdist_byTAD.boxplot.pdf"), w=7, h=7)


# expression correlation by subTAD and TAD distance bins
nDF <- data.frame(table(plotDFcis100HiC[, c("group", "distTadBin", "subTAD")], useNA="ifany"))
p = ggplot(plotDFcis100HiC, aes(x=subTAD, y=expCor^2)) + 
    geom_boxplot(aes(colour = group), lwd=1.5) + 
    facet_grid(.~distTadBin) + 
    scale_color_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + 
    guides(fill=guide_legend(title="")) +
    labs(y="Expression Correlation [R^2]", x="Linear distance bin [kb]") + 
    geom_text(aes(label=paste0("n=",Freq), vjust=(as.numeric(group)-1)*2-1,  y=0, angle=90), data=nDF)    
ggsave(p, file=paste0(outPrefix, ".cis100_pairs.expCor.by_TADdist_subTAD.boxplot.pdf"), w=7, h=7)

# expression correlation by subTAD 

#~ ggplot(nDF, aes(x=subTAD, y=Freq, fill=group)) +
#~     geom_bar(stat="identity", position="dodge") + scale_fill_manual(values=COL, guide_legend(title = ""))
nDF <- data.frame(table(plotDFcis100HiC[plotDFcis100HiC$distTadBin == 10, c("group", "subTAD")], useNA="ifany"))
p = ggplot(plotDFcis100HiC[plotDFcis100HiC$distTadBin == 10,], aes(x=group, y=expCor^2)) + 
    geom_boxplot(aes(colour = group), lwd=1.5) + 
    facet_grid(.~subTAD) + 
    scale_color_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="none") + 
    guides(fill=guide_legend(title="")) +
    labs(y="Expression Correlation [R^2]", x="") + 
    geom_text(aes(label=paste0("n=",Freq),  y=1.05), data=nDF)    
ggsave(p, file=paste0(outPrefix, ".cis100_pairs.10_1000kb.expCor.by_subTAD.boxplot.pdf"), w=7, h=7)


# HiC by subTAD 
p = ggplot(plotDFcis100HiC[plotDFcis100HiC$distTadBin == 10,], aes(x=group, y=HiC)) + 
    geom_boxplot(aes(colour = group), lwd=1.5) + scale_y_log10() +
    facet_grid(.~subTAD) + 
    scale_color_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="none") + 
    guides(fill=guide_legend(title="")) +
    labs(x="") + 
    geom_text(aes(label=paste0("n=",Freq),  y=110), data=nDF)    
ggsave(p, file=paste0(outPrefix, ".cis100_pairs.10_1000kb.HiC.by_subTAD.boxplot.pdf"), w=7, h=7)


# captureC_raw by subTAD 
p = ggplot(plotDFcis100HiC[plotDFcis100HiC$distTadBin == 10,], aes(x=group, y=captureC_raw)) + 
    geom_boxplot(aes(colour = group), lwd=1.5) + scale_y_log10() +
    facet_grid(.~subTAD) + 
    scale_color_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="none") + 
    guides(fill=guide_legend(title="")) +
    labs(x="") + 
    geom_text(aes(label=paste0("n=",Freq),  y=1000), data=nDF)    
ggsave(p, file=paste0(outPrefix, ".cis100_pairs.10_1000kb.captureC_raw.by_subTAD.boxplot.pdf"), w=7, h=7)
    
#-----------------------------------------------------------------------
# HIPPIE by TAD and distTadBin
nDF <- ddply(plotDFcis100HiC, .(group, inTAD, distTadBin), summarize, Freq=sum(!is.na(HIPPIE)))
p = ggplot(plotDFcis100HiC, aes(x=group, y=HIPPIE)) + 
    geom_boxplot(aes(colour = group), lwd=1.5) +
    facet_grid(.~distTadBin*inTAD) + 
    scale_color_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="none") + 
    guides(fill=guide_legend(title="")) +
    labs(x="") + 
    geom_text(aes(label=paste0("n=",Freq),  y=1), data=nDF)    
ggsave(p, file=paste0(outPrefix, ".cis100_pairs.HIPPIEscore.by_TAD_and_distTadBin.boxplot.pdf"), w=7, h=7)

# HIPPIE by TAD and groups [10kb - 1000kb]
nDF <- ddply(plotDFcis100HiC[plotDFcis100HiC$distTadBin == 10,], .(group, inTAD), summarize, Freq=sum(!is.na(HIPPIE)))
p = ggplot(plotDFcis100HiC[plotDFcis100HiC$distTadBin == 10,], aes(x=group, y=HIPPIE)) + 
    geom_boxplot(aes(colour = group), lwd=1.5) +
    facet_grid(.~inTAD) + 
    scale_color_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="none") + 
    guides(fill=guide_legend(title="")) +
    labs(x="") + 
    geom_text(aes(label=paste0("n=",Freq),  y=1), data=nDF)    
ggsave(p, file=paste0(outPrefix, ".cis100_pairs.10_1000kb.HIPPIEscore.by_group_and_TAD.boxplot.pdf"), w=3.5, h=7)


# HIPPIE by distTadBin and group
nDF <- ddply(plotDFcis100HiC, .(group, distTadBin), summarize, Freq=sum(!is.na(HIPPIE)))
p = ggplot(plotDFcis100HiC, aes(x=group, y=HIPPIE)) + 
    geom_boxplot(aes(colour = group), lwd=1.5) +
    facet_grid(.~distTadBin) + 
    scale_color_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="none") + 
    guides(fill=guide_legend(title="")) +
    labs(x="") + 
    geom_text(aes(label=paste0("n=",Freq),  y=1), data=nDF)    
ggsave(p, file=paste0(outPrefix, ".cis100_pairs.HIPPIEscore.by_distTadBin_and_group.boxplot.pdf"), w=3.5, h=7)


# PPI by TAD and group [10kb - 1000kb]
dc <- ddply(plotDFcis100HiC[plotDFcis100HiC$distTadBin == 10,], .(group, inTAD, replicate, PPI), summarize, count=length(PPI), .drop=FALSE)
summedCounts <- rep(ddply(dc, .(group, inTAD, replicate), summarize, s=sum(count))$s, each=length(levels(dc$PPI))+1)
freqRepDF <- cbind(dc, cbind(summedCounts, percent=dc$count * 100 / summedCounts))
freqDF <- ddply(freqRepDF, .(group, inTAD, PPI), summarize, avg=sum(count), avgPercent=mean(percent,na.rm=TRUE), sd=sd(percent,na.rm=TRUE))

p <- ggplot(freqDF, aes(x=PPI, y=avgPercent, fill=group)) +
    geom_bar(stat="identity", position="dodge", colour="black") + 
    geom_errorbar(aes(ymax = avgPercent + sd , ymin=avgPercent - sd), position=position_dodge(width=0.9), width=.25) +
    facet_grid(.~inTAD) +
    geom_text(aes(label=signif(avgPercent, 2), y= avgPercent), position=position_dodge(width=1), hjust=-0.25, size=5, angle=90) +
    scale_fill_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + labs(y="Gene pairs [%]", x="")
ggsave(p, file=paste0(outPrefix, ".cis100_pairs.10_1000kb.PPI.by_TAD_and_group.barplot.pdf"), w=3.5, h=7)

# PPI by TAD and group [10kb - 1000kb] only non NA
p <- ggplot(freqDF[freqDF$PPI == "PPI" & !is.na(freqDF$PPI),], aes(x=PPI, y=avgPercent, fill=group)) +
    geom_bar(stat="identity", position="dodge", colour="black") + 
    geom_errorbar(aes(ymax = avgPercent + sd , ymin=avgPercent - sd), position=position_dodge(width=0.9), width=.25) +
    facet_grid(.~inTAD) +
    geom_text(aes(label=signif(avgPercent, 2), y= avgPercent), position=position_dodge(width=1), vjust=1.25, size=5) +
    scale_fill_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + labs(y="Gene pairs with PPI  [%]", x="")
ggsave(p, file=paste0(outPrefix, ".cis100_pairs.10_1000kb.PPI_noNA.by_TAD_and_group.barplot.pdf"), w=3.5, h=7)

#-----------------------------------------------------------------------
# PPI by subTAD and group [10kb - 1000kb]
dc <- ddply(plotDFcis100HiC[plotDFcis100HiC$distTadBin == 10,], .(group, subTAD, replicate, PPI), summarize, count=length(PPI), .drop=FALSE)
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
ggsave(p, file=paste0(outPrefix, ".cis100_pairs.10_1000kb.PPI.by_subTAD_and_group.barplot.pdf"), w=7, h=7)

# PPI by subTAD and group [10kb - 1000kb] only non NA
p <- ggplot(freqDF[freqDF$PPI == "PPI" & !is.na(freqDF$PPI),], aes(x=PPI, y=avgPercent, fill=group)) +
    geom_bar(stat="identity", position="dodge", colour="black") + 
    geom_errorbar(aes(ymax = avgPercent + sd , ymin=avgPercent - sd), position=position_dodge(width=0.9), width=.25) +
    facet_grid(.~subTAD) +
    geom_text(aes(label=signif(avgPercent, 2), y= avgPercent), position=position_dodge(width=1), vjust=1.25, size=5) +
    scale_fill_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + labs(y="Gene pairs with PPI [%]", x="")
ggsave(p, file=paste0(outPrefix, ".cis100_pairs.10_1000kb.PPI_noNA.by_subTAD_and_group.barplot.pdf"), w=7, h=7)



#-----------------------------------------------------------------------
# PPI by subTAD and distTadBin
dc <- ddply(plotDFcis100HiC, .(group, subTAD, distTadBin, replicate, PPI), summarize, count=length(PPI), .drop=FALSE)
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
ggsave(p, file=paste0(outPrefix, ".cis100_pairs.PPI.by_subTAD_and_group.barplot.pdf"), w=14, h=7)

# PPI by subTAD and group [10kb - 1000kb] only non NA
p <- ggplot(freqDF[freqDF$PPI == "PPI" & !is.na(freqDF$PPI),], aes(x=PPI, y=avgPercent, fill=group)) +
    geom_bar(stat="identity", position="dodge", colour="black") + 
    geom_errorbar(aes(ymax = avgPercent + sd , ymin=avgPercent - sd), position=position_dodge(width=0.9), width=.25) +
    facet_grid(.~distTadBin*subTAD) +
    geom_text(aes(label=signif(avgPercent, 2), y= avgPercent), position=position_dodge(width=1), vjust=1.25, size=5) +
    scale_fill_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + labs(y="Gene pairs with PPI [%]", x="")
ggsave(p, file=paste0(outPrefix, ".cis100_pairs.PPI_noNA.by_subTAD_and_group.barplot.pdf"), w=7, h=7)


#-------------------------------

# iterate over all contact measurements (in HiCcolumns):
#~ for (HiCcol in HiCcolumns){
for (HiCcol in c("HiC", "captureC_raw")){
    
    message(paste("INFO: Plotting:", HiCcol))
    
    #------------------------------------------------------------------------
    # Close pairs:
    xlabels=paste0(c("Paralogs", "Sampled"), 
        "\n n=", applyToSubset(plotDFcloseHiC, function(v) sum(!is.na(v)), HiCcol, "group"), 
        "\n med=", signif(applyToSubset(plotDFcloseHiC, median, HiCcol, "group", na.rm=TRUE), 3),
        "\n avg=", signif(applyToSubset(plotDFcloseHiC, mean, HiCcol, "group", na.rm=TRUE), 3)
        )
    ws.test = wilcox.test(plotDFcloseHiC[,HiCcol] ~ plotDFcloseHiC[,"group"])
    
    p = ggplot(plotDFcloseHiC, aes_string(x="group", y=HiCcol, color="group")) + 
        geom_boxplot(lwd=1.5) + scale_y_log10() + annotation_logticks(sides="l") +
        scale_color_manual(values=COL, guide_legend(title = "")) +
        theme_bw() + theme(text = element_text(size=20), legend.position="none") + 
        guides(fill=guide_legend(title="")) +
        labs(y=HiCcol, x="", title=paste0("p = ", signif(ws.test$p.value, 3))) +
        scale_x_discrete(labels=xlabels )
    ggsave(p, file=paste0(outPrefix, ".close_pairs.", HiCcol, ".ggboxplot.pdf"), width=3.5)    

    # At least 100kb distance pairs:
    for (D in c(10, 100)){
        plotDFcloseHiCMinDist <- plotDFcloseHiC[plotDFcloseHiC$dist >= D,]
        xlabels=paste0(c("Paralogs", "Sampled"), 
            "\n n=", applyToSubset(plotDFcloseHiCMinDist, function(v) sum(!is.na(v)), HiCcol, "group"), 
            "\n med=", signif(applyToSubset(plotDFcloseHiCMinDist, median, HiCcol, "group", na.rm=TRUE), 3),
            "\n avg=", signif(applyToSubset(plotDFcloseHiCMinDist, mean, HiCcol, "group", na.rm=TRUE), 3)
            )
        ws.test = wilcox.test(plotDFcloseHiCMinDist[,HiCcol] ~ plotDFcloseHiCMinDist[,"group"])
    
        p = ggplot(plotDFcloseHiCMinDist, aes_string(x="group", y=HiCcol, color="group")) + 
            geom_boxplot(lwd=1.5)  + scale_y_log10() + annotation_logticks(sides="l") +
            scale_color_manual(values=COL, guide_legend(title = "")) +
            theme_bw() + theme(text = element_text(size=20), legend.position="none") + 
            guides(fill=guide_legend(title="")) +
            labs(y=HiCcol, x="", title=paste0("p = ", signif(ws.test$p.value, 3))) +
            scale_x_discrete(labels=xlabels )
        ggsave(p, file=paste0(outPrefix, ".close_pairs.min", D, "kb.", HiCcol, ".ggboxplot.pdf"), width=3.5)    
    }
    
    nDF <- data.frame(table(plotDFcloseHiC[!is.na(plotDFcloseHiC$distBin), c("group", "distBin")]))
    p = ggplot(plotDFcloseHiC[!is.na(plotDFcloseHiC$distBin), ], aes_string(x="group", y=HiCcol)) + 
        scale_y_log10() + annotation_logticks(sides="l") +
        geom_boxplot(aes(x=group, colour = group), lwd=1.5) + 
        #geom_jitter(aes(x=group, colour = group), alpha=0.1, width=.5) + 
        facet_grid(.~distBin) + 
        scale_color_manual(values=COL, guide_legend(title = "")) +
        theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="none") + 
        guides(fill=guide_legend(title="")) +
        labs(y=HiCcol, x="Linear distance bin [kb]") + 
        geom_text(aes(label=paste0("n=",Freq), y=max(plotDFcloseHiC[,HiCcol], na.rm=TRUE)), data=nDF)
        
    ggsave(p, file=paste0(outPrefix, ".close_pairs.", HiCcol, ".by_dist.boxplot.pdf"))

    p = ggplot(plotDFcloseHiC[!is.na(plotDFcloseHiC$distBin), ], aes_string(x="inTAD", y=HiCcol)) + 
        scale_y_log10() + annotation_logticks(sides="l") +
        geom_boxplot(aes(colour = group), lwd=1.5) + 
        facet_grid(.~distBin) + 
        scale_color_manual(values=COL, guide_legend(title = "")) +
        theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + 
        guides(fill=guide_legend(title="")) +
        labs(y=HiCcol, x="Linear distance bin [kb]")
    ggsave(p, file=paste0(outPrefix, ".close_pairs.", HiCcol, ".by_dist_byTAD.boxplot.pdf"), w=7, h=7)

    
    nDF <- data.frame(table(plotDFcloseHiC[!is.na(plotDFcloseHiC$distBin), c("group", "distBin", "sameStrand")]))
    p = ggplot(plotDFcloseHiC[!is.na(plotDFcloseHiC$distBin), ], aes_string(x="sameStrand", y=HiCcol)) + 
        scale_y_log10() + annotation_logticks(sides="l") +
        geom_boxplot(aes(colour = group), lwd=1.5) + 
        facet_grid(.~distBin) + 
        scale_color_manual(values=COL, guide_legend(title = "")) +
        theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + 
        guides(fill=guide_legend(title="")) +
        labs(y=HiCcol, x="Linear distance bin [kb]") +
        geom_text(aes(label=paste0("n=",Freq), vjust=(as.numeric(group)-1)*2-1,  y=max(plotDFcloseHiC[,HiCcol], na.rm=TRUE), hadjust=1, angle=90), data=nDF)

    ggsave(p, file=paste0(outPrefix, ".close_pairs.", HiCcol, ".by_dist_byStrand.boxplot.pdf"), w=7, h=7)
    
    nDF <- data.frame(table(plotDFcloseHiC[, c("group", "expBin")], useNA="ifany"))
    p = ggplot(plotDFcloseHiC, aes_string(x="group", y=HiCcol)) + 
        scale_y_log10() + annotation_logticks(sides="l") +
        geom_boxplot(aes(x=group, colour = group), lwd=1.5) + 
        #geom_jitter(aes(x=group, colour = group), alpha=0.1, width=.5) + 
        facet_grid(.~expBin) + 
        scale_color_manual(values=COL, guide_legend(title = "")) +
        theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="none") + 
        guides(fill=guide_legend(title="")) +
        labs(y=HiCcol, x="Average Expression") + 
        geom_text(aes(label=paste0("n=",Freq), y=max(plotDFcloseHiC[,HiCcol], na.rm=TRUE)), data=nDF)
    ggsave(p, file=paste0(outPrefix, ".close_pairs.", HiCcol, ".by_expBin.boxplot.pdf"))

    p <- dotplotWithDensityLogXY(revDF(plotDFcloseHiC), "avgExp+1", HiCcol, "group", COL, fit=TRUE, ALPHA=.2)
    pdf(paste0(outPrefix, ".close_pairs.", HiCcol, ".vs_avgExp.dotplot.pdf"))
        grid::grid.draw(p)
    dev.off()

    #------------------------------------------------------------------------
    # Distal pairs:
    xlabels=paste0(c("Paralogs", "Sampled"), 
        "\n n=", applyToSubset(plotDFdistalHiC, function(v) sum(!is.na(v)), HiCcol, "group"), 
        "\n med=", signif(applyToSubset(plotDFdistalHiC, median, HiCcol, "group", na.rm=TRUE), 3),
        "\n avg=", signif(applyToSubset(plotDFdistalHiC, mean, HiCcol, "group", na.rm=TRUE), 3)
        )
    ws.test = wilcox.test(plotDFdistalHiC[,HiCcol] ~ plotDFdistalHiC[,"group"])
    
    p = ggplot(plotDFdistalHiC, aes_string(x="group", y=HiCcol, color="group")) + 
        geom_boxplot(lwd=1.5) + scale_y_log10() + annotation_logticks(sides="l") +
        scale_color_manual(values=COL, guide_legend(title = "")) +
        theme_bw() + theme(text = element_text(size=20), legend.position="none") + 
        guides(fill=guide_legend(title="")) +
        labs(y=HiCcol, x="", title=paste0("p = ", signif(ws.test$p.value, 3))) +
        scale_x_discrete(labels=xlabels )
    ggsave(p, file=paste0(outPrefix, ".distal_pairs.", HiCcol, ".ggboxplot.pdf"), width=3.5)    

    p = ggplot(plotDFdistalHiC, aes_string(x=as.numeric("distBin"), y=HiCcol)) + 
        scale_y_log10() + annotation_logticks(sides="l") +
        geom_boxplot(aes(x=distBin, colour = group), lwd=1.5) + 
        scale_color_manual(values=COL, guide_legend(title = "")) +
        theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + 
        guides(fill=guide_legend(title="")) +
        labs(y=HiCcol, x="Linear distance bin [kb]")
    ggsave(p, file=paste0(outPrefix, ".distal_pairs.", HiCcol, ".by_dist.boxplot.pdf"))

    nDF <- data.frame(table(plotDFdistalHiC[, c("group", "expBin")], useNA="ifany"))
    p = ggplot(plotDFdistalHiC, aes_string(x="group", y=HiCcol)) + 
        scale_y_log10() + annotation_logticks(sides="l") +
        geom_boxplot(aes(x=group, colour = group), lwd=1.5) + 
        #geom_jitter(aes(x=group, colour = group), alpha=0.1, width=.5) + 
        facet_grid(.~expBin) + 
        scale_color_manual(values=COL, guide_legend(title = "")) +
        theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="none") + 
        guides(fill=guide_legend(title="")) +
        labs(y=HiCcol, x="Average Expression") + 
        geom_text(aes(label=paste0("n=",Freq), y=max(plotDFdistalHiC[,HiCcol], na.rm=TRUE)), data=nDF)
    ggsave(p, file=paste0(outPrefix, ".distal_pairs.", HiCcol, ".by_expBin.boxplot.pdf"))


    p <- dotplotWithDensityLogXY(revDF(plotDFdistalHiC), "avgExp+1", HiCcol, "group", COL, fit=TRUE, ALPHA=.2)
    pdf(paste0(outPrefix, ".distal_pairs.", HiCcol, ".vs_avgExp.dotplot.pdf"))
        grid::grid.draw(p)
    dev.off()
    
    #------------------------------------------------------------------------
    # All pairs:
    xlabels=paste0(c("Paralogs", "Sampled"), 
        "\n n=", applyToSubset(plotDFallHiC, function(v) sum(!is.na(v)), HiCcol, "group"), 
        "\n med=", signif(applyToSubset(plotDFallHiC, median, HiCcol, "group", na.rm=TRUE), 3),
        "\n avg=", signif(applyToSubset(plotDFallHiC, mean, HiCcol, "group", na.rm=TRUE), 3)
        )
    ws.test = wilcox.test(plotDFallHiC[,HiCcol] ~ plotDFallHiC[,"group"])
    
    p = ggplot(plotDFallHiC, aes_string(x="group", y=HiCcol, color="group")) + 
        geom_boxplot(lwd=1.5) + scale_y_log10() + annotation_logticks(sides="l") +
        scale_color_manual(values=COL, guide_legend(title = "")) +
        theme_bw() + theme(text = element_text(size=20), legend.position="none") + 
        guides(fill=guide_legend(title="")) +
        labs(y=HiCcol, x="", title=paste0("p = ", signif(ws.test$p.value, 3))) +
        scale_x_discrete(labels=xlabels )
    ggsave(p, file=paste0(outPrefix, ".all_pairs.", HiCcol, ".ggboxplot.pdf"), width=3.5)    
    
    p = ggplot(plotDFallHiC, aes_string(x="distBin", y=HiCcol)) + 
        scale_y_log10() + annotation_logticks(sides="l") +
        geom_boxplot(aes(x=distBin, colour = group), lwd=1.5) + 
        scale_color_manual(values=COL, guide_legend(title = "")) +
        theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + 
        guides(fill=guide_legend(title="")) +
        labs(y=HiCcol, x="Linear distance bin [kb]")
    ggsave(p, file=paste0(outPrefix, ".all_pairs.", HiCcol, ".by_dist.boxplot.pdf"))

    #------------------------------------------------------------------------
    # Cis pairs:
    xlabels=paste0(c("Paralogs", "Sampled"), 
        "\n n=", applyToSubset(plotDFcisHiC, function(v) sum(!is.na(v)), HiCcol, "group"), 
        "\n med=", signif(applyToSubset(plotDFcisHiC, median, HiCcol, "group", na.rm=TRUE), 3),
        "\n avg=", signif(applyToSubset(plotDFcisHiC, mean, HiCcol, "group", na.rm=TRUE), 3)
        )
    ws.test = wilcox.test(plotDFcisHiC[,HiCcol] ~ plotDFcisHiC[,"group"])
    
    p = ggplot(plotDFcisHiC, aes_string(x="group", y=HiCcol, color="group")) + 
        geom_boxplot(lwd=1.5) + scale_y_log10() + annotation_logticks(sides="l") +
        scale_color_manual(values=COL, guide_legend(title = "")) +
        theme_bw() + theme(text = element_text(size=20), legend.position="none") + 
        guides(fill=guide_legend(title="")) +
        labs(y=HiCcol, x="", title=paste0("p = ", signif(ws.test$p.value, 3))) +
        scale_x_discrete(labels=xlabels )
    ggsave(p, file=paste0(outPrefix, ".cis_pairs.", HiCcol, ".ggboxplot.pdf"), width=3.5)    

    # Minimum X kb distance pairs:
    for (D in c(10, 100)){
        plotDFcisHiCMinDist <- plotDFcisHiC[plotDFcisHiC$dist >= D,]
        xlabels=paste0(c("Paralogs", "Sampled"), 
            "\n n=", applyToSubset(plotDFcisHiCMinDist, function(v) sum(!is.na(v)), HiCcol, "group"), 
            "\n med=", signif(applyToSubset(plotDFcisHiCMinDist, median, HiCcol, "group", na.rm=TRUE), 3),
            "\n avg=", signif(applyToSubset(plotDFcisHiCMinDist, mean, HiCcol, "group", na.rm=TRUE), 3)
            )
        ws.test = wilcox.test(plotDFcisHiCMinDist[,HiCcol] ~ plotDFcisHiCMinDist[,"group"])
    
        p = ggplot(plotDFcisHiCMinDist, aes_string(x="group", y=HiCcol, color="group")) + 
            geom_boxplot(lwd=1.5)  + scale_y_log10() + annotation_logticks(sides="l") +
            scale_color_manual(values=COL, guide_legend(title = "")) +
            theme_bw() + theme(text = element_text(size=20), legend.position="none") + 
            guides(fill=guide_legend(title="")) +
            labs(y=HiCcol, x="", title=paste0("p = ", signif(ws.test$p.value, 3))) +
            scale_x_discrete(labels=xlabels )
        ggsave(p, file=paste0(outPrefix, ".cis_pairs.min", D, "kb.", HiCcol, ".ggboxplot.pdf"), width=3.5)    
    }
    
    p = ggplot(plotDFcisHiC, aes_string(x="distBin", y=HiCcol)) + 
        scale_y_log10() + annotation_logticks(sides="l") +
        geom_boxplot(aes(x=distBin, colour = group), lwd=1.5) + 
        scale_color_manual(values=COL, guide_legend(title = "")) +
        theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + 
        guides(fill=guide_legend(title="")) +
        labs(y=HiCcol, x="Linear distance bin [kb]")
    ggsave(p, file=paste0(outPrefix, ".cis_pairs.", HiCcol, ".by_dist.boxplot.pdf"))
    
    p = ggplot(plotDFcisHiC, aes_string(x="inTAD", y=HiCcol)) + 
        scale_y_log10() + annotation_logticks(sides="l") +
        geom_boxplot(aes(colour = group), lwd=1.5) + 
        facet_grid(.~distBin) + 
        scale_color_manual(values=COL, guide_legend(title = "")) +
        theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + 
        guides(fill=guide_legend(title="")) +
        labs(y=HiCcol, x="Linear distance bin [kb]")
    ggsave(p, file=paste0(outPrefix, ".cis_pairs.", HiCcol, ".by_dist_byTAD.boxplot.pdf"), w=14, h=7)

    #------------------------------------------------------------------------
    # Cis pairs with 100 bin sampling:
    xlabels=paste0(c("Paralogs", "Sampled"), 
        "\n n=", applyToSubset(plotDFcis100HiC, function(v) sum(!is.na(v)), HiCcol, "group"), 
        "\n med=", signif(applyToSubset(plotDFcis100HiC, median, HiCcol, "group", na.rm=TRUE), 3),
        "\n avg=", signif(applyToSubset(plotDFcis100HiC, mean, HiCcol, "group", na.rm=TRUE), 3)
        )
    ws.test = wilcox.test(plotDFcis100HiC[,HiCcol] ~ plotDFcis100HiC[,"group"])
    
    p = ggplot(plotDFcis100HiC, aes_string(x="group", y=HiCcol, color="group")) + 
        geom_boxplot(lwd=1.5) + scale_y_log10() + annotation_logticks(sides="l") +
        scale_color_manual(values=COL, guide_legend(title = "")) +
        theme_bw() + theme(text = element_text(size=20), legend.position="none") + 
        guides(fill=guide_legend(title="")) +
        labs(y=HiCcol, x="", title=paste0("p = ", signif(ws.test$p.value, 3))) +
        scale_x_discrete(labels=xlabels )
    ggsave(p, file=paste0(outPrefix, ".cis100_pairs.", HiCcol, ".ggboxplot.pdf"), width=3.5)    

    # Minimum X kb distance pairs:
    for (D in c(10, 100)){
        plotDFcis100HiCMinDist <- plotDFcis100HiC[plotDFcis100HiC$dist >= D,]
        xlabels=paste0(c("Paralogs", "Sampled"), 
            "\n n=", applyToSubset(plotDFcis100HiCMinDist, function(v) sum(!is.na(v)), HiCcol, "group"), 
            "\n med=", signif(applyToSubset(plotDFcis100HiCMinDist, median, HiCcol, "group", na.rm=TRUE), 3),
            "\n avg=", signif(applyToSubset(plotDFcis100HiCMinDist, mean, HiCcol, "group", na.rm=TRUE), 3)
            )
        ws.test = wilcox.test(plotDFcis100HiCMinDist[,HiCcol] ~ plotDFcis100HiCMinDist[,"group"])
    
        p = ggplot(plotDFcis100HiCMinDist, aes_string(x="group", y=HiCcol, color="group")) + 
            geom_boxplot(lwd=1.5)  + scale_y_log10() + annotation_logticks(sides="l") +
            scale_color_manual(values=COL, guide_legend(title = "")) +
            theme_bw() + theme(text = element_text(size=20), legend.position="none") + 
            guides(fill=guide_legend(title="")) +
            labs(y=HiCcol, x="", title=paste0("p = ", signif(ws.test$p.value, 3))) +
            scale_x_discrete(labels=xlabels )
        ggsave(p, file=paste0(outPrefix, ".cis100_pairs.min", D, "kb.", HiCcol, ".ggboxplot.pdf"), width=3.5)    
    }
    
    # average expression
    nDF <- data.frame(table(plotDFcis100HiC[, c("group", "expBin")], useNA="ifany"))
    p = ggplot(plotDFcis100HiC, aes_string(x="group", y=HiCcol)) + 
        scale_y_log10() + annotation_logticks(sides="l") +
        geom_boxplot(aes(x=group, colour = group), lwd=1.5) + 
        #geom_jitter(aes(x=group, colour = group), alpha=0.1, width=.5) + 
        facet_grid(.~expBin) + 
        scale_color_manual(values=COL, guide_legend(title = "")) +
        theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="none") + 
        guides(fill=guide_legend(title="")) +
        labs(y=HiCcol, x="Average Expression") + 
        geom_text(aes(label=paste0("n=",Freq), y=max(plotDFcis100HiC[,HiCcol], na.rm=TRUE)), data=nDF)
    ggsave(p, file=paste0(outPrefix, ".cis100_pairs.", HiCcol, ".by_expBin.boxplot.pdf"))

    p <- dotplotWithDensityLogXY(revDF(plotDFcis100HiC), "avgExp+1", HiCcol, "group", COL, fit=TRUE, ALPHA=.2)
    pdf(paste0(outPrefix, ".cis100_pairs.", HiCcol, ".vs_avgExp.dotplot.pdf"))
        grid::grid.draw(p)
    dev.off()

    # expression correlation
    nDF <- data.frame(table(plotDFcis100HiC[, c("group", "expCorBin")], useNA="ifany"))
    p = ggplot(plotDFcis100HiC, aes_string(x="group", y=HiCcol)) + 
        scale_y_log10() + annotation_logticks(sides="l") +
        geom_boxplot(aes(x=group, colour = group), lwd=1.5) + 
        facet_grid(.~expCorBin) + 
        scale_color_manual(values=COL, guide_legend(title = "")) +
        theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="none") + 
        guides(fill=guide_legend(title="")) +
        labs(y=HiCcol, x="Expression Correlation") + 
        geom_text(aes(label=paste0("n=",Freq), y=max(plotDFcis100HiC[,HiCcol], na.rm=TRUE)), data=nDF)
    ggsave(p, file=paste0(outPrefix, ".cis100_pairs.", HiCcol, ".by_expCorBin.boxplot.pdf"))

    p <- dotplotWithDensityLogXY(revDF(plotDFcis100HiC), "expCor", HiCcol, "group", COL, fit=TRUE, ALPHA=.2, xlog=FALSE, ylog=TRUE)
    pdf(paste0(outPrefix, ".cis100_pairs.", HiCcol, ".vs_expCor.dotplot.pdf"))
        grid::grid.draw(p)
    dev.off()
    
    # by TAD and TAD distance bins
    nDF <- data.frame(table(plotDFcis100HiC[, c("group", "distTadBin", "inTAD")], useNA="ifany"))
    p = ggplot(plotDFcis100HiC, aes_string(x="inTAD", y=HiCcol)) + 
        scale_y_log10() + annotation_logticks(sides="l") +
        geom_boxplot(aes(colour = group), lwd=1.5) + 
        facet_grid(.~distTadBin) + 
        scale_color_manual(values=COL, guide_legend(title = "")) +
        theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + 
        guides(fill=guide_legend(title="")) +
        labs(y=HiCcol, x="Linear distance bin [kb]") + 
        geom_text(aes(label=paste0("n=",Freq), vjust=(as.numeric(group)-1)*2-1,  y=10^-1, angle=90), data=nDF)
    ggsave(p, file=paste0(outPrefix, ".cis100_pairs.", HiCcol, ".by_TADdist_byTAD.boxplot.pdf"), w=7, h=7)

}


#------------------------------------------------------------------------
# plot Hi-C contacts versus genomic distance
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
p0a <- dotplotWithDensityLogXY(revDF(plotDFcloseHiC), "dist", "HiCRaw", "group", COL, fit=TRUE, ALPHA=.2)
p0b <- dotplotWithDensityLogXY(revDF(plotDFcloseHiC), "dist", "HiCRawNoZero", "group", COL, fit=TRUE, ALPHA=.2)
p1 <- dotplotWithDensityLogXY(revDF(plotDFcloseHiC), "dist", "HiC", "group", COL, fit=TRUE, ALPHA=.2)
p2 <- dotplotWithDensityLogXY(revDF(plotDFcloseHiC), "dist", "HiCobsExp", "group", COL, fit=TRUE, ALPHA=.2)
p3 <- dotplotWithDensityLogXY(revDF(plotDFcloseHiC), "dist", "captureC_raw", "group", COL, fit=TRUE, ALPHA=.2)
p4 <- dotplotWithDensityLogXY(revDF(plotDFcloseHiC), "dist", "captureC_ObsExp", "group", COL, fit=TRUE, ALPHA=.2)

pdf(paste0(outPrefix, ".close_pairs.Hi-C_vs_dist_all.ggboxplot.pdf"), w=12, h=18)
    do.call(grid.arrange, list(p0a, p0b, p1,p2,p3,p4)) 
dev.off()

#-----------------------------------------------------------------------
p0a <- dotplotWithDensityLogXY(revDF(plotDFdistalHiC), "dist", "HiCRaw", "group", COL, fit=TRUE, ALPHA=.2)
p0b <- dotplotWithDensityLogXY(revDF(plotDFdistalHiC), "dist", "HiCRawNoZero", "group", COL, fit=TRUE, ALPHA=.2)
p1 <- dotplotWithDensityLogXY(revDF(plotDFdistalHiC), "dist", "HiC", "group", COL, fit=TRUE, ALPHA=.2)
p2 <- dotplotWithDensityLogXY(revDF(plotDFdistalHiC), "dist", "HiCobsExp", "group", COL, fit=TRUE, ALPHA=.2)
p3 <- dotplotWithDensityLogXY(revDF(plotDFdistalHiC), "dist", "captureC_raw", "group", COL, fit=TRUE, ALPHA=.2)
p4 <- dotplotWithDensityLogXY(revDF(plotDFdistalHiC), "dist", "captureC_ObsExp", "group", COL, fit=TRUE, ALPHA=.2)
pdf(paste0(outPrefix, ".distal_pairs.Hi-C_vs_dist_all.ggboxplot.pdf"), w=12, h=18)
    do.call(grid.arrange, list(p0a, p0b, p1,p2,p3,p4)) 
dev.off()

#-----------------------------------------------------------------------
p0a <- dotplotWithDensityLogXY(revDF(plotDFallHiC), "dist", "HiCRaw", "group", COL, fit=TRUE, ALPHA=.2)
p0b <- dotplotWithDensityLogXY(revDF(plotDFallHiC), "dist", "HiCRawNoZero", "group", COL, fit=TRUE, ALPHA=.2)
p1 <- dotplotWithDensityLogXY(revDF(plotDFallHiC), "dist", "HiC", "group", COL, fit=TRUE, ALPHA=.2)
p2 <- dotplotWithDensityLogXY(revDF(plotDFallHiC), "dist", "HiCobsExp", "group", COL, fit=TRUE, ALPHA=.2)
p3 <- dotplotWithDensityLogXY(revDF(plotDFallHiC), "dist", "captureC_raw", "group", COL, fit=TRUE, ALPHA=.2)
p4 <- dotplotWithDensityLogXY(revDF(plotDFallHiC), "dist", "captureC_ObsExp", "group", COL, fit=TRUE, ALPHA=.2)

pdf(paste0(outPrefix, ".all_pairs.Hi-C_vs_dist_all.ggboxplot.pdf"), w=12, h=18)
    do.call(grid.arrange, list(p0a, p0b, p1,p2,p3,p4)) 
dev.off()

#-----------------------------------------------------------------------
p0a <- dotplotWithDensityLogXY(revDF(plotDFcisHiC), "dist", "HiCRaw", "group", COL, fit=TRUE, ALPHA=.2)
p0b <- dotplotWithDensityLogXY(revDF(plotDFcisHiC), "dist", "HiCRawNoZero", "group", COL, fit=TRUE, ALPHA=.2)
p1 <- dotplotWithDensityLogXY(revDF(plotDFcisHiC), "dist", "HiC", "group", COL, fit=TRUE, ALPHA=.2)
p2 <- dotplotWithDensityLogXY(revDF(plotDFcisHiC), "dist", "HiCobsExp", "group", COL, fit=TRUE, ALPHA=.2)
p3 <- dotplotWithDensityLogXY(revDF(plotDFcisHiC), "dist", "captureC_raw", "group", COL, fit=TRUE, ALPHA=.2)
p4 <- dotplotWithDensityLogXY(revDF(plotDFcisHiC), "dist", "captureC_ObsExp", "group", COL, fit=TRUE, ALPHA=.2)

pdf(paste0(outPrefix, ".cis_pairs.HiC.bin20.vs_dist.dotplot.pdf"))
    grid::grid.draw(p1)
dev.off()

pdf(paste0(outPrefix, ".cis_pairs.captureC_raw.bin20.vs_dist.dotplot.pdf"))
    grid::grid.draw(p3)
dev.off()

pdf(paste0(outPrefix, ".cis_pairs.Hi-C_vs_dist_all.ggboxplot.pdf"), w=12, h=18)
    do.call(grid.arrange, list(p0a, p0b, p1,p2,p3,p4)) 
dev.off()

#-----------------------------------------------------------------------
p0a <- dotplotWithDensityLogXY(revDF(plotDFcis100HiC[plotDFcis100HiC$dist > 1,]), "dist", "HiCRaw", "group", COL, fit=TRUE, ALPHA=.2)
p0b <- dotplotWithDensityLogXY(revDF(plotDFcis100HiC[plotDFcis100HiC$dist > 1,]), "dist", "HiCRawNoZero", "group", COL, fit=TRUE, ALPHA=.2)
p1 <- dotplotWithDensityLogXY(revDF(plotDFcis100HiC[plotDFcis100HiC$dist > 1,]), "dist", "HiC", "group", COL, fit=TRUE, ALPHA=.2)
p2 <- dotplotWithDensityLogXY(revDF(plotDFcis100HiC[plotDFcis100HiC$dist > 1,]), "dist", "HiCobsExp", "group", COL, fit=TRUE, ALPHA=.2)
p3 <- dotplotWithDensityLogXY(revDF(plotDFcis100HiC[plotDFcis100HiC$dist > 1,]), "dist", "captureC_raw", "group", COL, fit=TRUE, ALPHA=.2)
p4 <- dotplotWithDensityLogXY(revDF(plotDFcis100HiC[plotDFcis100HiC$dist > 1,]), "dist", "captureC_ObsExp", "group", COL, fit=TRUE, ALPHA=.2)

pdf(paste0(outPrefix, ".cis100_pairs.HiC.vs_dist.dotplot.pdf"))
    grid::grid.draw(p1)
dev.off()

pdf(paste0(outPrefix, ".cis100_pairs.captureC_raw.vs_dist.dotplot.pdf"))
    grid::grid.draw(p3)
dev.off()

pdf(paste0(outPrefix, ".cis100_pairs.Hi-C_vs_dist_all.ggboxplot.pdf"), w=12, h=18)
    do.call(grid.arrange, list(p0a, p0b, p1,p2,p3,p4)) 
dev.off()

#------------------------------------------------------------------------
# plot expression correlation vs. genomic distance
#-----------------------------------------------------------------------
p <- dotplotWithDensityLogXY(revDF(plotDFcis100HiC), "dist", "expCor", "group", COL, fit=TRUE, ALPHA=.2, xlog=TRUE, ylog=FALSE)

pdf(paste0(outPrefix, ".cis100_pairs.expCor_vs_dist.dotplot.pdf"))
    grid::grid.draw(p)
dev.off()

#-----------------------------------------------------------------------
# Expression analysis
#-----------------------------------------------------------------------

paraExpCorSummaryList=list()

# iterate over all expression data sets
for (expName in names(expDFlist)){

    expDF = expDFlist[[expName]]

    # plot distribution of correlation coefficient
    paraExpCor = closePairs[,paste0(expName, "_expCor")]
    sampledExpCor = sampClosePairsCombined[,paste0(expName, "_expCor")]
    
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
    commonEhGroups = getGroupsByCounts(closePairs$commonEnhancer, maxCount)
    
    # make ggplot 
    p = ggplot(closePairs, aes(x=commonEhGroups, y=paraExpCor)) + 
        geom_boxplot(fill=NA, aes(colour = commonEhGroups), lwd=1.5) + 
        geom_jitter(alpha=.3, size=2, aes(colour = commonEhGroups)) + 
        scale_colour_manual(values= colorRampPalette(c("gray", COL_EH))(maxCount+1), guide=FALSE) + theme_bw() + 
        theme(text = element_text(size=20)) + labs(y="Pearson correlation coefficient", x = "Number of shared enhancers")
    ggsave(p, file=paste0(outPrefix,".expCor_by_shared_enhancers.", expName, ".boxplot.pdf"), width=7, height=7)

    #-------------------------------------------------------------------
    # Compare expression against localization in common TAD
    #-------------------------------------------------------------------
    for (tadName in names(allTADs)){
        inTADgroup <- factor(closePairs[, tadName], levels=c(FALSE, TRUE), labels=c("Not in same TAD", "In same TAD"))
        
        # test difference with wilcoxon test
        ws.test = wilcox.test(paraExpCor ~  inTADgroup)
    
        p = ggplot(closePairs, aes(x=inTADgroup, y=paraExpCor)) + 
            geom_boxplot(fill=NA, aes(colour = inTADgroup), lwd=1.5) + 
            geom_jitter(alpha=.3, size=2, aes(colour = inTADgroup)) + 
            scale_colour_manual(values= colorRampPalette(c("gray", COL_TAD[1]))(2), guide=FALSE) + 
            theme_bw() + theme(text = element_text(size=20), axis.text.x = element_text(angle = 45, hjust = 1)) + labs(y="Pearson correlation coefficient", x = "", title=paste0("p = ",signif(ws.test$p.value,2)))
        ggsave(p, file=paste0(outPrefix,".expCor_", expName, ".by_sahred_TAD_", tadName, ".boxplot.pdf"), width=3.5, height=7)
    }

    # plot for all cell types: 
    expPlotDF = data.frame(
        expCor = rep(paraExpCor, length(allTADs)),
        withinTAD = unlist(lapply(names(allTADs), function(tadName) factor(closePairs[, tadName], levels=c(TRUE, FALSE), labels=c("In same TAD", "Not in same TAD")))),
        TadSource = rep(factor(gsub("_", " ", names(allTADs)), gsub("_", " ", names(allTADs))) , each=length(paraExpCor))
    )
    p = ggplot(expPlotDF, aes(x=withinTAD, y=expCor, colour=withinTAD)) + geom_violin(aes(colour = withinTAD), lwd=1, alpha=.25) + geom_boxplot(aes(color=withinTAD), fill=NA, width=.25, lwd=1)  + 
    facet_wrap(~ TadSource) + scale_colour_manual(values= colorRampPalette(c(COL_TAD[1], "gray"))(2), guide=FALSE) + theme_bw() + theme(text = element_text(size=20), axis.text.x = element_text(angle = 45, hjust = 1)) + labs(y="Pearson correlation coefficient", x = "") +  theme(legend.position="bottom")

    # save ggplot as pdf
    ggsave(p, file=paste0(outPrefix,".expCor_", expName, ".by_sahred_TAD.allTADs.boxplot.pdf"), w=7, h=7)
        
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

pdf(paste0(outPrefix,".sampled_genes.Expression_correlation.all_data.violinplot.pdf"))
    ggplot(paraExpCorSummaryDF, na.rm=TRUE, aes(x=genepairs, y = expCor, fill=genepairs)) + geom_violin(adjust = .2) + theme_bw() + labs(y="Pearson correlation coefficient", x = "", title= "Co-expression correlation over tissues") + scale_x_discrete(labels="")  + facet_grid(.~dataset) + theme(legend.position = "bottom") + theme(text = element_text(size=15)) + scale_fill_manual(values=COL) + geom_boxplot(aes(fill=NULL), width=.2, colour="black", outlier.shape = NA)+guides(fill=guide_legend(title=""))
dev.off()

pdf(paste0(outPrefix,".sampled_genes.Expression_correlation.all_data.violinplot.rotated.pdf"))

    ggplot(paraExpCorSummaryDF, na.rm=TRUE, aes(x=genepairs, y = expCor, fill=genepairs)) + geom_violin(adjust = .2) + theme_bw() + labs(y="Pearson correlation coefficient", x = "", title= "Co-expression correlation over tissues") + scale_x_discrete(labels="")  + facet_grid(dataset~.) + coord_flip() + theme(legend.position = "bottom") + scale_fill_manual(values=COL) + geom_boxplot(aes(fill=NULL), width=.2, colour="black", outlier.shape = NA)+guides(fill=guide_legend(title=""))
   
dev.off()

#=======================================================================
# 5.) Compare one2one orthologs of human paralogs in mouse and dog
#=======================================================================

# iterate over all species:
for (orgStr in names(orgStr2Name)){
    
    orgName = orgStr2Name[orgStr]
    
    # check who many pairs have one-to-one orthologs
    pdf(paste0(outPrefix, ".paralogPairs_One2One_orthologs.", orgName, ".barplot.pdf"), width=3.5)
        randPercent = sapply(randPairsInCis, function(gP) percentTrue(gP[,paste0(orgStr, "_one2one")]))
        height = c(percentTrue(allCisPairs[,paste0(orgStr, "_one2one")]), mean(randPercent))
        par(lwd=2, cex=1.5)
        bp = my.barplot(height, 
            names=c("Paralog pairs\n(on same chrom)", "Random pairs\n(on same chrom)"), 
            addValues=TRUE, col=COL_RAND,
            main=paste("One-to-one\n orthologs\n in ", orgName), 
            ylab="Pairs with one-to-one orthologs [%]")
        
        error.bar(bp,height, c(NA, sd(randPercent)))
    dev.off()
    
    # check who many of the ortholog pairs are located on the same chrom.
    pdf(paste0(outPrefix, ".paralogPairs_orthologs_onSameChrom.", orgName, ".barplot.pdf"), width=3.5)

        paraPercent = percentTrue(allCisPairs[allCisPairs[,paste0(orgStr, "_one2one")], paste0(orgStr, "_sameChrom")])
        
        randPercent = sapply(randPairsInCis, function(gP) percentTrue(gP[gP[,paste0(orgStr, "_one2one")], paste0(orgStr, "_sameChrom")]))        
        height = c(paraPercent , mean(randPercent))

        par(lwd=2, cex=1.5)
        bp = my.barplot(height, 
            names=c("Orthologs of\n paralog pairs\n", "Orthologs of\n random pairs\n"), 
            addValues=TRUE, col=c(COL_ORTHO[1], COL_RAND[2]),
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
    
    pdf(paste0(outPrefix, ".paralogPairs_orthologs_dist_", orgName, ".hist.pdf"))
        par(lwd=2, cex=1.5, mfrow=c(2,1))

        hist(orthoDist, 50, col=COL_ORTHO[1],
            main="Distance between orthologs of paralogs", xlab="Distance (kb)")
        hist(randDist, 50, col=COL_RAND[2],
            main="Distance between orthologs of random gene pairs", xlab="Distance (kb)")    
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
        
    pdf(paste0(outPrefix, ".paralogPairs_orthologs_dist_sampled_", orgName, ".hist.pdf"))
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
    pdf(paste0(outPrefix, ".paralogPairs_orthologs_dist_sampled_", orgName, ".boxplot.pdf"))
        par(lwd=2, cex=1.5)
        my.boxplot(list(paraDist, orthoDist, randDist, randOrthoDist), 
        ylab="Distance [kb]", border=c(COL, COL_ORTHO)[c(1,3,2,4)], names=c("paralogs", "orthologs of paralogs", "sampled pairs", "orthologs of\n sampled pairs"))
    dev.off()

    # correlate linear distances of paralogs in human with mouse orthologs
    pdf(paste0(outPrefix, ".paralogPairs_orthologs_dist_", orgName, ".dotplot.pdf"))

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
    pdf(paste0(outPrefix, ".paralogPairs_orthologs_sampled_dist_", orgName, ".dotplot.pdf"))
        
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
        
    pdf(paste0(outPrefix, ".paralogPairs_orthologs_inSameTAD", orgName, ".barplot.pdf"), width=3.5)

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

        pdf(paste0(outPrefix, ".paralogPairs_orthologs_inSameTAD", orgName, ".from_", tadName,".barplot.pdf"), width=3.5)

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
runBasicParalogAnalysis(paste0(outPrefix, ".Mouse"), paralogPairsMouse, "mmusculus_paralog_ds", speciesTssGR[["mmusculus"]], speciesGenesGR[["mmusculus"]], speciesTADs[["mmusculus"]], tissueName="Mouse", HiClist=speciesHiC[["mmusculus"]][[1]], HiClistNorm=speciesHiC[["mmusculus"]][[2]])
message("Finish Mouse paralog analysis!")

runBasicParalogAnalysis(paste0(outPrefix, ".Dog"), paralogPairsDog, "cfamiliaris_paralog_ds", speciesTssGR[["cfamiliaris"]], speciesGenesGR[["cfamiliaris"]], speciesTADs[["cfamiliaris"]], tissueName="Dog", HiClist=speciesHiC[["cfamiliaris"]][[1]], HiClistNorm=speciesHiC[["cfamiliaris"]][[2]])
message("Finish Dog paralog analysis!")

# save temporary image file
save.image(paste0(WORKIMAGE_FILE, ".after_species_analysis.Rdata"))
#~ stop("INFO: Stopped script HERE!")

#=======================================================================
# Compare sahred enhancer positions
#=======================================================================
nEhPos <- 3
nParaPairs <- nrow(closePairs)
nRandPairs <- nrow(sampEhClosePairsCombined)
ehColNames <- c("eh_up", "eh_cent", "eh_down")
ehPosFactor <- factor(c("upstream", "center", "downstream"), levels=c("upstream", "center", "downstream"), labels=c("upstream", "center", "downstream"))

ehPosDF = data.frame(
        "group"=c(rep("paralogs", nParaPairs*nEhPos), rep("sampled", nRandPairs*nEhPos)),
        "enahncer_position"= unlist(list(rep(ehPosFactor, each=nParaPairs), rep(ehPosFactor, each=nRandPairs))),
        "eh_counts"=c(unlist(closePairs[,ehColNames]), unlist(sampEhClosePairsCombined[,ehColNames])),
        "sameIMR90_TAD"=factor(c(rep(closePairs[,"Rao_IMR90"], nEhPos), rep(sampEhClosePairsCombined[,"Rao_IMR90"], nEhPos)), levels=c(TRUE, FALSE), labels=c("same TAD (Rao IMR90)", "not same TAD"))
)

pdf(paste0(outPrefix,"enhancer_position_same_strand.paralogs_sampled.1-5_enhancers.barplot.pdf"))
    ggplot(ehPosDF[!is.na(ehPosDF$eh_count) & ehPosDF$eh_count <= 6 & ehPosDF$eh_count >=1,]) + 
        geom_bar(aes(factor(eh_counts), fill=enahncer_position), width=.6, position="dodge") +
        theme_bw() + facet_grid(group~sameIMR90_TAD, scales="free_y") + 
        labs(y="Counts", x="Number of enhancers", title= "Position of shared enhancer") + scale_fill_manual(values=COL_EH_POS) 
dev.off()

pdf(paste0(outPrefix,"enhancer_position_same_strand.paralogs_sampled.number_of_enhancers.barplot.pdf"))
    ggplot(ehPosDF[!is.na(ehPosDF$eh_count) & ehPosDF$eh_count >=1,]) + 
        geom_bar(aes(enahncer_position, fill=factor(eh_counts), weight=eh_counts)) +
        theme_bw() + facet_grid(group~sameIMR90_TAD, scales="free_y", margins=TRUE) + 
        labs(y="Summed counts", x="Enhancer position", title= "Shared enhancer position relative to genes pairs") + 
        scale_fill_manual(values=colorRampPalette(brewer.pal(9, "Oranges"))(length(unique(ehPosDF$eh_count))), guide = guide_legend("Enhancer\ncounts")) + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf(paste0(outPrefix,"enhancer_position_same_strand.paralogs_sampled.all_summed.barplot.pdf"))
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
nRandPairs <- sum(sapply(sampClosePairsGR, length))


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
        dist = rep(c(closePairs[,"dist"], sampClosePairsCombined[,"dist"]), nExp),
        DupAge = factor(rep(c(closePairs[,AGE_MARK], sampClosePairsCombined[,AGE_MARK]), nExp), c(TRUE, FALSE), c("Young", "Old"))
    )

    # make data.frame with additional columns indicating the groups
    paraEhDF <- data.frame(
        commonEnhancer = c(closePairs[,"commonEnhancer"], sampEhClosePairsCombined[,"commonEnhancer"]),
        genepairs = rep(
            c("Paralog", "Sampled"),
            c(nParaPairs, nRandPairs)
            ),
        DupAge = factor(c(closePairs[,AGE_MARK], sampEhClosePairsCombined[,AGE_MARK]), c(TRUE, FALSE), c("Young", "Old"))
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
    g <- ggplot(paraEhDF, aes(commonEnhancer, group=DupAge, fill=DupAge)) + 
        geom_bar(aes(y = ..density..), width=.5, position="dodge") + xlim(0,10) +
        theme_bw() + 
        labs(y="Density", title= "Number of shared enhancer")  + theme(legend.position = "bottom") + scale_color_manual(values=COL) + scale_fill_manual(values=COL_AGE) + 
        facet_grid(genepairs~.)
    ggsave(paste0(outPrefix,"old_vs_jung.", AGE_MARK, ".sharedEnhancer.pdf"))

    g <- ggplot(paraEhDF, aes(HasCommonEh, group=DupAge, fill=DupAge)) + 
        geom_bar(aes(y = ..density..), width=.5, position="dodge")  +
        theme_bw() + 
        labs(y="Density", title= "Number of shared enhancer")  + theme(legend.position = "bottom") + scale_color_manual(values=COL) + scale_fill_manual(values=COL_AGE) + 
        facet_grid(genepairs~.)
    ggsave(paste0(outPrefix,"old_vs_jung.", AGE_MARK, ".sharedEnhancer.pdf"))

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
            HiCraw = c(distalPairs[,"HiCfreq"], sampDistalPairsCombined[,"HiCfreq"]),
            HiCobsExp = c(distalPairs[,"HiCobsExp"], sampDistalPairsCombined[,"HiCobsExp"]),
            captureC_raw = c(distalPairs[,"captureC_raw"], sampDistalPairsCombined[,"captureC_raw"]),
            captureC_ObsExp = c(distalPairs[,"captureC_ObsExp"], sampDistalPairsCombined[,"captureC_ObsExp"]),
            dist=abs(c(distalPairs[,"dist"], sampDistalPairsCombined[,"dist"]))/10^3,
            distBin=as.factor(breaks[.bincode(abs(c(distalPairs[,"dist"], sampDistalPairsCombined[,"dist"]))/10^3, breaks)])
    )
    hicDF$HiCrawNoZero = hicDF$HiCraw
    hicDF$HiCrawNoZero[hicDF$HiCraw == 0] = NA
    hicDF$HiCobsExpNoZero = hicDF$HiCobsExp
    hicDF$HiCobsExpNoZero[hicDF$HiCobsExp == 0] = NA

    g <- ggplot(hicDF, aes(x=group, y=HiCraw, colour = DupAge))  +
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
# LRRC8 Gene family
#=======================================================================

# get genomic interval arround LRRC8 geens:
lrrc8GR = tssGR[tssGR$hgnc_symbol=="LRRC8C",]
lrrc8GR = resize(lrrc8GR, 2*10^6, fix="center")

# get gene names with TSS inside this interval
#geneSymbols = c("LRRC8B", "LRRC8C", "LRRC8D")
#~ geneSymbols = subsetByOverlaps(tssGR,lrrc8GR, ignore.strand=TRUE)$hgnc_symbol
exampleGenes = subsetByOverlaps(tssGR,lrrc8GR, ignore.strand=TRUE)
geneSymbols = exampleGenes$hgnc_symbol

# define colors for paralog groups
geneCols = rep("black", length(geneSymbols))
geneCols[grep("^LRRC8.*", geneSymbols)] = COL_FAMILY[3]
geneCols[grep("^GBP.*", geneSymbols)] = COL_FAMILY[5]

matEnhancer = sharedEnhancerMatrix(names(exampleGenes), gene2ehID, geneSymbols=geneSymbols)
matDist = pairwiseDistMatrix(exampleGenes, geneSymbols=geneSymbols)/10^3


matContacts = pairwiseContacstMatrix(exampleGenes, HiClist)
diag(matContacts) = NA
matContactsNorm = pairwiseContacstMatrix(exampleGenes, HiClistNorm)
diag(matContactsNorm) = NA

# query pairwise promoter-promoter contacts
matCaptureCContacts = as.matrix(captureHiC[["raw"]][names(exampleGenes),names(exampleGenes)])
diag(matCaptureCContacts) = NA
matCaptureCContactsObsExp = as.matrix(captureHiC[["obsExp"]][names(exampleGenes),names(exampleGenes)])
diag(matCaptureCContactsObsExp) = NA

pdf(paste0(outPrefix, ".paralogPairs_LRRC8_common_enhancer.heatmap.pdf"))
    my.heatmap.2(matEnhancer, cellnote=matEnhancer, Rowv=FALSE, Colv=FALSE, 
        revC=FALSE, dendrogram="none", 
        notecex=1.5, notecol="black", trace="none",
        ClabSide=3, RlabSide=2, margins=c(2, 2), key.xlab="Shared enhancers", main="Shared enhancers", col=brewer.pal(9,"Blues"), ClabColor=geneCols, RlabColor=geneCols)
dev.off()

pdf(paste0(outPrefix, ".paralogPairs_LRRC8_dist.heatmap.pdf"))
    my.heatmap.2(matDist, Rowv=FALSE, Colv=FALSE, 
        revC=FALSE, dendrogram="none", 
        notecex=1, notecol="black", trace="none",
        ClabSide=3, RlabSide=2, margins=c(2, 2), key.xlab="Distance (kb)", main="Linear distance", col=colorRampPalette(rev(brewer.pal(9,"Blues"))), ClabColor=geneCols, RlabColor=geneCols 
    )
dev.off()

# Hi-C contacts has pairwise heatmaps
pdf(paste0(outPrefix, ".paralogPairs_LRRC8_Hi-C_raw_contact.heatmap.pdf"))
    my.heatmap.2(matContacts, Rowv=FALSE, Colv=FALSE, 
        revC=FALSE, dendrogram="none", 
        notecex=1, notecol="black", trace="none",
        ClabSide=3, RlabSide=2, margins=c(2, 2), key.xlab="Hi-C contacs", main="Hi-C contacts of TSS pairs", col=colorRampPalette(brewer.pal(9,"Reds")), na.color="darkgray", ClabColor=geneCols, RlabColor=geneCols
    )
dev.off()

pdf(paste0(outPrefix, ".paralogPairs_LRRC8_Hi-C_normalized_contact.heatmap.pdf"))
    my.heatmap.2(log2(matContactsNorm), Rowv=FALSE, Colv=FALSE, 
        revC=FALSE, dendrogram="none", 
        notecex=1, notecol="black", trace="none",
        ClabSide=3, RlabSide=2, margins=c(2, 2), key.xlab="log2(O/E) Hi-C", main="Normalized Hi-C contacts\n of TSS pairs", col=colorRampPalette(rev(brewer.pal(11,"RdBu"))), na.color="darkgray", ClabColor=geneCols, RlabColor=geneCols
    )
dev.off()

# plot Capture Hi-C contacts
pdf(paste0(outPrefix, ".paralogPairs_LRRC8_CaptureHiC_raw_contact.heatmap.pdf"))
    my.heatmap.2(matCaptureCContacts, Rowv=FALSE, Colv=FALSE, 
        revC=FALSE, dendrogram="none", 
        notecex=1, notecol="black", trace="none", labRow=geneSymbols, labCol=geneSymbols,
        ClabSide=3, RlabSide=2, margins=c(2, 2), key.xlab="Capture Hi-C contacs", main="Capture Hi-C contacts of TSS pairs", col=colorRampPalette(brewer.pal(9,"Reds")), na.color="darkgray", ClabColor=geneCols, RlabColor=rev(geneCols)
    )
dev.off()

pdf(paste0(outPrefix, ".paralogPairs_LRRC8_CaptureHiC_ObsExp_contact.heatmap.pdf"))
    my.heatmap.2(matCaptureCContactsObsExp, Rowv=FALSE, Colv=FALSE, 
        revC=FALSE, dendrogram="none", labRow=geneSymbols, labCol=geneSymbols,
        notecex=1, notecol="black", trace="none",
        ClabSide=3, RlabSide=2, margins=c(2, 2), key.xlab="log2(O/E) Capture Hi-C", main="Normalized Capture Hi-C contacts\n of TSS pairs", col=colorRampPalette(rev(brewer.pal(11,"RdBu"))), na.color="darkgray", ClabColor=geneCols, RlabColor=geneCols
    )
dev.off()

# iterate over all expression data sets
for (i in seq(nExp)){

    expName = names(expDFlist)[i]
    expDF = expDFlist[[i]]

    matExpression = getCor(cbind(names(exampleGenes), names(exampleGenes)), expDF)
    dimnames(matExpression)<- list(geneSymbols, geneSymbols)
    diag(matExpression) = NA
    
    pdf(paste0(outPrefix, ".paralogPairs_LRRC8.", expName, ".expression_correlation.heatmap.pdf"))
        my.heatmap.2(matExpression, Rowv=FALSE, Colv=FALSE, 
            revC=FALSE, dendrogram="none", breaks = seq(-1,1,len=21),
            notecex=1, notecol="black", trace="none",
            ClabSide=3, RlabSide=2, margins=c(2, 2), key.xlab="Pearson R", main=paste0("Gene expression correlation\n over n=", ncol(expDF), " tissues/cells in \n", expName), col=colorRampPalette(rev(brewer.pal(11,"RdBu"))),tracecol="darkgreen", na.color="darkgray", ClabColor=geneCols, RlabColor=geneCols)
    dev.off()
    
    pdf(paste0(outPrefix, ".paralogPairs_LRRC8.", expName, ".expression_correlation_dendrogram.heatmap.pdf"))
        my.heatmap.2(matExpression, revC=TRUE,
            notecex=1, notecol="black", trace="none", breaks = seq(-1,1,len=21),
            margins=c(8, 8), key.xlab="Pearson R", main=paste0("Gene expression correlation\n over n=", ncol(expDF), " tissues/cells in \n", expName), col=colorRampPalette(rev(brewer.pal(11,"RdBu"))),tracecol="darkgreen", na.color="darkgray", ClabColor=geneCols, RlabColor=geneCols)
    dev.off()
    
    # write subset of expression to output table
    expSubset = expDF[names(exampleGenes),]
    row.names(expSubset) = tssGR[names(exampleGenes)]$hgnc_symbol
    write.table(expSubset, paste0(outPrefix, ".paralogPairs_LRRC8.", expName, ".expression_values.csv"), sep="\t", quote=FALSE, col.names=NA)
    
    # heatmap of expression over tissues
    pdf(paste0(outPrefix, ".paralogPairs_LRRC8.", expName, ".expression_tissues.heatmap.pdf"))
        
        # get breaks:
        if (max(expSubset)>50) {
            breaks = seq(0,50,len=21)
        }else{
            breaks = 21
        }
        my.heatmap.2(as.matrix(expSubset), 
            Rowv=FALSE, Colv=FALSE, revC=FALSE, trace="none", dendrogram="none", 
            col=colorRampPalette(brewer.pal(9,"Blues"))(20), breaks=breaks,
            notecex=1, notecol="black",
            RlabColor=geneCols, ClabSide=3, RlabSide=2, margins=c(2, 2),
            key.xlab="Expression level\n[FPKM]", main=paste("Expression level in\n",gsub('_', ' ', expName))
        )
    dev.off()
}

#
# Hi-C plot of the LRRC8C region

getHiCregion <- function(HiClist, gr, MARGIN=c(1,2), ...){
    
    chr = as.character(seqnames(gr))
    mapName = paste0(chr, chr)
    
    extractRegion(HiClist[[mapName]], MARGIN, chr, start(gr), end(gr), ...)
    
}

# get subset of the interaction maps for the target region
subMap = getHiCregion(HiClist, lrrc8GR)
subMapNorm = getHiCregion(HiClistNorm, lrrc8GR)


# write browser track in long-range interaction format for the subMap
subMapLarge =  getHiCregion(HiClist, resize(lrrc8GR, 6*10^6, fix="center"))
writeHiCinteractions(list(subMapLarge), paste0(outPrefix, ".paralogPairs_LRRC8_Hi-C_IMR90_50kb.map.track.txt"))

# binn the data for plotting
subMap.binned = binningC(subMap, binsize=HIC_RESOLUTION, bin.adjust=FALSE, optimize.by="speed")
subMapNorm.binned = binningC(subMapNorm, binsize=HIC_RESOLUTION, bin.adjust=FALSE, optimize.by="speed")
intdata(subMapNorm.binned) = log2(intdata(subMapNorm.binned))

# plot the interaction map with annotation tracks


alteringStrand <- function(gr){
    gr = sort(gr)
    strand(gr) = rep(c('+', '-'), times=length(gr))[1:length(gr)]    
    return(gr)
}

MAXCOUNT = 1000
plotDomains = alteringStrand(allTADs[["Rao_IMR90"]])
pdf(paste0(outPrefix, ".paralogPairs_LRRC8_Hi-C_IMR90_50kb.map.pdf"))
    mapC(subMap.binned, subMapNorm.binned, 
        minrange=0, maxrange=MAXCOUNT, grid=TRUE,
        tracks=list("TSS"=tssGR, "Enhancers"=ehGR, "TADs"=plotDomains),
        title="Chromatin contacs at the LRRC8 locus \nIMR90 in situ Hi-C at 50kb resolution from Rao et al. 2014"
        )
    legend("bottomright", legend=MAXCOUNT, fill="red", border="red", cex=1.5)
dev.off()

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
# Olfactory receptors:
#=======================================================================
orPairIDX <- which(closePairs[,1] %in% ORids & closePairs[,2] %in% ORids)
orDistalIDX <- which(distalPairs[,1] %in% ORids & distalPairs[,2] %in% ORids)


# get subset of pairs taht are OR
orPairs <- closePairs[orPairIDX,]
orPairGR <- closePairsGR[orPairIDX]

orDistalPairs <- distalPairs[orDistalIDX,]

# Make data.frame for plotting with ggplot2
breaks = seq(DISTAL_MIN_DIST, DISTAL_MAX_DIST, length.out=10) / 10^3
plotDF = data.frame(
        group = c(rep("OR", nrow(orDistalPairs)), rep("paralogs", nrow(distalPairs)), rep("sampled", nrow(sampDistalPairsCombined))),
        HiCraw = c(orDistalPairs[,"HiCfreq"], distalPairs[,"HiCfreq"], sampDistalPairsCombined[,"HiCfreq"]),
        HiCobsExp = c(orDistalPairs[,"HiCobsExp"], distalPairs[,"HiCobsExp"], sampDistalPairsCombined[,"HiCobsExp"]),
        captureC_raw = c(orDistalPairs[,"captureC_raw"], distalPairs[,"captureC_raw"], sampDistalPairsCombined[,"captureC_raw"]),
        captureC_ObsExp = c(orDistalPairs[,"captureC_ObsExp"],  distalPairs[,"captureC_ObsExp"], sampDistalPairsCombined[,"captureC_ObsExp"]),
        dist=abs(c(orDistalPairs[,"dist"], distalPairs[,"dist"], sampDistalPairsCombined[,"dist"]))/10^3,
        distBin=as.factor(breaks[.bincode(abs(c(orDistalPairs[,"dist"], distalPairs[,"dist"], sampDistalPairsCombined[,"dist"]))/10^3, breaks)])
    )


plotDF$HiCrawNoZero = plotDF$HiCraw
plotDF$HiCrawNoZero[plotDF$HiCraw == 0] = NA
plotDF$HiCobsExpNoZero = plotDF$HiCobsExp
plotDF$HiCobsExpNoZero[plotDF$HiCobsExp == 0] = NA

groupCol = c(COL_FAMILY[6], COL)
p =  ggplot(plotDF, aes(x=group, y=HiCobsExp, colour = group)) + scale_y_log10() +
    geom_boxplot(lwd=1.5) + 
    scale_color_manual(values=groupCol, guide=FALSE) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + 
    guides(fill=guide_legend(title="")) +
    labs(y="Normalized Hi-C [obs/exp]", x="")

# save ggplot as pdf
ggsave(p, file=paste0(outPrefix, "OR_genes.distal_pairs.Hi-C_normalized_contacts.boxplot.pdf"), w=3.5, h=7)

#-----------------------------------------------------------------------
# analysie close OR pairs
#-----------------------------------------------------------------------

write.table(orPairs, col.names=TRUE, row.names=FALSE, file=paste0(outPrefix, ".close_pairs_OR.annotation.txt"), sep="\t", quote=FALSE)


#~ plot(start(orPairGR[seqnames(orPairGR) == "chr1"]), ylim=c(1.58*10^8, 1.6*10^8))
#~ plot(start(orPairGR[seqnames(orPairGR) == "chr11"]), ylim=c(0, 1*10^7))

exampleLoci <- list(
    "OR_cluster_1.4"=GRanges("chr1", IRanges(1.58*10^8, 1.6*10^8), seqinfo=seqInfo),
    "OR_cluster_1.5"=GRanges("chr1", IRanges(2.47*10^8, 2.49*10^8), seqinfo=seqInfo),
    "OR_cluster_11.3"=GRanges("chr11", IRanges(4*10^6, 1*10^7), seqinfo=seqInfo),
    "OR_cluster_11.11"=GRanges("chr11", IRanges(5.4*10^7, 6*10^7), seqinfo=seqInfo)
)

for (exampleName in names(exampleLoci)){

    exampleGR <- exampleLoci[[exampleName]]

    # get gene names with TSS inside this interval
    exampleTSS = subsetByOverlaps(tssGR, exampleGR, ignore.strand=TRUE)
    exampleGenes = subsetByOverlaps(genesGR, exampleGR, ignore.strand=TRUE)
    geneSymbols = exampleTSS$hgnc_symbol
    
    # define colors for paralog groups
    geneCols = rep("black", length(geneSymbols))
    #~ geneCols[grep("^OR1.*", geneSymbols)] = COL_FAMILY[6]
    #~ geneCols[grep("^OR2.*", geneSymbols)] = COL_FAMILY[7]
    candidateIDs <- which(names(exampleTSS) %in% ORids)
    geneCols[candidateIDs] = COL_FAMILY[6]
    candidateGenes <- exampleGenes[candidateIDs]
    
    matEnhancer = sharedEnhancerMatrix(names(exampleTSS), gene2ehID, geneSymbols=geneSymbols)
    matDist = pairwiseDistMatrix(exampleTSS, geneSymbols=geneSymbols)/10^3
    
    # get Hi-C data
    matContacts = pairwiseContacstMatrixSameChrom(exampleTSS, HiClist)
    matContactsNorm = pairwiseContacstMatrixSameChrom(exampleTSS, HiClistNorm)
    
    dimNames <- exampleTSS$hgnc_symbol
    dimNames[dimNames == ""] <- names(exampleTSS)[dimNames == ""]
    dimnames(matContacts) <- list(dimNames, dimNames)
    diag(matContacts) <- NA
    dimnames(matContactsNorm) <- list(dimNames, dimNames)
    diag(matContactsNorm) = NA
    
    # query pairwise promoter-promoter contacts
    matCaptureCContacts = as.matrix(captureHiC[["raw"]][names(exampleTSS),names(exampleTSS)])
    diag(matCaptureCContacts) = NA
    matCaptureCContactsObsExp = as.matrix(captureHiC[["obsExp"]][names(exampleTSS),names(exampleTSS)])
    diag(matCaptureCContactsObsExp) = NA
    
    pdf(paste0(outPrefix, ".paralogPairs.", exampleName, ".common_enhancer.heatmap.pdf"))
        if(nrow(matEnhancer)<20){
            labelValues <- matEnhancer
        }else{
            labelValues <- matrix(rep(NA, nrow(matEnhancer)^2), nrow(matEnhancer))
        }
        my.heatmap.2(matEnhancer, cellnote=labelValues, Rowv=FALSE, Colv=FALSE, 
            revC=FALSE, dendrogram="none", 
            notecex=1.5, notecol="black", trace="none",
            ClabSide=3, RlabSide=2, margins=c(2, 2), key.xlab="Shared enhancers", main="Shared enhancers", col=brewer.pal(9,"Blues"), ClabColor=geneCols, RlabColor=geneCols)
    dev.off()
    
    pdf(paste0(outPrefix, ".paralogPairs.", exampleName, ".dist.heatmap.pdf"))
        my.heatmap.2(matDist, Rowv=FALSE, Colv=FALSE, 
            revC=FALSE, dendrogram="none", 
            notecex=1, notecol="black", trace="none",
            ClabSide=3, RlabSide=2, margins=c(2, 2), key.xlab="Distance (kb)", main="Linear distance", col=colorRampPalette(rev(brewer.pal(9,"Blues"))), ClabColor=geneCols, RlabColor=geneCols 
        )
    dev.off()
    
    # Hi-C contacts has pairwise heatmaps
    pdf(paste0(outPrefix, ".paralogPairs.", exampleName, ".Hi-C_raw_contact.heatmap.pdf"))
        my.heatmap.2(matContacts, Rowv=FALSE, Colv=FALSE, 
            revC=FALSE, dendrogram="none", 
            notecex=1, notecol="black", trace="none",
            ClabSide=3, RlabSide=2, margins=c(2, 2), key.xlab="Hi-C contacs", main="Hi-C contacts of TSS pairs", col=colorRampPalette(brewer.pal(9,"Reds")), na.color="darkgray", ClabColor=geneCols, RlabColor=geneCols
        )
    dev.off()
    
    pdf(paste0(outPrefix, ".paralogPairs.", exampleName, ".Hi-C_normalized_contact.heatmap.pdf"))
        my.heatmap.2(log2(matContactsNorm+.01), Rowv=FALSE, Colv=FALSE, 
            revC=FALSE, dendrogram="none", 
            notecex=1, notecol="black", trace="none",
            ClabSide=3, RlabSide=2, margins=c(2, 2), key.xlab="log2(O/E) Hi-C", main="Normalized Hi-C contacts\n of TSS pairs", col=colorRampPalette(rev(brewer.pal(11,"RdBu"))), na.color="darkgray", ClabColor=geneCols, RlabColor=geneCols, cexCol=30/nrow(matContactsNorm)
        )
    dev.off()
    
    # plot Capture Hi-C contacts
    pdf(paste0(outPrefix, ".paralogPairs.", exampleName, ".CaptureHiC_raw_contact.heatmap.pdf"))
        my.heatmap.2(matCaptureCContacts, Rowv=FALSE, Colv=FALSE, 
            revC=FALSE, dendrogram="none", 
            notecex=1, notecol="black", trace="none", labRow=geneSymbols, labCol=geneSymbols,
            ClabSide=3, RlabSide=2, margins=c(2, 2), key.xlab="Capture Hi-C contacs", main="Capture Hi-C contacts of TSS pairs", col=colorRampPalette(brewer.pal(9,"Reds")), na.color="darkgray", ClabColor=geneCols, RlabColor=rev(geneCols)
        )
    dev.off()
    
    pdf(paste0(outPrefix, ".paralogPairs.", exampleName, ".CaptureHiC_ObsExp_contact.heatmap.pdf"))
        my.heatmap.2(matCaptureCContactsObsExp, Rowv=FALSE, Colv=FALSE, 
            revC=FALSE, dendrogram="none", 
            notecex=1, notecol="black", trace="none", labRow=geneSymbols, labCol=geneSymbols,
            ClabSide=3, RlabSide=2, margins=c(2, 2), key.xlab="log2(O/E) Capture Hi-C", main="Normalized Capture Hi-C contacts\n of TSS pairs", col=colorRampPalette(rev(brewer.pal(11,"RdBu"))), na.color="darkgray", ClabColor=geneCols, RlabColor=rev(geneCols)
        )
    dev.off()
    
    # iterate over all expression data sets
    for (i in seq(nExp)){
    
        expName = names(expDFlist)[i]
        expDF = expDFlist[[i]]
    
        matExpression = getCor(cbind(names(exampleTSS), names(exampleTSS)), expDF)
        dimnames(matExpression)<- list(geneSymbols, geneSymbols)
        diag(matExpression) = NA
        
        pdf(paste0(outPrefix, ".paralogPairs.", exampleName, ".", expName, ".expression_correlation.heatmap.pdf"))
            my.heatmap.2(matExpression, Rowv=FALSE, Colv=FALSE, 
                revC=FALSE, dendrogram="none", breaks = seq(-1,1,len=21),
                notecex=1, notecol="black", trace="none",
                ClabSide=3, RlabSide=2, margins=c(2, 2), key.xlab="Pearson R", main=paste0("Gene expression correlation\n over n=", ncol(expDF), " tissues/cells in \n", expName), col=colorRampPalette(rev(brewer.pal(11,"RdBu"))),tracecol="darkgreen", na.color="darkgray", ClabColor=geneCols, RlabColor=geneCols)
        dev.off()
        
        if(!any(is.na(matExpression))){
            pdf(paste0(outPrefix, ".paralogPairs.", exampleName, ".", expName, ".expression_correlation_dendrogram.heatmap.pdf"))
                my.heatmap.2(matExpression, revC=TRUE,
                    notecex=1, notecol="black", trace="none", breaks = seq(-1,1,len=21),
                    margins=c(8, 8), key.xlab="Pearson R", main=paste0("Gene expression correlation\n over n=", ncol(expDF), " tissues/cells in \n", expName), col=colorRampPalette(rev(brewer.pal(11,"RdBu"))),tracecol="darkgreen", na.color="darkgray", ClabColor=geneCols, RlabColor=geneCols)
            dev.off()
        }
        
        # write subset of expression to output table
        expSubset = expDF[names(exampleTSS),]
        rowNames <- tssGR[names(exampleTSS)]$hgnc_symbol
        rowNames[rowNames == ""] <- names(exampleTSS)[rowNames == ""]
        row.names(expSubset) = rowNames
        write.table(expSubset, paste0(outPrefix, ".paralogPairs.", exampleName, ".", expName, ".expression_values.csv"), sep="\t", quote=FALSE, col.names=NA)
        
        # heatmap of expression over tissues
        pdf(paste0(outPrefix, ".paralogPairs.", exampleName, ".", expName, ".expression_tissues.heatmap.pdf"))
            
            # get breaks:
            if (max(expSubset, na.rm=TRUE)>50) {
                breaks = seq(0,50,len=21)
            }else{
                breaks = 21
            }
            my.heatmap.2(as.matrix(expSubset), 
                Rowv=FALSE, Colv=FALSE, revC=FALSE, trace="none", dendrogram="none", 
                col=colorRampPalette(brewer.pal(9,"Blues"))(20), breaks=breaks,
                notecex=1, notecol="black",
                RlabColor=geneCols, ClabSide=3, RlabSide=2, margins=c(2, 2),
                key.xlab="Expression level\n[FPKM]", main=paste("Expression level in\n",gsub('_', ' ', expName))
            )
        dev.off()
    }

    #
    # Hi-C plot of the LRRC8C region

    # get subset of the interaction maps for the target region
    subMap = getHiCregion(HiClist, exampleGR)
    subMapNorm = getHiCregion(HiClistNorm, exampleGR)


    # write browser track in long-range interaction format for the subMap
#~     subMapLarge =  getHiCregion(HiClist, exampleGRresize(exampleGR, 6*10^6, fix="center"))
    writeHiCinteractions(list(subMap), paste0(outPrefix, ".paralogPairs.", exampleName, ".Hi-C_IMR90_50kb.map.track.txt"))
    
    # binn the data for plotting
    subMap.binned = binningC(subMap, binsize=HIC_RESOLUTION, bin.adjust=FALSE, optimize.by="speed")
    subMapNorm.binned = binningC(subMapNorm, binsize=HIC_RESOLUTION, bin.adjust=FALSE, optimize.by="speed")
    intdata(subMapNorm.binned) = log2(intdata(subMapNorm.binned))

    # plot the interaction map with annotation tracks
    MAXCOUNT = 1000
    plotDomains = alteringStrand(allTADs[["Rao_IMR90"]])
    pdf(paste0(outPrefix, ".paralogPairs.", exampleName, ".Hi-C_IMR90_50kb.map.pdf"), w=14, h=7)
        mapC(subMap.binned, 
            minrange=0, maxrange=MAXCOUNT, grid=TRUE,
            tracks=list("candidates"=candidateGenes, "genes"=exampleGenes, "Enhancers"=ehGR, "TADs"=plotDomains),
            title=paste("Chromatin contacs at", exampleName, "locus:", getLocStr(exampleGR), "\nIMR90 in situ Hi-C at 50kb resolution from Rao et al. 2014")
            )
#~         legend("topright", legend=MAXCOUNT, fill="red", border="red", cex=1.5)
    dev.off()

    pdf(paste0(outPrefix, ".paralogPairs.", exampleName, ".Hi-C_IMR90_50kb.map_norm.pdf"))
        mapC(subMap.binned, subMapNorm.binned, 
            minrange=0, maxrange=MAXCOUNT, grid=TRUE,
            tracks=list("candidates"=candidateGenes, "genes"=exampleGenes, "Enhancers"=ehGR, "TADs"=plotDomains),
            title=paste("Chromatin contacs at", exampleName, "locus \nIMR90 in situ Hi-C at 50kb resolution from Rao et al. 2014")
            )
        legend("bottomright", legend=MAXCOUNT, fill="red", border="red", cex=1.5)
    dev.off()

}

#=======================================================================
# make data.frame with all annotations
#=======================================================================

# combine all parlogs and all randomly sampled pairs
paralogPairsUniq$g1 <- paralogPairsUniq[,1]
paralogPairsUniq$g2 <- paralogPairsUniq[,2]

commonCols <- intersect(names(paralogPairsUniq), names(randPairsCombined))


groupFactor <- factor(rep(c("paralog", "random"), c(nrow(paralogPairsUniq), nrow(randPairsCombined))), levels=c("paralog", "random"))

allPairs <- rbind.fill(paralogPairsUniq[,commonCols], randPairsCombined[,commonCols])
allPairs$group <- groupFactor
# add flag for cis pairs
allPairs$sameChrom <- !is.na(allPairs[,"dist"])

#=======================================================================
# combine all parlogs and all randomly sampled pairs
#=======================================================================
# DIMENSION:
# - paralog vs sampled      (2)
# - close vs distal         (2)
# - TAD source              (10)
# - Expression              (4)
# - other species           (2)


closePairs$g1 <- closePairs[,1]
closePairs$g2 <- closePairs[,2]
distalPairs$g1 <- distalPairs[,1]
distalPairs$g2 <- distalPairs[,2]


allCombined <- rbind.fill(closePairs, distalPairs, sampClosePairsCombined, sampDistalPairsCombined)

# add noZero HI-C 
allCombined <- addNoZero(allCombined)


groupFactor <- factor(rep(c("paralog", "sampled"), c(nrow(closePairs)+nrow(distalPairs), nrow(sampClosePairsCombined)+nrow(sampDistalPairsCombined))), levels=c("paralog", "sampled"))
distGroupFactor <- factor(c(  
    rep(c("close", "distal"), c(nrow(closePairs), nrow(distalPairs))),
    rep(c("close", "distal"), c(nrow(sampClosePairsCombined), nrow(sampDistalPairsCombined)))
))
tadFactor <- factor(names(allTADs), levels=names(allTADs))
expFactor <- factor(names(expDFlist), levels=names(expDFlist))
speciesFactor <- factor(orgStr2Name, orgStr2Name)

nSpecies  <- length(orgStr2Name) 
nTAD <- length(allTADs)

allDF <- data.frame(
    group = rep(groupFactor, nTAD*nExp*nSpecies),
    distGroup = rep(distGroupFactor, nTAD*nExp*nSpecies),
    dist = rep(abs(allCombined[,"dist"])/10^3, nTAD*nExp*nSpecies),
    strand = factor(rep(allCombined[,"sameStrand"], nTAD*nExp*nSpecies), c(TRUE, FALSE), c("same", "opposite")),
    tadSource = rep(rep(tadFactor, each=nrow(allCombined)), nExp*nSpecies),
    inTAD = rep(unlist(lapply(names(allTADs), function(t) allCombined[,t])), nExp*nSpecies),
    expCor = rep(unlist(lapply(names(expDFlist), function(e) rep(allCombined[,paste0(e,"_expCor")], nTAD))), nSpecies),
    expSource = rep(rep(expFactor, each=nrow(allCombined)*nTAD), nSpecies),
    species = rep(speciesFactor, each=nrow(allCombined)*nTAD*nExp)
)

speciesCols <- c("one2one", "sameChrom", "dist", "TAD", "HiC", "HiCnorm", "commonOrtholg", "NotOne2one")

speciesDF <- do.call("rbind", lapply(names(orgStr2Name), function(spec){
 spDF <- allCombined[,paste0(spec, "_", speciesCols)]
 names(spDF) <- paste0("ortholog_", speciesCols)
 return(spDF[rep(1:nrow(allCombined), nTAD*nExp),])
 }))

commonCols <- c("g1", "g2", "HGNC_g1", "HGNC_g2", "hsapiens_paralog_perc_id", "hsapiens_paralog_perc_id_r1", "hsapiens_paralog_dn", "hsapiens_paralog_ds", 
"sameStrand",
HiCcolumns,
"comp_combination", "same_comp_region", "common_comp",
"subcomp_combination", "same_subcomp_region", "common_subcomp", "age", "g1_exp_IMR90", "g2_exp_IMR90")

commonDF <- allCombined[rep(1:nrow(allCombined), nTAD*nExp*nSpecies), commonCols]

allDF <- cbind(allDF, commonDF, speciesDF)

write.table(allDF, file=paste0(outPrefix, ".allDF.csv"),
    sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#~ allDF <- read.table(paste0(outPrefix, ".allDF.csv"), header=TRUE, sep="\t")

#=======================================================================
# combine all parlogs and all randomly sampled pairs
#=======================================================================

ehCombined <- rbind.fill(closePairs, sampEhClosePairsCombined)

groupFactor <- factor(rep(c("paralog", "sampled"), c(nrow(closePairs), nrow(sampEhClosePairsCombined))), levels=c("paralog", "sampled"))
tadFactor <- factor(names(allTADs), levels=names(allTADs))
expFactor <- factor(names(expDFlist), levels=names(expDFlist))

nTAD <- length(allTADs)

ehDF <- data.frame(
    group = rep(groupFactor, nTAD*nExp),
    dist = rep(abs(ehCombined[,"dist"])/10^3, nTAD*nExp),
    strand = factor(rep(ehCombined[,"sameStrand"], nTAD*nExp), c(TRUE, FALSE), c("same", "opposite")),
    tadSource = rep(rep(tadFactor, each=nrow(ehCombined)), nExp),
    inTAD = rep(unlist(lapply(names(allTADs), function(t) ehCombined[,t])), nExp),
    expCor = unlist(lapply(names(expDFlist), function(e) rep(ehCombined[,paste0(e,"_expCor")], nTAD))),
    expSource = rep(expFactor, each=nrow(ehCombined)*nTAD)
)

commonEhCols <- intersect(names(closePairs), names(sampEhClosePairsCombined))
commonEhDF <- ehCombined[rep(1:nrow(ehCombined), nTAD*nExp), commonEhCols]

ehDF <- cbind(ehDF, commonEhDF)

write.table(ehDF, file=paste0(outPrefix, ".ehDF.csv"),
    sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#=======================================================================
#=======================================================================

# debug pot
#ggplot(allDF, aes(x=abs(dist))) + geom_histogram() + facet_grid(group ~ distGroup, scales="free_y") + scale_x_log10()
#~ subDF <- subset(allDF, distGroup=="distal" & tadSource=="stable_TADs" & species=="mouse")
subDF <- subset(allDF, tadSource=="stable_TADs" & species=="mouse")

g <- ggplot(subDF, aes(x=group, y = expCor, fill=group)) + 
    geom_violin(aes(color = group), adjust = .4)  + geom_boxplot(fill=NA, width=.25, lwd=1, outlier.shape = NA) + 
    facet_grid(distGroup~expSource) + scale_fill_manual(values=COL) + scale_color_manual(values=COL) + 
    theme_bw()+ theme(legend.position = "none", axis.text.x=element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=14)) + 
    labs(y="Pearson correlation coefficient", x = "", title= "Co-expression correlation over tissues")  #+ theme(legend.position = "bottom") + guides(fill=guide_legend(title=""))

#~ + geom_boxplot(aes(fill=NULL), width=.2, colour="black", outlier.shape = NA) 

ggsave(paste0(outPrefix,".allDF.expCor.by_close_and_distal.violinplot.pdf"), g)

# close paralogs Hi-C by inTAD
subDF <- subset(allDF, distGroup=="close" & expSource=="GTEx" & species=="mouse")
pdf(paste0(outPrefix,".allDF.close.Hi-C_norm.byTAD.boxplot.pdf"), w=14,h=7)
    ggplot(subDF, aes(x=inTAD, fill=group, y=HiCobsExp)) + 
        geom_boxplot() + labs(y="Normalized Hi-C") + scale_y_log10()  + 
        facet_grid(.~tadSource) + theme_bw() + theme(legend.position = "bottom") + theme(text = element_text(size=15)) + scale_fill_manual(values=COL)
dev.off() 
       
pdf(paste0(outPrefix,".allDF.close.Captrue_HiC_norm.byTAD.boxplot.pdf"), w=14,h=7)
    ggplot(subDF, aes(x=inTAD, fill=group, y=captureC_ObsExp)) + 
        geom_boxplot() + theme_bw() + labs(y="Captrue HiC (obs/exp)") + scale_y_log10()  + 
        facet_grid(.~tadSource) + theme(legend.position = "bottom") + theme(text = element_text(size=15)) + scale_fill_manual(values=COL)
dev.off() 

# in age by TAD

subDF <- subset(allDF, group=="paralog" & distGroup=="close" & expSource=="GTEx" & species=="mouse")
subDF <- subDF[subDF$age %in% 0:10,]
pdf(paste0(outPrefix,".allDF.close.age_by_TAD.boxplot.pdf"), w=14,h=7)

    p.vals <- sapply(names(allTADs), function(s) wilcox.test(age ~ inTAD, data=subDF[subDF$tadSource==s,])$p.value)
    p.val.df <- data.frame(tadSource=names(allTADs), pval=p.vals, age=.9*max(subDF$age, na.rm=TRUE), inTAD=TRUE)

    ggplot(subDF, aes(x=tadSource, y=age, fill=inTAD)) +
        geom_boxplot() + theme_bw() + theme(legend.position = "bottom") + theme(text = element_text(size=15)) + scale_fill_manual(values=rev(COL_TAD)) + geom_text(aes(label=paste0("p=", signif(pval, 3))), data=p.val.df) + ylab("Duplication age")
dev.off()

require(colorspace)
# as barplot
d <- ddply(subDF, .(tadSource, age), summarize, percentInTAD=percentTrue(inTAD), n=length(inTAD))

pdf(paste0(outPrefix,".allDF.close.inTAD_by_age.stable_TADs.barplot.pdf"), w=7,h=3.5)
    ggplot(d[d$tadSource=="stable_TADs",], aes(x=factor(age), y=percentInTAD, fill=factor(age))) +
        geom_bar(stat="identity") + theme_bw() + theme(text = element_text(size=15)) + scale_fill_manual(values=terrain_hcl(12)) + 
        geom_text(aes(label=paste0(c("n=", rep("", 10)), n), y=-2)) +
        geom_text(aes(label=paste0(signif(percentInTAD, 2), "%")), , vjust=-0.25) + ylab("% in same TAD")
        
dev.off()

pdf(paste0(outPrefix,".allDF.close.inTAD_by_age.Dixon_hESC.barplot.pdf"), w=7,h=3.5)
    ggplot(d[d$tadSource=="Dixon_hESC",], aes(x=factor(age), y=percentInTAD, fill=factor(age))) +
        geom_bar(stat="identity") + theme_bw() + theme(text = element_text(size=15)) + scale_fill_manual(values=terrain_hcl(12)) + 
        geom_text(aes(label=paste0(c("n=", rep("", 10)), n), y=-2)) +
        geom_text(aes(label=paste0(signif(percentInTAD, 2), "%")), , vjust=-0.25) + ylab("% in same TAD")
dev.off()

pdf(paste0(outPrefix,".allDF.close.inTAD_by_age.barplot.pdf"), w=3.5,h=7)
    ggplot(d, aes(x=factor(age), y=percentInTAD, fill=factor(age))) +
        geom_bar(stat="identity") + facet_grid(tadSource~.) + theme_bw() + theme(legend.position = "bottom") + theme(text = element_text(size=15)) + scale_fill_manual(values=terrain_hcl(12)) + ylab("% in same TAD")
#~         + geom_text(aes(label=paste0("p=", signif(pval, 3))), data=p.val.df) + ylab("Duplication age")
dev.off()

# dist by age

pdf(paste0(outPrefix,".allDF.close.dist_by_age.boxplot.pdf"), w=7,h=3.5)
    ggplot(subDF[subDF$tadSource=="stable_TADs",], aes(x=factor(age), y=abs(dist), fill=factor(age))) +
        geom_boxplot() + theme_bw() + theme(text = element_text(size=15)) + scale_fill_manual(values=terrain_hcl(12)) + xlab("Duplication Age") + ylab("Genomic Distance [kb]") +
        geom_text(data=data.frame(age=levels(factor(subDF$age)), y=-50, lab=paste0(c("n=", rep("", 10)), d$n[1:11])), aes(label=lab, x=age, y=y))
dev.off()

# exp cor by age

subExpDF <- subset(allDF, group=="paralog" & distGroup=="close" & tadSource=="stable_TADs" & species=="mouse")
subExpDF <- subExpDF[subExpDF$age %in% 0:10,]
require(colorspace)

pdf(paste0(outPrefix,".allDF.close.expCor_by_age.GTEx.boxplot.pdf"), w=7,h=3.5)
    ggplot(subExpDF[subExpDF$expSource=="GTEx",], aes(x=factor(age), y=expCor^2, fill=factor(age))) +
    geom_boxplot() + theme_bw() + theme(text = element_text(size=15)) + scale_fill_manual(values=terrain_hcl(12)) + xlab("Duplication Age") + ylab("Expression Correlation [R^2]")
dev.off()

pdf(paste0(outPrefix,".allDF.close.expCor_by_age.boxplot.pdf"), w=7,h=7)
    ggplot(subExpDF, aes(x=factor(age), y=expCor^2, fill=factor(age))) + facet_grid(expSource~.) +
    geom_boxplot() + theme_bw() + theme(text = element_text(size=15)) + scale_fill_manual(values=terrain_hcl(12)) + xlab("Duplication Age") + ylab("Expression Correlation [R^2]")
dev.off()

# eh shared by age
subEhDF <- subset(allDF, group=="paralog" & distGroup=="close" & tadSource=="stable_TADs"  & expSource=="GTEx" & species=="mouse")
subEhDF <- subEhDF[subEhDF$age %in% 0:10,]
d <- ddply(subDF, .(age), summarize, sharedEh=percentTrue(inTAD), n=length(inTAD))

require(colorspace)


subDF <- subset(allDF, group=="paralog" & distGroup=="close" & expSource=="GTEx" & species=="mouse")

pdf(paste0(outPrefix,".allDF.close.ds_by_TAD.boxplot.pdf"), w=14,h=7)
    ggplot(subDF, aes(x=inTAD, y=hsapiens_paralog_ds, fill=inTAD)) +
        geom_boxplot() + scale_y_log10() +
        facet_grid(~tadSource) + theme_bw() + theme(legend.position = "bottom") + theme(text = element_text(size=15)) + scale_fill_manual(values=rev(COL_TAD))
dev.off()

pdf(paste0(outPrefix,".allDF.close.dn_by_TAD.boxplot.pdf"), w=14,h=7)
    ggplot(subDF, aes(x=inTAD, y=hsapiens_paralog_dn, fill=inTAD)) +
        geom_boxplot() + scale_y_log10() +
        facet_grid(~tadSource) + theme_bw() + theme(legend.position = "bottom") + theme(text = element_text(size=15)) + scale_fill_manual(values=rev(COL_TAD))
dev.off()

subDF <- subset(allDF, distGroup=="close" & expSource=="GTEx" & species=="mouse")
pdf(paste0(outPrefix,".allDF.close.age_commonOrth_by_TAD.barplot.pdf"), w=14,h=7)
    ggplot(subDF, aes(x=ortholog_commonOrtholg, fill=inTAD)) +
        geom_bar(aes(y = ..count..), position="dodge") +
        facet_grid(group~tadSource, scales="free_y") + theme_bw() + theme(legend.position = "bottom") + theme(text = element_text(size=15)) + scale_fill_manual(values=rev(COL_TAD))
dev.off()
    
# same strand vs. distance in paralogs and sampled
subDF <- subset(allDF, distGroup=="close" & tadSource=="stable_TADs" & expSource=="GTEx"  & species=="mouse")
pPara <- wilcox.test(dist ~ sameStrand, data=subDF[subDF$group=="paralog",])
pSamp <- wilcox.test(dist ~ sameStrand, data=subDF[subDF$group=="sampled",])

pVals <- data.frame(strand=rep(1.5,2), dist=1.1*MAX_DIST/10^3, group=c("paralog", "sampled"), pVal=paste0("p = ", signif(c(pPara$p.value, pSamp$p.value),3))) 

pdf(paste0(outPrefix,".allDF.close.dist_by_sameStrand_vs_sampled.boxplot.pdf"), w=3.5, h=7)
    ggplot(subDF, aes(x=strand, y=dist, fill=group)) +
        geom_boxplot() + ylim(0,1.1*MAX_DIST/10^3) + #scale_y_log10() +
        geom_text(data=pVals, aes(x=strand, y=dist, label=pVal)) + 
        facet_grid(~group) + theme_bw() + theme(legend.position = "bottom", axis.text.x=element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=15)) + scale_fill_manual(values=COL) + labs(x="Strand", y="Distance [kb]")
dev.off()

#~ fracSubDF <- ddply(subDF, .(group, ortholog_commonOrtholg, tadSource), summarize, frac=percentTrue(inTAD))
#~ 
#~ sumSubDF <- ddply(subDF, .(group, ortholog_commonOrtholg, tadSource), summarize, frac=sum(inTAD))

# do paralogs exp correlation vs in same TAD
subDF <- subset(allDF, distGroup=="close" & expSource=="GTEx"  & species=="mouse")

pdf(paste0(outPrefix,".allDF.close.expCor_GTEx_vs_inTAD.boxplot.pdf"), w=14, h=7)
    ggplot(subDF, aes(x=inTAD, y=expCor^2, fill=inTAD)) +
        geom_boxplot() +
        facet_grid(group~tadSource) + theme_bw() + theme(legend.position = "bottom", axis.text.x=element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=15)) + scale_fill_manual(values=rev(COL_TAD))
dev.off()

#-----------------------------------------------------------------------
# subset of pairs within same TAD: shared enhancers, expression correlation
#-----------------------------------------------------------------------
# do paralogs exp correlation vs in same TAD
#~ subDF <- subset(allDF, distGroup=="close" & tadSource=="stable_TADs"  & species=="mouse")
subDF <- subset(allDF, distGroup=="close" & tadSource=="Rao_IMR90"  & species=="mouse")
subDF$inTAD <- factor(subDF$inTAD, levels=c(TRUE, FALSE), labels=c("same TAD", "different TADs"))

pdf(paste0(outPrefix,".allDF.close.IMR90.expCor_vs_inTAD.boxplot.pdf"), w=7, h=7)
    ggplot(subDF, aes(x=group, y=expCor, fill=group)) +
        geom_boxplot() +
        facet_grid(inTAD~expSource) + theme_bw() + theme(legend.position = "bottom", axis.text.x=element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=15)) + scale_fill_manual(values=COL)
dev.off()

g <- ggplot(subDF, aes(x=group, y=expCor, fill=group)) +
        geom_violin(aes(color = group), adjust=.5, lwd=1) + geom_boxplot(fill=NA, width=.25, lwd=1) +
        facet_grid(inTAD~expSource) + theme_bw() + theme(legend.position = "none", axis.text.x=element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=14)) + scale_fill_manual(values=COL) + scale_color_manual(values=COL) +
        ylab("Expression Correlation [Pearson R]") + xlab("")
ggsave(paste0(outPrefix,".allDF.close.IMR90.expCor_vs_inTAD.violinplot.pdf"), w=7, h=7)

#-----------------------------------------------------------------------
# shared enhancer
#~ subDF <- subset(ehDF,  tadSource=="stable_TADs" & expSource=="GTEx")
subDF <- subset(ehDF,  tadSource=="Rao_IMR90" & expSource=="GTEx")
subDF$inTAD <- factor(subDF$inTAD, levels=c(TRUE, FALSE), labels=c("same TAD", "different TADs"))
#~ subDF$someSharedEH <- factor(subDF$commonEnhancer >= 1 , levels=c(TRUE, FALSE), labels=c("Shared Enhancer", "No sahred enhancer"))
subDF$someSharedEH <- subDF$commonEnhancer >= 1 

dc <- ddply(subDF, .(group, inTAD, someSharedEH), summarize, count=length(someSharedEH))
summedCounts <- rep(ddply(dc, .(group, inTAD), summarize, s=sum(count))$s, each=2)
subDFcomb <- cbind(dc, cbind(summedCounts, percent=dc$count * 100 / summedCounts))

subDFcombEh = subDFcomb[subDFcomb$someSharedEH == TRUE,]

# fisher test on contingency tables:
ct <- table(subDF[,c("group", "someSharedEH","inTAD")])
pVals <- sapply(levels(subDF$inTAD), function(s) fisher.test(ct[,,s])$p.value )

pValsDF <- data.frame(group=NA, pval=pVals, percent=1.1*max(subDFcombEh$percent, na.rm=TRUE), inTAD=levels(subDF$inTAD))

#~     ggplot(subDF, aes(x=tadSource, y=age, fill=inTAD)) +
#~         geom_boxplot() + theme_bw() + theme(legend.position = "bottom") + theme(text = element_text(size=15)) + scale_fill_manual(values=rev(COL_TAD)) + geom_text(aes(label=paste0("p=", signif(pval, 3))), data=p.val.df) + ylab("Duplication age")

g <- ggplot(subDFcombEh, aes(x=group, fill = group, y = percent)) +
 geom_bar(stat="identity") + facet_grid(.~inTAD) + theme_bw() + theme(legend.position = "bottom", axis.text.x=element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=15)) + scale_fill_manual("", values=COL) + ylab("Gene paris with shared enhancer [%]") + xlab("") +
 geom_text(aes(label = signif(percent,3), y = percent, vjust=-.2), size=5) + 
 geom_text(aes(label=paste0("p=", signif(pval, 3)), x=1.5), data=pValsDF, size=5)
  
ggsave(paste0(outPrefix,".allDF.close.IMR90.sahredEH_vs_inTAD.barplot.pdf"), g, w=3.5, h=7)

#-----------------------------------------------------------------------
# same compartment
#-----------------------------------------------------------------------
subDF <- subset(allDF, tadSource=="stable_TADs" & expSource=="GTEx"  & species=="mouse")
#~ subDF$comp_combination <- factor(subDF$comp_combination, levels=c("A/A", "B/B", "A/B", "unknown"))

#-----------------------------------------------------------------------
# Check how many pairs are in different compartment regions, i.e usable for analysis

subDFdiffRegionFrac <- ddply(subDF, .(group), summarize, 
    diffReg=percentTrue(!is.na(comp_combination) & !same_comp_region), 
    diffRegSum=sum(!is.na(comp_combination) & !same_comp_region), 
    diffSubReg=percentTrue(!is.na(common_subcomp) & !same_subcomp_region),
    diffSubRegSum=sum(!is.na(common_subcomp) & !same_subcomp_region),
    n=length(comp_combination)
    )

g <- ggplot(subDFdiffRegionFrac, aes(x=group, y=diffReg, fill=group)) +
    geom_bar(stat="identity") + 
    theme_bw() + theme(legend.position = "none", axis.text.x=element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=15)) + scale_fill_manual(values=COL) + ylab("Different compartment regions [%]") + xlab("") + 
    geom_text(aes(label=paste0(signif(diffReg, 3),"\n(n=",diffRegSum, ")")), vjust=1.5)
ggsave(paste0(outPrefix,".allDF.diff_comp_region.barplot.pdf"), g, w=3.5, h=7)

g <- ggplot(subDFdiffRegionFrac, aes(x=group, y=diffSubReg, fill=group)) +
    geom_bar(stat="identity") + 
    theme_bw() + theme(legend.position = "none", axis.text.x=element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=15)) + scale_fill_manual(values=COL) + ylab("Different compartment regions [%]") + xlab("") + 
    geom_text(aes(label=paste0(signif(diffSubReg, 3),"\n(n=",diffSubRegSum, ")")), vjust=1.5)
ggsave(paste0(outPrefix,".allDF.diff_subcomp_region.barplot.pdf"), g, w=3.5, h=7)

# separate by distGroup

subDFdiffRegionDistFrac <- ddply(subDF, .(group, distGroup), summarize, diffReg=percentTrue(!is.na(comp_combination) & !same_comp_region), diffSubReg=percentTrue(!is.na(common_subcomp) & !same_subcomp_region))

g <- ggplot(subDFdiffRegionDistFrac, aes(x=group, y=diffReg, fill=group)) +
    geom_bar(stat="identity") + facet_grid(.~distGroup) + 
    theme_bw() + theme(legend.position = "none", axis.text.x=element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=15)) + scale_fill_manual(values=COL) + ylab("Different compartment regions [%]") + xlab("") + 
    geom_text(aes(label=signif(diffReg, 3)), vjust=1.5)
ggsave(paste0(outPrefix,".allDF.diff_comp_region.by_dist.barplot.pdf"), g, w=3.5, h=7)


g <- ggplot(subDFdiffRegionDistFrac, aes(x=group, y=diffSubReg, fill=group)) +
    geom_bar(stat="identity") + facet_grid(.~distGroup) + 
    theme_bw() + theme(legend.position = "none", axis.text.x=element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=15)) + scale_fill_manual(values=COL) + ylab("Different compartment regions [%]") + xlab("") + 
    geom_text(aes(label=signif(diffSubReg, 3)), vjust=1.5)
ggsave(paste0(outPrefix,".allDF.diff_subcomp_region.by_dist.barplot.pdf"), g, w=3.5, h=7)

#-----------------------------------------------------------------------
# Distribution of compartment combinations including NA's

#mark all combinations as NA that have genes in the same compartment region
subDF$comp_combination[subDF$same_comp_region] <- NA
subDF$subcomp_combination[subDF$same_subcomp_region] <- NA

# make factors for sorted plots:
subDF$comp_combination <- factor(subDF$comp_combination, levels=c("A/A", "B/B", "A/B", "unknown"))

dc <- ddply(subDF, .(group, comp_combination), summarize, count=table(as.character(comp_combination), useNA="ifany"))
summedCounts <- rep(ddply(dc, .(group), summarize, s=sum(count))$s, each=length(unique(dc$comp_combination)))
subDFcomb <- cbind(dc, cbind(summedCounts, percent=dc$count * 100 / summedCounts))

# add percent position for labels:
subDFcomb <- ddply(subDFcomb, .(group), 
   transform, pos = cumsum(percent) - (0.5 * percent)
)

# mark pairs with unknown compartment in combination
subDFcomb$comp_combination[is.na(subDFcomb$comp_combination)] <- "unknown"

g <- ggplot(subDFcomb, aes(x=group, fill = comp_combination, y = percent)) +
 geom_bar(stat="identity") + theme_bw() + theme(legend.position = "bottom", axis.text.x=element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=15)) + scale_fill_manual("", values=COL_COMP) + ylab("Compartment combination [%]") + xlab("") +
 geom_text(aes(label = signif(percent,3), y = pos)) 
ggsave(paste0(outPrefix,".allDF.compartment_combination.pdf"), g, w=3.5, h=7)

# same but separated by dist group

dc <- ddply(subDF, .(distGroup, group, comp_combination), summarize, count=table(as.character(comp_combination), useNA="ifany"))
summedCounts <- rep(ddply(dc, .(distGroup, group), summarize, s=sum(count))$s, each=length(unique(dc$comp_combination)))
subDFcomb <- cbind(dc, cbind(summedCounts, percent=dc$count * 100 / summedCounts))

# add percent position for labels:
subDFcomb <- ddply(subDFcomb, .(distGroup, group), 
   transform, pos = cumsum(percent) - (0.5 * percent)
)

# mark pairs with unknown compartment in combination
subDFcomb$comp_combination[is.na(subDFcomb$comp_combination)] <- "unknown"

g <- ggplot(subDFcomb, aes(x=group, fill = comp_combination, y = percent)) +
 geom_bar(stat="identity") + facet_grid(.~distGroup) + theme_bw() + theme(legend.position = "bottom", axis.text.x=element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=15)) + scale_fill_manual("", values=COL_COMP) + ylab("Compartment combination [%]") + xlab("") +
 geom_text(aes(label = signif(percent,3), y = pos)) 
ggsave(paste0(outPrefix,".allDF.compartment_combination.by_dist.pdf"), g, w=3.5, h=7)

#-----------------------------------------------------------------------
# without NAs
dc <- ddply(subDF, .(group, comp_combination), summarize, count=table(as.character(comp_combination)))
summedCounts <- rep(ddply(dc, .(group), summarize, s=sum(count))$s, each=length(unique(dc$comp_combination)))
subDFcomb <- cbind(dc, cbind(summedCounts, percent=dc$count * 100 / summedCounts))

# add percent position for labels:
subDFcomb <- ddply(subDFcomb, .(group), transform, pos = cumsum(percent) - (0.5 * percent))

g <- ggplot(subDFcomb, aes(x=group, fill = comp_combination, y = percent)) +
 geom_bar(stat="identity") + theme_bw() + theme(legend.position = "bottom", axis.text.x=element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=15)) + scale_fill_manual("", values=COL_COMP) + ylab("Compartment combination [%]") + xlab("") +
 geom_text(aes(label = signif(percent,3), y = pos)) 
ggsave(paste0(outPrefix,".allDF.compartment_combination.withoutNA.pdf"), g, w=3.5, h=7)

# same but separated by dist group

dc <- ddply(subDF, .(distGroup, group, comp_combination), summarize, count=table(as.character(comp_combination)))
summedCounts <- rep(ddply(dc, .(distGroup, group), summarize, s=sum(count))$s, each=length(unique(dc$comp_combination)))
subDFcomb <- cbind(dc, cbind(summedCounts, percent=dc$count * 100 / summedCounts))

# add percent position for labels:
subDFcomb <- ddply(subDFcomb, .(distGroup, group), transform, pos = cumsum(percent) - (0.5 * percent))

g <- ggplot(subDFcomb, aes(x=group, fill = comp_combination, y = percent)) +
 geom_bar(stat="identity") + facet_grid(.~distGroup) + theme_bw() + theme(legend.position = "bottom", axis.text.x=element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=15)) + scale_fill_manual("", values=COL_COMP) + ylab("Compartment combination [%]") + xlab("") +
 geom_text(aes(label = signif(percent,3), y = pos)) 
ggsave(paste0(outPrefix,".allDF.compartment_combination.withoutNA.by_dist.pdf"), g, w=3.5, h=7)

#-----------------------------------------------------------------------
# distribution of subcompartment combinations

dc <- ddply(subDF, .(group, subcomp_combination), summarize, count=table(as.character(subcomp_combination), useNA="ifany"))
summedCounts <- rep(ddply(dc, .(group), summarize, s=sum(count))$s, each=length(unique(dc$subcomp_combination)))
subDFcomb <- cbind(dc, cbind(summedCounts, percent=dc$count * 100 / summedCounts))

# add percent position for labels:
subDFcomb <- ddply(subDFcomb, .(group), 
   transform, pos = cumsum(percent) - (0.5 * percent)
)

# mark pairs with unknown compartment in combination
subDFcomb$subcomp_combination[is.na(subDFcomb$subcomp_combination)] <- "unknown"

g <- ggplot(subDFcomb, aes(x=subcomp_combination, fill = group, y = percent)) +
 geom_bar(stat="identity", position="dodge") + theme_bw() + theme(legend.position = "bottom", axis.text.x=element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=15)) + scale_fill_manual("", values=COL) + ylab("%") + xlab("Compartment combination") +
 geom_text(aes(label = signif(percent,3)), position=position_dodge(width=1), vjust=.5, hjust=0, angle=90)
ggsave(paste0(outPrefix,".allDF.compartment_subcombination.pdf"), g, w=14, h=3.5)

# same but without NAs

dc <- ddply(subDF, .(group, subcomp_combination), summarize, count=table(as.character(subcomp_combination)))
summedCounts <- rep(ddply(dc, .(group), summarize, s=sum(count))$s, each=length(unique(dc$subcomp_combination)))
subDFcomb <- cbind(dc, cbind(summedCounts, percent=dc$count * 100 / summedCounts))

# add percent position for labels:
subDFcomb <- ddply(subDFcomb, .(group), transform, pos = cumsum(percent) - (0.5 * percent))

g <- ggplot(subDFcomb, aes(x=subcomp_combination, fill = group, y = percent)) +
 geom_bar(stat="identity", position="dodge") + theme_bw() + theme(legend.position = "bottom", axis.text.x=element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=15)) + scale_fill_manual("", values=COL) + ylab("%") + xlab("Compartment combination") +
 geom_text(aes(label = signif(percent,3)), position=position_dodge(width=1), vjust=.5, hjust=0, angle=90)
ggsave(paste0(outPrefix,".allDF.compartment_subcombination.withoutNA.pdf"), g, w=14, h=3.5)

#-----------------------------------------------------------------------
# Percent of pairs with same compartment

# mark common cmopartment as NA if combination is NA
subDF$common_comp[is.na(subDF$comp_combination)] <- NA
subDF$common_subcomp[is.na(subDF$subcomp_combination)] <- NA

tab <- table(subDF[, c("group", "common_comp")])
allCompP <- fisher.test(tab)$p.value

tabDF <- as.data.frame(tab)
tabDF$percent <- tabDF$Freq / rep(rowSums(tab), 2) * 100


g <- ggplot(tabDF[tabDF$common_comp==TRUE,], aes(x=group, y=percent, fill=group)) +
    geom_bar(stat="identity") + 
    theme_bw() + theme(legend.position = "none", axis.text.x=element_text(angle = 45, hjust = 1)) + theme(text = element_text(size=15), axis.title.x = element_blank()) + scale_fill_manual(values=COL) + ylab("Same Compartment [%]") + ggtitle(paste0("p=", signif(allCompP,3))) +
    geom_text(aes(label=signif(percent, 3)), vjust=1.5)
ggsave(paste0(outPrefix,".allDF.common_compartment.pdf"), g, w=1.75, h=3.5)



#=======================================================================
# save workspace image
#=======================================================================
save.image(WORKIMAGE_FILE)
#INFO: FANTOM enhancer map: 65953 of 66942 enhancer-associated genes could be mapped uniquelly to ENSG ID
