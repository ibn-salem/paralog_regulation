########################################################################
#
# A script to analyse the adjacent human gene pairs to correlate synteny 
# breakpoints with TAD boundaries.
#
########################################################################

require(stringr)        # for some string functionality
require(RColorBrewer)   # for nice colors
require(GenomicRanges)  # for genomic intervals and overlap calculation
require(rtracklayer)    # for import.bed
require(plyr)           # count() function
require(ggplot2)        # for nice plots
require(BiocParallel)   # for parallel computing

#-----------------------------------------------------------------------
# Load parameters from external script
#-----------------------------------------------------------------------
# read parameter script location from command-line argument
#~ args <- commandArgs(trailingOnly = TRUE)
#~ PARAM_SCRIPT=args[1]

PARAM_SCRIPT="R/adjacent_genes.param.v02.R"
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
source("R/parseHiC.R")
source("R/functions.Hi-C.R")
source("R/functions.GRanges.R")

# load ensemble data sets of genes and paralog pairs
source("R/data.ensembl.R")     # load ensambl data


#-----------------------------------------------------------------------
# Parse TADs
#-----------------------------------------------------------------------

# parse TAD data sets as list of GRanges
RaoTADs = lapply(RaoDomainFiles, parseDomainsRao, disjoin=FALSE, seqinfo=seqInfo)
DixonTADs <- lapply(DixonDomainFiles, import.bed, seqinfo=seqInfo)

stableTADs <- list(
#~     "stable_TADs_n3_f80"=getConservedTADs(RaoTADs, n=3, fraction=.8),
#~     "stable_TADs_n3_f90"=getConservedTADs(RaoTADs, n=3, fraction=.9),
#~     "stable_TADs_n4_f80"=getConservedTADs(RaoTADs, n=4, fraction=.8),
    "stable_TADs"=getConservedTADs(RaoTADs, n=4, fraction=.9)
)

#~ stableTADsRes <- list(
#~     "stable_TADs_n3_10kb"=getConservedByHits(resolutionOverlap(allTADsRaoGR, allTADsRaoGR, resolution=10^4), n=3),
#~     "stable_TADs_n3_50kb"=getConservedByHits(resolutionOverlap(allTADsRaoGR, allTADsRaoGR, resolution=5*10^4), n=3),
#~     "stable_TADs_n4_10kb"=getConservedByHits(resolutionOverlap(allTADsRaoGR, allTADsRaoGR, resolution=10^4), n=4),
#~     "stable_TADs_n4_50kb"=getConservedByHits(resolutionOverlap(allTADsRaoGR, allTADsRaoGR, resolution=10^4), n=4)
#~ )

allTADs = c(
    RaoTADs,
    DixonTADs,
    stableTADs
#~     stableTADsRes
)

message("INFO: Finshed parsing of TADs.")


orgStr2Name = c(mmusculus="mouse", cfamiliaris="dog")

# parse mouse and dog TADs
speciesSeqInfo = list("mmusculus"=seqInfoMouse, "cfamiliaris"=seqInfoDog)


#~ speciesTADs = lapply(1:2, function(i)
#~     parseRudanTADs(RudanFile, sheet=i, disjoin=FALSE, seqinfo=speciesSeqInfo[[i]])
#~     )
#~ names(speciesTADs) = names(speciesSeqInfo)


#-----------------------------------------------------------------------
# adjacent gene paris analysis for correlation of synteny breaks and TAD borders
# - annotate adjacent pairs with
#   - orthologs same chrom
#   - orthologs adjacent
#   - same TAD for all TADs
# - test if linear distance between adjacent paris is different for those
#   that are also adjacent (sameChrom) in mouse/dog
#-----------------------------------------------------------------------

# get all adjacent pairs in the human genome:
adjPairs = getAdjacentPairs(tssGR)

# add distance between adjacent genes in mouse and human
adjPairs = addPairDist(adjPairs, tssGR)

# annotate pairs with orthologs in mouse and dog genome
adjPairs = addOrthologAnnotation(adjPairs, orthologsSpeciesList[["mmusculus"]], "mmusculus", tssGRmouse)
adjPairs = addOrthologAnnotation(adjPairs, orthologsSpeciesList[["cfamiliaris"]], "cfamiliaris", tssGRdog)

# check if the one to one orthologs are also adjacent
adjPairs[,"mmusculus_orthologs_adjacent"] = orthologsAdjacent(adjPairs, "mmusculus", orthologsSpeciesList, tssGRmouse)
adjPairs[,"cfamiliaris_orthologs_adjacent"] = orthologsAdjacent(adjPairs, "cfamiliaris", orthologsSpeciesList, tssGRdog)

# convert to GR object
adjPairsGR = getPairAsGR(adjPairs, tssGR)

# annotate gene pairs if they overlap TADs
for(tadName in names(allTADs)){
    TAD = allTADs[[tadName]]
    # co-occurance within the same domain
    adjPairs[,paste0(tadName, "_TAD")] = getWithinSubject(adjPairsGR, TAD)
}

# make one data fram with all numbers for ggplot
nSpecies <- length(speciesSeqInfo)
nPairs <- nrow(adjPairs)
nTAD <- length(allTADs)

#equalBreaks = seq(10^0, 10^6, length.out=11) / 10^3
#quantileBreaks = quantile(adjPairs[,"dist"],  seq(0, 1, .25)) / 10^3
quants = quantile(adjPairs[,"dist"])
quantileBreaks = cut(adjPairs[,"dist"], breaks=quants, labels=names(quants)[2:length(quants)], include.lowest=TRUE)

adjPlotDF <- data.frame(
    dist = rep(adjPairs[,"dist"], nSpecies * nTAD),
    quantileBin=rep(quantileBreaks, nSpecies * nTAD),
    adj = unlist(lapply(names(speciesSeqInfo), function(s){rep(adjPairs[, paste0(s, "_orthologs_adjacent")], nTAD)})),
    inTAD = rep(unlist(lapply(names(allTADs), function(tad) adjPairs[, paste0(tad, "_TAD")])), nSpecies),
    tissue = rep(rep(names(allTADs), each=nPairs), nSpecies),
    species = rep(orgStr2Name, each=nPairs*nTAD)
)

# order factors by name:
adjPlotDF$adj <- factor(adjPlotDF$adj, levels=c(TRUE, FALSE), labels=c("syntenic", "rearranged"))
adjPlotDF$inTAD <- factor(adjPlotDF$inTAD, levels=c(TRUE, FALSE), labels=c("Same TAD", "Not same TAD"))

# count occurrences of combinations
adjPlotDFcout <- count(adjPlotDF[!is.na(adjPlotDF[,"adj"]),], var=c("adj", "inTAD", "tissue", "species", "quantileBin"))

# get fisher test p-values and odds ratios:
# iterate over pairs of conditions:
pVals <- sapply(orgStr2Name, function(org){
    sapply(names(allTADs), function(tad){
        fisher.test(table(
            adjPlotDF[adjPlotDF$species == org & adjPlotDF$tissue == tad, c("inTAD", "adj")]
            ))$p.value
    })
})

ors <- sapply(orgStr2Name, function(org){
    sapply(names(allTADs), function(tad){
        fisher.test(table(
            adjPlotDF[adjPlotDF$species == org & adjPlotDF$tissue == tad, c("inTAD", "adj")]
            ))$estimate
    })
})

yMax = max(adjPlotDFcout$freq)

# make data frame for plot annotations
annotDF <- data.frame(
    freq=rep(yMax, nTAD * nSpecies),
    inTAD=rep(factor(FALSE, levels=c(TRUE, FALSE), labels=c("Same TAD", "Not same TAD")), nTAD * nSpecies),
    adj=rep(factor(TRUE, levels=c(TRUE, FALSE), labels=c("syntenic", "rearranged")), nTAD*nSpecies),
    tissue=rep(names(allTADs), nSpecies),
    species=rep(orgStr2Name, each=nTAD),
    pVals=as.vector(pVals),
    ors=as.vector(ors)
)


#-----------------------------------------------------------------------
# Distance between adjacent gene pairs
#-----------------------------------------------------------------------

pdf(paste0(outPrefix, ".adjacent_genes.all.distance_vs_sameTAD.boxplot.pdf"), w=14, h=7)
    ggplot(adjPlotDF[!is.na(adjPlotDF[,"adj"]),], na.rm=TRUE, aes(x=adj, y = log10(dist), fill=inTAD)) + geom_boxplot() + theme_bw() + facet_grid(species~tissue) + scale_fill_manual(values=COL_TAD) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    labs(y="Distance between adjacent genes [log10(bp)]", x="")    
dev.off()

pdf(paste0(outPrefix, ".adjacent_genes.all.distance_vs_sameTAD.byDistQuant.boxplot.pdf"), w=14, h=7)
    ggplot(adjPlotDF[!is.na(adjPlotDF[,"adj"]),], na.rm=TRUE, aes(x=adj, y = log10(dist), fill=inTAD)) + geom_boxplot() + theme_bw() + facet_grid(species+quantileBin~tissue, scales="free_y") + scale_fill_manual(values=COL_TAD) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    labs(y="Distance between adjacent genes [log10(bp)]", x="")    
dev.off()

#-----------------------------------------------------------------------
# Number of rearranged pairs 
#-----------------------------------------------------------------------
pdf(paste0(outPrefix, ".adjacent_genes.all.adjacent_vs_sameTAD.barplot.pdf"), w=14, h=7)    
    ggplot(adjPlotDFcout, aes(x=adj, y=freq, fill=inTAD)) + geom_bar(position="dodge",stat="identity") +  
    facet_grid(species~tissue) +
    theme_bw() + scale_fill_manual(values=COL_TAD) + ylim(0, 1.2*yMax) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_text(aes(label=freq), position=position_dodge(width=1), vjust=-0.25, size=3) + 
    geom_text(data=annotDF, aes(x=1.5, y=freq+.15*yMax, label=paste0("OR=", signif(ors, 2))), size=3 ) +
    geom_text(data=annotDF, aes(x=1.5, y=freq+.1*yMax, label=paste0("p=", signif(pVals, 2))), size=3) + 
    labs(y="Adjacent gene pairs", x="")
dev.off()

pdf(paste0(outPrefix, ".adjacent_genes.all.adjacent_vs_sameTAD.byDistQuant.barplot.pdf"), w=14, h=7)

    ggplot(adjPlotDFcout, aes(x=adj, y=freq, fill=inTAD)) + geom_bar(position="dodge",stat="identity") +  
    facet_grid(species+quantileBin~tissue) +
    theme_bw() + scale_fill_manual(values=COL_TAD) 

dev.off()
    
#~     + ylim(0, 1.2*yMax) +
#~     theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#~     geom_text(aes(label=freq), position=position_dodge(width=1), vjust=-0.25, size=3) + 
#~     geom_text(data=annotDF, aes(x=1.5, y=freq+.15*yMax, label=paste0("OR=", signif(ors, 2))), size=3 ) +
#~     geom_text(data=annotDF, aes(x=1.5, y=freq+.1*yMax, label=paste0("p=", signif(pVals, 2))), size=3) + 
#~     labs(y="Adjacent gene pairs", x="")


#=======================================================================
# save workspace
#=======================================================================
save.image(WORKIMAGE_FILE)

#=======================================================================
# OLD version with iteration over all combinations of species and TADs
#=======================================================================
# iterate over species
pValsChromSpecies = c()
oddsRatioChromSpecies = c()
pValsAdjacentSpecies = c()
oddsRatioAdjacentSpecies = c()

for (orgStr in names(speciesSeqInfo)){
    
    orgName = orgStr2Name[orgStr]
    
    # get subset vector for bothe genes having orthologs
    bothOrtholog <- adjPairs[,paste0(orgStr, "_one2one")] & !is.na(adjPairs[,paste0(orgStr, "_orthologs_adjacent")])
    sameChromSpecies = adjPairs[bothOrtholog, paste0(orgStr, "_sameChrom")]    
    orthologAdjacent = adjPairs[bothOrtholog, paste0(orgStr, "_orthologs_adjacent")]
    
    message(orgName)
    message(paste("INFO: bothOrtholog NA?:", sum(is.na(bothOrtholog))))
    message(paste("INFO: sameChromSpecies NA?:", sum(is.na(sameChromSpecies))))
    message(paste("INFO: orthologAdjacent NA?:", sum(is.na(orthologAdjacent))))
    
    pValsSameChrom = c()
    oddsRatioSameChrom = c()
    pValsAdjacent = c()
    oddsRatioAdjacent = c()
    
    # iterate over TAD types and test for correlation
    for(tadName in names(allTADs)){
        
        # check for difference in lienar distance
        sameTAD <- adjPairs[bothOrtholog, paste0(tadName, "_TAD")]
        distSub <- abs(adjPairs[bothOrtholog, "dist"])
        pWilcox = wilcox.test(distSub[sameTAD], distSub[!sameTAD])$p.value
        
        pdf(paste0(outPrefix, ".adjacent_genes.", orgName, ".sameTAD_", tadName, "dist.boxplot.pdf"))
            boxplot(log10(distSub[sameTAD]), log10(distSub[!sameTAD]) )
        dev.off()
        
        # freq matrix
        freqMatChrom = matrix(c(
            percentTrue(sameChromSpecies & sameTAD),
            percentTrue(sameChromSpecies & !sameTAD),
            percentTrue(!sameChromSpecies & sameTAD),
            percentTrue(!sameChromSpecies & !sameTAD)
            ), 2)
        chiChrom = chisq.test(sameChromSpecies, sameTAD)
        fisherChrom = fisher.test(sameChromSpecies, sameTAD)

        pdf(paste0(outPrefix, ".adjacent_genes.sameChrom_", orgName, ".sameTAD_", tadName, ".pdf"))
            par(cex=2, lwd=1.5)
            yMax = 1.75*max(freqMatChrom)
            bp = my.barplot(freqMatChrom, beside=TRUE, yMax=yMax, col=COL_TAD, 
                ylab="Adjacent human gene pairs [%]", names=paste0(c("Same chr in", "Not same chr in"), "\n", orgName), main=paste("P-value =", signif(chiChrom$p.value, 3)), addValues=TRUE)
            legend("topright", c("Same TAD", "Not same TAD"), fill=COL_TAD)
        dev.off()

        # freq matrix
        freqMatAdj = matrix(c(
            sum(orthologAdjacent & sameTAD, na.rm=TRUE),
            sum(orthologAdjacent & !sameTAD, na.rm=TRUE),
            sum(!orthologAdjacent & sameTAD, na.rm=TRUE),
            sum(!orthologAdjacent & !sameTAD, na.rm=TRUE)
            ), 2)
        chiAdj = chisq.test(orthologAdjacent, sameTAD)
        fisherAdj = fisher.test(orthologAdjacent, sameTAD)
        pdf(paste0(outPrefix, ".adjacent_genes.adjacent_", orgName, ".sameTAD_", tadName, ".pdf"))
            par(cex=2, lwd=1.5)
            yMax = 1.75*max(freqMatAdj)
            bp = my.barplot(freqMatAdj, beside=TRUE, col=COL_TAD, 
                ylab="Adjacent human gene pairs [%]", names=paste0(c("Adjacent in", "Not adjacent in"), "\n", orgName), main=paste("P-value =", signif(chiAdj$p.value, 3)), addValues=TRUE, yMax=yMax)
            legend("topright", c("Same TAD", "Not same TAD"), fill=COL_TAD)
        dev.off()
        
        # add p-values and odds ratio to list
        pValsSameChrom = c(pValsSameChrom, chiChrom$p.value)
        oddsRatioSameChrom = c(oddsRatioSameChrom, fisherChrom$estimate)
        pValsAdjacent = c(pValsAdjacent, chiAdj$p.value)
        oddsRatioAdjacent = c(oddsRatioAdjacent, fisherAdj$estimate)
    }
    
    pValsChromSpecies = rbind(pValsChromSpecies, pValsSameChrom)
    oddsRatioChromSpecies = rbind(oddsRatioChromSpecies, oddsRatioSameChrom)
    pValsAdjacentSpecies = rbind(pValsAdjacentSpecies, pValsAdjacent)
    oddsRatioAdjacentSpecies = rbind(oddsRatioAdjacentSpecies, oddsRatioAdjacent)
    
}

pdf(paste0(outPrefix, ".adjacent_genes.sameChrom_vs_TAD.oddsRatio.barplot.pdf"))
    par(cex=1.5, lwd=1.5)
    bp = my.barplot(oddsRatioChromSpecies, beside=TRUE, ylab="Odds ratio", names=names(allTADs), col=COL_SPECIES)
    abline(h=1, col="red")
    legend("topleft", orgStr2Name[names(speciesSeqInfo)], fill=COL_SPECIES)
dev.off()

pdf(paste0(outPrefix, ".adjacent_genes.sameChrom_vs_TAD.pValues.barplot.pdf"))
    par(cex=1.5, lwd=1.5)
    starsNotation = matrix(pValToStars(pValsChromSpecies, ""), nrow(pValsChromSpecies))
    
    bp = my.barplot(-log10(pValsChromSpecies), offset=0.05, beside=TRUE, ylab="-log_10 P-Value", names=names(allTADs), col=COL_SPECIES, addValues=TRUE, customValues=starsNotation, digits=0, main="Re-arrengments versus TAD")
    legend("topleft", orgStr2Name[names(speciesSeqInfo)], fill=COL_SPECIES)
dev.off()


pdf(paste0(outPrefix, ".adjacent_genes.adjacentOrtho_vs_TAD.oddsRatio.barplot.pdf"))
    par(cex=1.5, lwd=1.5)
    bp = my.barplot(oddsRatioAdjacentSpecies, beside=TRUE, ylab="Odds ratio", names=names(allTADs), col=COL_SPECIES)
    abline(h=1, col="red")
    legend("topleft", orgStr2Name[names(speciesSeqInfo)], fill=COL_SPECIES)
dev.off()

pdf(paste0(outPrefix, ".adjacent_genes.adjacentOrtho_vs_TAD.pValues.barplot.pdf"))
    par(cex=1.5, lwd=1.5)
    starsNotation = matrix(pValToStars(pValsAdjacentSpecies, ""), nrow(pValsAdjacentSpecies))
    
    bp = my.barplot(-log10(pValsAdjacentSpecies), offset=0.05, beside=TRUE, ylab="-log_10 P-Value", names=names(allTADs), col=COL_SPECIES, addValues=TRUE, customValues=starsNotation, main="Conserved neighbour genes vs. same TAD")
    legend("topleft", orgStr2Name[names(speciesSeqInfo)], fill=COL_SPECIES)
dev.off()


