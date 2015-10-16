#!/usr/bin/Rscript
#=======================================================================
#
#   Check the distribution of strand specific genes and pseudo genes
#   with respect to TADs.
# 
#=======================================================================

#-----------------------------------------------------------------------
# load some useful libraries
#-----------------------------------------------------------------------
require(stringr)        # for functions like paste0()
require(GenomicRanges)  # for GRanges objects
require(ggplot2)        # for nice plots
require(xtable)         # to directly output latex tables

#-----------------------------------------------------------------------
# Load parameters from external script
#-----------------------------------------------------------------------
# read parameter script location from command-line argument
#~ args <- commandArgs(trailingOnly = TRUE)
#~ PARAM_SCRIPT=args[1]
PARAM_SCRIPT="R/TAD_gene_distribution.param.v02.R"
source(PARAM_SCRIPT)

#-----------------------------------------------------------------------
# Options for parallel computation
#-----------------------------------------------------------------------
# use all available cores but generate random number streams on each worker
multicorParam <- MulticoreParam(RNGseed=RANDOM_SEED)   
# set options
register(multicorParam)  
# bpparam() # to print current options

# load costume scripts
source("R/functions.GRanges.R")
source("R/functions.plot.R")
source("R/parseHiC.R")

source("R/data.ensembl.R")

#=======================================================================
# parse input data
#=======================================================================

# parse TAD data sets as list of GRanges
allTADs = lapply(RaoDomainFiles, parseDomainsRao, disjoin=FALSE, seqinfo=seqInfo)

TAD = allTADs[["Rao_GM12878"]]
selfOverlap <- countOverlaps(TAD, TAD)


allTADs <- c(allTADs, list("Rao_GM12878_nonOverlapping"=TAD[selfOverlap == 1]))

countOverlaps(TAD, TAD)

parsePsg <- function(inFile){
    
    d = read.table(inFile, header=TRUE)
    gr = GRanges(
        d$Chromosome,
        IRanges(d$Start, d$End),
        strand = d$Strand, 
        names = d$Pseudogene_id,
        d[,c("Biotype", "Parent_gene", "Parent_transcript", "Parent_name")],
        seqinfo=seqInfo
        )    
}

# parse the pseudogenes
PsgList = lapply(PseudoGenesFiles, parsePsg)

#~ regionList = list(
#~     "Non-coding exon" = ncExonGR,
#~     "Coding exon" = codingExonGR,    
#~     "Intron" = intronGR,    
#~     "Intergenic" = intergenicGR
#~ )

# list of GRs of TSS of gene types of interest to be counted inside and outside of the TADs
geneTypeTssList = c(
    list("protGenes"=tssGR), 
    lapply(c(list("lincRNA"=lincRNAsGR), PsgList), resize, width=1, fix="start")
    )

# separate by chrom
geneTypePerChormList <- bplapply(geneTypeTssList,
    function(reg){
        chrList <- GRangesList(lapply(seqlevels(reg), function(chr){
            reg[seqnames(reg) == chr]
        }))
        names(chrList) <- seqlevels(reg)
        return(chrList)
    })

#-----------------------------------------------------------------------
# Transpose a list of n GRAnges with m equal sized bins into list of all the ith bin 
#-----------------------------------------------------------------------
transposeBinList <- function(TADbins, m=20){
    allBins = unlist(TADbins)
    n = length(allBins)
    lapply(1:m, function(i) allBins[seq(i,n,m)])    
}   

fractionOfGenome <- function(gr, seqinfo=seqInfo){
    genomeLength = sum(as.numeric(seqlengths(seqInfo)))
    sum(as.numeric(width(reduce(gr, ignore.strand=TRUE)))) / genomeLength
}

#-----------------------------------------------------------------------
# Summary table of TAD data
#-----------------------------------------------------------------------

TADstats = round(data.frame(
    "n" = sapply(allTADs, length), 
    "mean_size" = sapply(sapply(allTADs, width), mean) / 10^3,
    "sd_size" = sapply(sapply(allTADs, width), sd) / 10^3, 
    "median_size" = sapply(sapply(allTADs, width), median) / 10^3, 
    "percent_genome" = sapply(allTADs, fractionOfGenome),
    "avg_genes" = sapply(allTADs, function(tadGR) mean(countOverlaps(tadGR, tssGR)))
    ), 2)
row.names(TADstats) = gsub("Rao_", "", names(allTADs))

write.table(TADstats, paste0(outPrefix, ".TAD_statistics.csv"), sep="\t", quote=FALSE, col.names=NA)
print.xtable(xtable(TADstats), type="latex", file=paste0(outPrefix, ".TAD_statistics.tex"))

#-----------------------------------------------------------------------
# Summary table of genes
#-----------------------------------------------------------------------
RegStats = round(data.frame(
    "n" = sapply(geneTypeTssList, length), 
    "mean_size" = sapply(sapply(geneTypeTssList, width), mean) / 10^3,
    "sd_size" = sapply(sapply(geneTypeTssList, width), sd) / 10^3, 
    "median_size" = sapply(sapply(geneTypeTssList, width), median) / 10^3, 
    "percent_genome" = sapply(geneTypeTssList, fractionOfGenome)
    ), 2)
    
write.table(RegStats, paste0(outPrefix, ".geneTypeTssList_statistics.csv"), sep="\t", quote=FALSE, col.names=NA)


#=======================================================================
# Divide TADs and surrounding regions into bins and compute coverage and score along TAD body
#=======================================================================

# iterate over bin numbers:

for (nBin in N_BIN_LIST){
    
    # nBin=40

    # Iterate over each set of TADs
    #~ for (TADname in names(allTADs)) {
#~     for (TADname in c("Rao_GM12878")) {
    for (TADname in c("Rao_GM12878_nonOverlapping")) {
        
        #TADname = "Rao_GM12878"
        TAD = allTADs[[TADname]]
        
        # resize TADs by +-50% and divide them into equal sized bins
        TADreg = TAD
        ranges(TADreg) = IRanges(start(TAD) - .5 * width(TAD), end(TAD) + .5 * width(TAD))
        
        # drop TADs that extends out of the chromosome space
        withinChrom = width(TADreg) == width(trim(TADreg))
        TADreg = TADreg[withinChrom]
        
        # make nBin equal sized bins for each TAD
        TADbins = tile(TADreg, n=nBin)
        seqinfo(TADbins) <- seqInfo
        
        # transpose the list to have always the i-th bin together
        TADbinsList = transposeBinList(TADbins, nBin)
        TADbinsSize = width(TADbinsList[[1]])
    
        # get bins arround TAD boundaries
        BoundaryBins = binsAroundTADboundaries(TAD, windowSize=WINDOW_SIZE, nbins=nBin)
        message("INFO: Finished calculation of bins.")
    
        # transform bin list to get each i-th bin together
        boundaryBinList = transposeBinList(BoundaryBins, nBin)
    
        #=======================================================================
        # count regions in Bins of relative size along window around TAD
        #=======================================================================
        
        regInBinSum = rbind()
        regInBinSumPlus = rbind()
        regInBinSumMinus = rbind()
        
        for (typeName in names(geneTypeTssList)){
            
            regType = geneTypeTssList[[typeName]]
            
            chromsOfInterest <- c("chr4", "chr16", "chr19")
            regionOfInterest <- GRangesList(c( 
                list(regType),
                lapply(chromsOfInterest, 
                    function(chr) geneTypePerChormList[[typeName]][[chr]]
                )
            ))
            names(regionOfInterest) <- c(typeName,  paste0(typeName, "_", chromsOfInterest))
            
            for (regName in names(regionOfInterest)){
                
                reg = regionOfInterest[[regName]]
                
                # count for each bin the number of regions overlapping the bin
                countsInBins = lapply(TADbinsList, countOverlaps, reg)
                countsInBinsPlus = lapply(TADbinsList, countOverlaps, reg[strand(reg) == "+"])
                countsInBinsMinus = lapply(TADbinsList, countOverlaps, reg[strand(reg) == "-"])
                
                names(countsInBins) = paste0("bin_", 1:length(countsInBins))
                names(countsInBinsPlus) = paste0("bin_", 1:length(countsInBinsPlus))
                names(countsInBinsMinus) = paste0("bin_", 1:length(countsInBinsMinus))
                                
                # get the mean/sd count over all i-th bins
                sumCount = apply(as.data.frame(countsInBins), 2, sum)
                sumCountPlus = apply(as.data.frame(countsInBinsPlus), 2, sum)
                sumCountMinus = apply(as.data.frame(countsInBinsMinus), 2, sum)
                
                # test total sum of positve/negative strand counts in first and second half of the TAD
                #For 4 bin case: contTab <- rbind(sumCountPlus[2:3], sumCountMinus[2:3])
                firstHalfBins = (1:nBin)[(nBin/4+1):(nBin/2)]
                secondHalfBins = (1:nBin)[(nBin/2+1):(nBin*3/4)]
                contTab <- matrix(c(
                    sum(sumCountPlus[firstHalfBins]),sum(sumCountPlus[secondHalfBins]),
                    sum(sumCountMinus[firstHalfBins]),sum(sumCountMinus[secondHalfBins])
                    ),2, byrow=TRUE)
                fisherTest <- fisher.test(contTab)
                
                # add values to matrix
                regInBinSum = rbind(regInBinSum, sumCount)
                regInBinSumPlus = rbind(regInBinSumPlus, sumCountPlus)
                regInBinSumMinus = rbind(regInBinSumMinus, sumCountMinus)
                
                colReg = COL_PAIRED_2[which(typeName == names(geneTypeTssList))]        
                colRegPlus = COL_PAIRED_2[which(typeName == names(geneTypeTssList))]        
                colRegMinus = COL_PAIRED_1[which(typeName == names(geneTypeTssList))]        
        
                outFile = paste0(outPrefix, ".nBin_", nBin, ".", TADname, ".", regName, ".counts_along_TADs.pdf")
                plotCountsAroundTADwindow(sumCount, outFile, nbins=nBin, col=colReg, main=regName, ylab="Counts per bin")
        
                outFilePlus = paste0(outPrefix, ".nBin_", nBin, ".", TADname, ".", regName, ".pos_counts_along_TADs.pdf")
                plotCountsAroundTADwindow(sumCountPlus, outFilePlus, nbins=nBin, col=colRegPlus, main=regName, ylab="Counts per bin", type="b", pch=-8853)
        
                outFileMinus = paste0(outPrefix, ".nBin_", nBin, ".", TADname, ".", regName, ".neg_counts_along_TADs.pdf")
                plotCountsAroundTADwindow(sumCountMinus, outFileMinus, nbins=nBin, col=colRegMinus, main=paste(regName, "p =", signif(fisherTest$p.value, 3)), ylab="Counts per bin", type="b", pch=-8854)
                
            }
        }
        
        row.names(regInBinSum) = paste0(rep(names(geneTypeTssList), each=length(names(regionOfInterest))), "_", rep(names(regionOfInterest), length(geneTypeTssList)))
        row.names(regInBinSumPlus) = paste0(rep(names(geneTypeTssList), each=length(names(regionOfInterest))), "_", rep(names(regionOfInterest), length(geneTypeTssList)))
        row.names(regInBinSumMinus) = paste0(rep(names(geneTypeTssList), each=length(names(regionOfInterest))), "_", rep(names(regionOfInterest), length(geneTypeTssList)))
    #~     row.names(regInBinSD) = names(geneTypeTssList)
        
        # wirte values to table output
        write.table(regInBinSum, paste0(outPrefix, ".nBin_", nBin, ".", TADname, ".all.sum_counts_along_TADs.csv"), sep="\t", quote=FALSE, col.names=NA)
        write.table(regInBinSumPlus, paste0(outPrefix, ".nBin_", nBin, ".", TADname, ".all.sum_pos_counts_along_TADs.csv"), sep="\t", quote=FALSE, col.names=NA)
        write.table(regInBinSumMinus, paste0(outPrefix, ".nBin_", nBin, ".", TADname, ".all.sum_neg_counts_along_TADs.csv"), sep="\t", quote=FALSE, col.names=NA)
    

        #=======================================================================
        # count regions in Bins around TAD boundaries
        #=======================================================================
        for (regName in names(geneTypeTssList)){
        
            reg = geneTypeTssList[[regName]]

            # count for each bin the number of regions overlapping the bin
            countsInBins = lapply(boundaryBinList, countOverlaps, reg)
            countsInBinsPlus = lapply(boundaryBinList, countOverlaps, reg[strand(reg) == "+"])
            countsInBinsMinus = lapply(boundaryBinList, countOverlaps, reg[strand(reg) == "-"])

            names(countsInBins) = paste0("bin_", 1:length(countsInBins))
            names(countsInBinsPlus) = paste0("bin_", 1:length(countsInBinsPlus))
            names(countsInBinsMinus) = paste0("bin_", 1:length(countsInBinsMinus))
                
            # get the sum/mean/sd count over all i-th bins
            meanCount = apply(as.data.frame(countsInBins), 2, mean)
            meanCountPlus = apply(as.data.frame(countsInBinsPlus), 2, mean)
            meanCountMinus = apply(as.data.frame(countsInBinsMinus), 2, mean)

                
            colReg = COL_PAIRED_2[which(regName == names(geneTypeTssList))]        
            colRegPlus = COL_PAIRED_2[which(regName == names(geneTypeTssList))]        
            colRegMinus = COL_PAIRED_1[which(regName == names(geneTypeTssList))]        
    
            outFile = paste0(outPrefix, ".nBin_", nBin, ".", TADname, ".", regName, ".counts_arround_boundary.pdf")
            plotCountsAroundBoundaries(meanCount, outFile, windowSize=WINDOW_SIZE, nbins=nBin, col=colReg,
            main=regName)

            outFilePlus = paste0(outPrefix, ".nBin_", nBin, ".", TADname, ".", regName, ".pos_counts_arround_boundary.pdf")
            plotCountsAroundBoundaries(meanCountPlus, outFilePlus, windowSize=WINDOW_SIZE, nbins=nBin, col=colRegPlus,
            main=regName, type="b", pch=-8853)

            outFileMinus = paste0(outPrefix, ".nBin_", nBin, ".", TADname, ".", regName, ".neg_counts_arround_boundary.pdf")
            plotCountsAroundBoundaries(meanCountMinus, outFileMinus, windowSize=WINDOW_SIZE, nbins=nBin, col=colRegMinus,
            main=regName, type="b", pch=-8854)

        }
    


    }
}

#=======================================================================
# save workspace image
#=======================================================================
save.image(WORKIMAGE_FILE)


########################################################################
# OLD STUFF:
########################################################################

#~     #-------------------------------------------------------------------
#~     # now do the same for the conservation score
#~     #-------------------------------------------------------------------
#~ 
#~     # compute the mean conservation score for each bin for all TADs
#~     binScore = lapply(TADbins, function(tb) mean(scoreConsReg[tb]))
#~     binScoreDF = data.frame(matrix(unlist(binScore), nrow=length(binScore), byrow=TRUE))
#~     
#~     binScoreMean = colMeans(binScoreDF)
#~     binScoreSE = standardError(binScoreDF)
#~     
#~     # plot it along the TAD body
#~     pdf(paste0(outPrefix, ".", TADname, ".conservation_score_along_TADs.pdf"))
#~     
#~         par(cex=1.5, lwd=2)
#~         plot(binScoreMean, type="n" , xaxt = "n", xlab="", ylab="Average conservation score")
#~     
#~         xlab = c("-50%", "start", "domain", "end", "+50%")
#~         xlabAT = seq(0,N_BINS, N_BINS/4)+.5
#~         axis(1, at=xlabAT, labels=xlab)
#~     
#~         
#~         addWiskers(1:N_BINS, binScoreMean, binScoreSE, length=.05, col="gray")
#~         addPolygon(1:N_BINS, binScoreMean, binScoreSE, col="gray")
#~         
#~         lines(binScoreMean, type="o", pch=20)
#~         
#~         par(xpd=TRUE)
#~         valueRange = max(binScoreMean) - min(binScoreMean)
#~         yMin = min(binScoreMean) - .04 * valueRange
#~         yLow = min(binScoreMean) - .07 * valueRange
#~         polygon(xlabAT[c(2, 2,4, 4)], c(yLow, yMin, yMin, yLow), col="black")
#~         #rasterImage(as.raster(cols_down)[sort(km$cluster)], -31, 0, -29, 1, interpolate=FALSE)
#~     
#~     dev.off()


########################################################################

#~ getPairwiseCor <- function(geneIDs, expDF){
#~     
#~     # get expression matrix with genes as columns and tissues as rows
#~     expMat = t(expDF[geneIDs,])
#~     
#~     # get pearson correlation matrix for all gene paris
#~     m = cor(expMat, method="pearson")
#~ 
#~     # return the upper triangualr matrix (all uniq pairs)
#~     return(m[upper.tri(m)])
#~ }
