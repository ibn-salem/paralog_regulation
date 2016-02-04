########################################################################
#
# Script to work on a fair sampling of gene pairs
#
# Try out the following to fix sampling problems:
#   X remove gene pairs closer than 1kb
#   X remove gene overlapping gene pairs
#   - remove gene pairs closer than 10kb if on same strand
#   X sample by distance in log space
#   X remove Hi-C contact of genes located on the same bin
#   - remove paralog pairs form smapled pairs to have negative-control
#
########################################################################


require(biomaRt)        # to retrieve human paralogs from Ensembl
require(stringr)        # for some string functionality
require(RColorBrewer)   # for nice colors
require(GenomicRanges)  # for genomic intervals and overlap calculation
require(rtracklayer)    # for import.bed
require(plyr)           # count() function
require(gplots)         # heatmap.2 function
require(data.table)     # for data.table object
require(ggplot2)        # for nice plots
require(scales)         # for proper logarithmic scales in ggplot
require(BiocParallel)   # for parallel computing

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#-----------------------------------------------------------------------
# External data file paths
#-----------------------------------------------------------------------
REGMAP_FILE_FANTOM5 = "data/Andersson2014/enhancer_tss_associations.bed"
EH_FILE_FANTOM5 = "data/Andersson2014/permissive_enhancers.bed"


RaoDomainFiles = c(
    Rao_HeLa="data/Rao2014/GSE63525_HeLa_Arrowhead_domainlist.txt",
    Rao_HUVEC="data/Rao2014/GSE63525_HUVEC_Arrowhead_domainlist.txt",
    Rao_K562="data/Rao2014/GSE63525_K562_Arrowhead_domainlist.txt",
    Rao_KBM7="data/Rao2014/GSE63525_KBM7_Arrowhead_domainlist.txt",
    Rao_NHEK="data/Rao2014/GSE63525_NHEK_Arrowhead_domainlist.txt",
    Rao_IMR90="data/Rao2014/GSE63525_IMR90_Arrowhead_domainlist.txt",
    Rao_GM12878="data/Rao2014/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt"
)

DixonDomainFiles = c(
    Dixon_hESC="data/Dixon2012/hESC.hg18.bed.hg19.bed",
    Dixon_IMR90="data/Dixon2012/IMR90.hg18.bed.hg19.bed"
)

# Mouse and dog TADs from Rudan et al. 2015
RudanFile = "data/Rudan2015/mmc2.xlsx"

#-----------------------------------------------------------------------
# Rao et al. 2014 Hi-C sample
#-----------------------------------------------------------------------
CELL="IMR90"
HIC_RESOLUTION=50*10^3 # 50kb
HIC_DATA_DIR="data/Rao2014"

#-----------------------------------------------------------------------
# Colors used for plotting
#-----------------------------------------------------------------------
# two colors for paralog and non-paralog genes
COL=brewer.pal(9, "Set1")[c(1,9)]   # for paralog vs. sampled genes
COL_RAND=c(COL[1], brewer.pal(8, "Dark2")[8])   # random genes
COL2=brewer.pal(9, "Set1")[c(1,2)]  # for paralog vs. non-paralog
COL_DOMAIN=brewer.pal(12, "Set3")
COL_FAMILY=brewer.pal(8, "Dark2")
COL_EH = brewer.pal(9, "Set1")[5]
COL_TAD = brewer.pal(8, "Set1")[c(3,5)]
COL_ORTHO = c(brewer.pal(12, "Paired")[5], brewer.pal(8, "Pastel2")[8])
COL_SPECIES = brewer.pal(8, "Accent")[1:2]
COL_STRAND = brewer.pal(8, "Pastel1")[2:1]
COL_AGE = brewer.pal(9, "YlOrRd")[c(6,9)]
COL_EH_POS=brewer.pal(9, "Set3")[c(5,9,6)]

#-----------------------------------------------------------------------
# Parameters for data loading
#-----------------------------------------------------------------------

# use local data or downlaod data from ensemble
USE_LOCAL = TRUE

#-----------------------------------------------------------------------
# Parameters critical to the analysis
#-----------------------------------------------------------------------

# number of random permutations of whole data sets
N_RAND=10               

# parameter to adjust bandwidth in density estimation for sampling
DENSITY_BW_ADJUST=0.1   

RANDOM_SEED=13521

MAX_DIST=10^6
DISTAL_MAX_DIST=10^9
DISTAL_MIN_DIST=10^6

DIST_TH=c(10^6, 10^5)

#-----------------------------------------------------------------------
# Version of the analysis and output file paths
#-----------------------------------------------------------------------

VERSION="vSamp"

outDataPrefix = "results/paralog_regulation/EnsemblGRCh37_paralog_genes"
outPrefix = paste0("results/paralog_regulation/", VERSION, "_maxDist_", MAX_DIST, "_nrand_", N_RAND, "/EnsemblGRCh37_paralog_genes")

# create directory, if not exist
dir.create(dirname(outPrefix), recursive=TRUE, showWarnings = FALSE)

# image file to save temporary and final R session with all data
WORKIMAGE_FILE = paste0(outPrefix, "workspace.Rdata")

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
source("R/functions.regMap.R")
source("R/functions.GRanges.R")
source("R/functions.genePairs.R")
source("R/functions.genePairs.randomization.R")
source("R/functions.genePairs.paralog_analysis.R")
source("R/functions.Hi-C.R")
source("R/parseHiC.R")

# load ensemble data sets of genes and paralog pairs
source("R/data.ensembl.R")     # load ensambl data
source("R/data.expression.R")  # load expression data from EBI expression atlas
source("R/data.captureHiC.R")  # load capture Hi-C data between promoters from Mifsud et al. 2015

#=======================================================================
# 1.) Parse data
#=======================================================================

#-----------------------------------------------------------------------
# Parse enhancers in correlation based maps
#-----------------------------------------------------------------------
ehGR = rtracklayer::import.bed(EH_FILE_FANTOM5, seqinfo=seqInfoRealChrom)

# parse regulatory map from FANTOM5
regMap = parseFANTOM5RegMap(REGMAP_FILE_FANTOM5, tssGR, ehGR, refSeqToENSG)

# add distance between enhancer and gene to map
regMap$dist = getMapDist(regMap, tssGR, ehGR)

# make GR for associations
mapGR = getMapAsGR(regMap, tssGR, ehGR, strand.as.direction=TRUE)

# map gene symbols to linked enhancer IDs 
gene2ehID = getGenetoEhIDmapping(names(tssGR)[regMap[,1]], regMap[,2])

# add number of linked enhancer to tssGR
tssGR$linked_enhancer =  sapply(gene2ehID[names(tssGR)], length)

#-----------------------------------------------------------------------
# Load Hi-C data from Rao et al. 2014
#-----------------------------------------------------------------------
if ( !USE_LOCAL) {

    HiClist = parseRaoHiC(CELL, HIC_RESOLUTION, HIC_DATA_DIR, seqInfo)
    
    # normalize by expected counts based on linear distance
    HiClistNorm = parseRaoHiC(CELL, HIC_RESOLUTION, HIC_DATA_DIR, seqInfo, normalizeByExpected=TRUE)
    
    # save or load downloaded data 
    save(HiClist, HiClistNorm, file=paste0(outDataPrefix, ".Hi-C.", CELL, ".", HIC_RESOLUTION, ".RData"))
    
    # write Hi-C interaction in track.txt format for visuallization in Epi-Genome browser
#~     writeHiCinteractions(HiClist[1], paste0(outDataPrefix, ".Hi-C.", CELL, ".", HIC_RESOLUTION, ".map.track.txt"))

}else{
    load(paste0(outDataPrefix, ".Hi-C.", CELL, ".", HIC_RESOLUTION, ".RData"))
}
#=======================================================================
# 2.) Annotate and filter paralog gene pairs
#=======================================================================

# FILTERING ORDER:
# paralogPairs                      w dups        n=
# +-- paralogPairsUniqP             w/o dups      n=
# +-- paralogPairsUniqG             w dups        n=
#     +-- paralogPairsUniqGuniqP
#         +-- allCisPairs           
#             +-- closePairs           
#             +-- distalPairs           


# the data.frame "paralogPairs" is loaded from data.ensembl.R

# remove double entries of the form A-B and B-A
paralogPairsUniqP = uniquePair(paralogPairs)

# get for each gene only one unique pair, the one with highest similarity
# this is computed by an maximum weight matching
paralogPairsWithDS = paralogPairs[!is.na(paralogPairs[,"hsapiens_paralog_ds"]),]
paralogPairsUniqG = uniquePairPerGeneBySim(paralogPairsWithDS, -1*paralogPairsWithDS[,"hsapiens_paralog_ds"])
    
# get only a unique pair order (one of A-B, B-A) form the unique pairs
paralogPairsUniqGuniqP = uniquePair(paralogPairsUniqG)

# subset of paralog pairs that are located on the same chromosome
allCisPairs = getCisPairs(paralogPairsUniqGuniqP, tssGR)


#-----------------------------------------------------------------------
# Annotate:
#-----------------------------------------------------------------------
# add same Strand info to all paralogs
paralogPairsUniqGuniqP = addSameStrand(paralogPairsUniqGuniqP, tssGR)

# add linear distance between TSS
allCisPairs = addPairDist(allCisPairs, tssGR)

#~ # add number of common enhancers
#~ allCisPairs = addCommonEnhancer(allCisPairs, gene2ehID)
#~ 
#~ # add position of enhancers:
#~ allCisPairs = addRelativeEnhancerPosition(allCisPairs, tssGR, gene2ehID, ehGR)


# Adds Hi-C contact frequencies to a gene pair data set
allCisPairs <- addHiCfreq(allCisPairs, tssGR, HiClist, ignoreSameBin=TRUE)
allCisPairs <- addHiCfreq(allCisPairs, tssGR, HiClistNorm, label="HiCnorm", ignoreSameBin=TRUE)

#-----------------------------------------------------------------------
# Filter for close and distal pairs
#-----------------------------------------------------------------------

# get close cis pairs
closePairsOrg = allCisPairs[abs(allCisPairs$dist) <= MAX_DIST & abs(allCisPairs$dist) > 0,] 
closePairs = allCisPairs[abs(allCisPairs$dist) <= MAX_DIST & abs(allCisPairs$dist) > 10^3,] 

closeNonOVL <- nonOverlappingGenePairs(closePairsOrg, genesGR)
closePairsNonOVL <- closePairsOrg[closeNonOVL,]

# get distal pairs
distalPairs = allCisPairs[abs(allCisPairs$dist) > MAX_DIST,] 

# pairs with zero distance (artefact of annotation?)
zereDistPairs = allCisPairs[abs(allCisPairs$dist) == 0,]


#-----------------------------------------------------------------------
# paralog pair filtering numbers
#-----------------------------------------------------------------------
nPairs = c(
    "paralogPairs"=nrow(paralogPairs), 
    "paralogPairsUniqP"=nrow(paralogPairsUniqP), 
    "paralogPairsWithDS"=nrow(paralogPairsWithDS), 
    "paralogPairsUniqG"=nrow(paralogPairsUniqG), 
    "paralogPairsUniqGuniqP"=nrow(paralogPairsUniqGuniqP),
    "allCisPairs"=nrow(allCisPairs),
    "zereDistPairs"=nrow(zereDistPairs),
    "closePairsOrg"=nrow(closePairsOrg),
    "closePairs"=nrow(closePairs),
    "closePairsNonOVL"=nrow(closePairsNonOVL),
    "distalPairs"=nrow(distalPairs)
    )
    
write.table(nPairs, file=paste0(outPrefix, ".paralog_pairs_filtering.txt"),
    sep="\t", quote=FALSE, col.names=FALSE)

# save before sampling
#~ save.image(paste0(WORKIMAGE_FILE, ".only_annotation.Rdata"))

#=======================================================================
# 3.) Sample random control/background data sets
#=======================================================================
# There will be the following types of sampled gene pairs:
# - randPairs           := sampled completely pairs randomly from all genes
# - sampClosePairs      := sampled by distance of closePairs
# - sampEhClosePairs    := sampled by distance and num. of enhancers from close Pairs
# - sampDistalPairs     := sampled by distance of distalPairs


#-----------------------------------------------------------------------
# Sample cis pairs according to linear distance in paralog gene pairs
#-----------------------------------------------------------------------
# get all possible gene pairs within MAX_DIST bp
allCloseGenePairsOrg <- getAllGenePairs(tssGR, maxDist=MAX_DIST, minDist=1)

# get sample weights according to distance
closeWeightsOrg <- getSampleWeightsByDist(allCloseGenePairsOrg, sourcePairs=closePairsOrg, adjust=DISTAL_MIN_DIST)

# sample close pairs according to distance weight
sampCloseOrgPairs <- bplapply(1:N_RAND, function(x){ 
    sampleFromAllPairsByWeight(n=nrow(closePairs), hitDF=allCloseGenePairs, tssGR, weight=closeWeights)
    })
Sys.sleep(1) # hack to fix problems with bplapply on MOGON

#-----------------------------------------------------------------------
# Sample cis pairs according to linear distance in paralog gene pairs
#-----------------------------------------------------------------------
# get all possible gene pairs within MAX_DIST bp
allCloseGenePairs <- getAllGenePairs(tssGR, maxDist=MAX_DIST, minDist=10^3)

# get sample weights according to distance
#~ closeWeights <- getSampleWeightsByDist(allCloseGenePairs, sourcePairs=closePairs, adjust=DENSITY_BW_ADJUST)
#~ closeWeights <- getSampleWeightsByDistLog(allCloseGenePairs, sourcePairs=closePairs, adjust=1)
#~ closeWeights <- weightsByBin(log10(abs(closePairs$dist)), log10(abs(allCloseGenePairs$dist)), breaks=50)
N_SAMPLING_BIN <- 20
breaks = seq(log10(10^3), log10(MAX_DIST), length.out=N_SAMPLING_BIN+1)
closeWeights <- weightsByBin(log10(abs(closePairs$dist)), log10(abs(allCloseGenePairs$dist)), breaks=breaks)

# sample close pairs according to distance weight
sampClosePairs <- bplapply(1:N_RAND, function(x){ 
    sampleFromAllPairsByWeight(n=nrow(closePairs), hitDF=allCloseGenePairs, tssGR, weight=closeWeights)
    })
Sys.sleep(1) # hack to fix problems with bplapply on MOGON

#-----------------------------------------------------------------------
# Non-overlapping close pairs
#-----------------------------------------------------------------------
# get all possible gene pairs within MAX_DIST bp
allNonOVL <- nonOverlappingGenePairs(allCloseGenePairsOrg, genesGR, useIDs=TRUE)
allCloseGenePairsNonOVL <- allCloseGenePairsOrg[allNonOVL,]

# get sample weights according to distance
#~ closeWeights <- getSampleWeightsByDist(allCloseGenePairs, sourcePairs=closePairs, adjust=DENSITY_BW_ADJUST)
#~ closeWeights <- getSampleWeightsByDistLog(allCloseGenePairs, sourcePairs=closePairs, adjust=1)
#~ closeWeights <- weightsByBin(log10(abs(closePairs$dist)), log10(abs(allCloseGenePairs$dist)), breaks=50)
N_SAMPLING_BIN <- 20
breaks = seq(log10(1), log10(MAX_DIST), length.out=N_SAMPLING_BIN+1)
closeWeightsNonOVL <- weightsByBin(log10(abs(closePairsNonOVL$dist)), log10(abs(allCloseGenePairsNonOVL$dist)), breaks=breaks)
#~ closeWeightsNonOVL <- getSampleWeightsByDistLogNew(allCloseGenePairsNonOVL, sourcePairs=closePairsNonOVL, adjust=0.1)

# sample close pairs according to distance weight
sampClosePairsNonOVL <- bplapply(1:N_RAND, function(x){ 
    sampleFromAllPairsByWeight(n=nrow(closePairsNonOVL), hitDF=allCloseGenePairsNonOVL, tssGR, weight=closeWeightsNonOVL)
    })
Sys.sleep(1) # hack to fix problems with bplapply on MOGON


#-----------------------------------------------------------------------
# Non-overlapping close pairs
#-----------------------------------------------------------------------
# get all possible gene pairs within MAX_DIST bp
nonParalog <- containsGenePairs(allCloseGenePairsNonOVL, paralogPairs, gPidx=TRUE, nPidx=FALSE, tssGR)
allCloseGenePairsNonPar <- allCloseGenePairsNonOVL[!nonParalog,]

# get sample weights according to distance
#~ closeWeights <- getSampleWeightsByDist(allCloseGenePairs, sourcePairs=closePairs, adjust=DENSITY_BW_ADJUST)
#~ closeWeights <- getSampleWeightsByDistLog(allCloseGenePairs, sourcePairs=closePairs, adjust=1)
#~ closeWeights <- weightsByBin(log10(abs(closePairs$dist)), log10(abs(allCloseGenePairs$dist)), breaks=50)
N_SAMPLING_BIN <- 20
breaks = seq(log10(1), log10(MAX_DIST), length.out=N_SAMPLING_BIN+1)
closeWeightsNonPar <- weightsByBin(log10(abs(closePairsNonOVL$dist)), log10(abs(allCloseGenePairsNonPar$dist)), breaks=breaks)

# sample close pairs according to distance weight
sampClosePairsNonPar <- bplapply(1:N_RAND, function(x){ 
    sampleFromAllPairsByWeight(n=nrow(closePairsNonOVL), hitDF=allCloseGenePairsNonPar, tssGR, weight=closeWeightsNonPar)
    })
Sys.sleep(1) # hack to fix problems with bplapply on MOGON



#-----------------------------------------------------------------------
# Sample cis pairs according to enhancer number in paralogs and linear distance in paralog gene pairs
#-----------------------------------------------------------------------
 
# get sample weights according to enahncer number and distance
#~ closeAndEhWeights = getSampleWeightsByDistAndEnhancers(allCloseGenePairs, tssGR, closePairs, adjust=DENSITY_BW_ADJUST)

closeAndEhWeights = weightsByBinDistAndEnhancers(distWeight=closeWeightsNonPar, sourcePairs=closePairsNonOVL, hitDF=allCloseGenePairsNonPar, tssGR)

# sample according to enahncer number and distance
sampEhClosePairs = bplapply(1:N_RAND, function(x){ 
    sampleFromAllPairsByWeight(n=nrow(closePairs), hitDF=allCloseGenePairs, tssGR, weight=closeAndEhWeights)
    })
Sys.sleep(1) # hack to fix problems with bplapply on MOGON


#-----------------------------------------------------------------------
# Sample distal cis pairs by distance only
#-----------------------------------------------------------------------
# Now sample from all possible gene pairs within DISTAL_MIN_DIST - DISTAL_MAX_DIST bp
allDistalGenePairs = getAllGenePairs(tssGR, maxDist=DISTAL_MAX_DIST, minDist=1)

# get observed probabilty density of distances
#~ distalWeight = getSampleWeightsByDist(hitDF=allDistalGenePairs, sourcePairs=distalPairs, adjust=DENSITY_BW_ADJUST)
#~ distalWeight = getSampleWeightsByDistLog(hitDF=allDistalGenePairs, sourcePairs=distalPairs, adjust=DENSITY_BW_ADJUST)
#~ distalWeight = weightsByBin(log10(abs(distalPairs$dist)), log10(abs(allDistalGenePairs$dist)), breaks=50)
N_SAMPLING_BIN <- 20
breaks = seq(log10(DISTAL_MIN_DIST), log10(DISTAL_MAX_DIST), length.out=N_SAMPLING_BIN+1)
distalWeight = weightsByBin(log10(abs(distalPairs$dist)), log10(abs(allDistalGenePairs$dist)), breaks=breaks)
#~ distalWeight = getSampleWeightsByDistLogNew(hitDF=allDistalGenePairs, sourcePairs=distalPairs, adjust=1)

# sample according to distance
sampDistalPairs = bplapply(1:N_RAND, function(x){ 
    sampleFromAllPairsByWeight(n=nrow(distalPairs), hitDF=allDistalGenePairs, tssGR, weight=distalWeight)
    })
Sys.sleep(1) # hack to fix problems with bplapply on MOGON

#-----------------------------------------------------------------------
# Sample nonPar distal
distalPar <- containsGenePairs(allDistalGenePairs, paralogPairs, gPidx=TRUE, nPidx=FALSE, tssGR)
allDistalGenePairsNonPar <- allDistalGenePairs[!distalPar,]

# get observed probabilty density of distances
#~ distalWeight = getSampleWeightsByDist(hitDF=allDistalGenePairs, sourcePairs=distalPairs, adjust=DENSITY_BW_ADJUST)
#~ distalWeight = getSampleWeightsByDistLog(hitDF=allDistalGenePairs, sourcePairs=distalPairs, adjust=DENSITY_BW_ADJUST)
#~ distalWeight = weightsByBin(log10(abs(distalPairs$dist)), log10(abs(allDistalGenePairs$dist)), breaks=50)
N_SAMPLING_BIN <- 20
breaks = seq(log10(DISTAL_MIN_DIST), log10(DISTAL_MAX_DIST), length.out=N_SAMPLING_BIN+1)
distalWeightNonPar = weightsByBin(log10(abs(distalPairs$dist)), log10(abs(allDistalGenePairsNonPar$dist)), breaks=breaks)

# sample according to distance
sampDistalPairsNonPar = bplapply(1:N_RAND, function(x){ 
    sampleFromAllPairsByWeight(n=nrow(distalPairs), hitDF=allDistalGenePairsNonPar, tssGR, weight=distalWeightNonPar)
    })
Sys.sleep(1) # hack to fix problems with bplapply on MOGON

#-----------------------------------------------------------------------
# annotate sampled gene pairs
#-----------------------------------------------------------------------
#~ 
#~ # annotate common enhancers
#~ sampEhClosePairs = bplapply(sampEhClosePairs, addCommonEnhancer, gene2ehID)
#~ Sys.sleep(1) # hack to fix problems with bplapply on MOGON
#~ sampEhClosePairs = bplapply(sampEhClosePairs, addRelativeEnhancerPosition, tssGR, gene2ehID, ehGR)
#~ Sys.sleep(1) # hack to fix problems with bplapply on MOGON



# combine all sampling replicates to one data frame
sampClosePairsCombined <- do.call("rbind", sampClosePairs)
sampCloseOrgPairsCombined <- do.call("rbind", sampCloseOrgPairs)
sampCloseNonOVLPairsCombined <- do.call("rbind", sampClosePairsNonOVL)
sampCloseNonParPairsCombined <- do.call("rbind", sampClosePairsNonPar)
#~ sampEhClosePairsCombined <- do.call("rbind", sampEhClosePairs)
sampDistalPairsCombined <- do.call("rbind", sampDistalPairs)
sampDistalNonParPairsCombined <- do.call("rbind", sampDistalPairsNonPar)

#=======================================================================
# 4.) Run analysis
#=======================================================================

closeOrgDistDF <- data.frame(
        group = rep(c("paralogs", "sampled"), c(nrow(closePairsOrg), nrow(sampCloseOrgPairsCombined))),
        dist=abs(c(closePairsOrg[,"dist"], sampCloseOrgPairsCombined[,"dist"]))/10^3
)

closeDistDF <- data.frame(
        group = rep(c("paralogs", "sampled"), c(nrow(closePairs), nrow(sampClosePairsCombined))),
        dist=abs(c(closePairs[,"dist"], sampClosePairsCombined[,"dist"]))/10^3
)
closeNonOVLDistDF <- data.frame(
        group = rep(c("paralogs", "sampled"), c(nrow(closePairsNonOVL), nrow(sampCloseNonOVLPairsCombined))),
        dist=abs(c(closePairsNonOVL[,"dist"], sampCloseNonOVLPairsCombined[,"dist"]))/10^3
)
closeNonParDistDF <- data.frame(
        group = rep(c("paralogs", "sampled"), c(nrow(closePairsNonOVL), nrow(sampCloseNonParPairsCombined))),
        dist=abs(c(closePairsNonOVL[,"dist"], sampCloseNonParPairsCombined[,"dist"]))/10^3
)
distalDistDF <- data.frame(
        group = rep(c("paralogs", "sampled"), c(nrow(distalPairs), nrow(sampDistalPairsCombined))),
        dist=abs(c(distalPairs[,"dist"], sampDistalPairsCombined[,"dist"]))/10^3
)
distaNonParlDistDF <- data.frame(
        group = rep(c("paralogs", "sampled"), c(nrow(distalPairs), nrow(sampDistalNonParPairsCombined))),
        dist=abs(c(distalPairs[,"dist"], sampDistalNonParPairsCombined[,"dist"]))/10^3
)


# sampled close by dist org
p1 <- ggplot(closeOrgDistDF, aes(dist, fill=group)) + 
  geom_density(alpha=.5) + theme_bw() + 
  scale_fill_manual(values = COL)
  
p2 <- ggplot(closeOrgDistDF, aes(dist, fill=group)) + 
  geom_density(alpha=.5) + scale_x_log10() + theme_bw() + 
  scale_fill_manual(values = COL)

# sampled close by dist
p3 <- ggplot(closeDistDF, aes(dist, fill=group)) + 
  geom_density(alpha=.5) + theme_bw() + 
  scale_fill_manual(values = COL)
  
p4 <- ggplot(closeDistDF, aes(dist, fill=group)) + 
  geom_density(alpha=.5) + scale_x_log10() + theme_bw() + 
  scale_fill_manual(values = COL)

# nonOVL sampled close by dist
p5 <- ggplot(closeNonOVLDistDF, aes(dist, fill=group)) + 
  geom_density(alpha=.5) + theme_bw() + 
  scale_fill_manual(values = COL)
  
p6 <- ggplot(closeNonOVLDistDF, aes(dist, fill=group)) + 
  geom_density(alpha=.5) + scale_x_log10() + theme_bw() + 
  scale_fill_manual(values = COL)

# nonPar sampled close by dist
p7 <- ggplot(closeNonParDistDF, aes(dist, fill=group)) + 
  geom_density(alpha=.5) + theme_bw() + 
  scale_fill_manual(values = COL2)
  
p8 <- ggplot(closeNonParDistDF, aes(dist, fill=group)) + 
  geom_density(alpha=.5) + scale_x_log10() + theme_bw() + 
  scale_fill_manual(values = COL2)

# sampled distal by dist
p9 <- ggplot(distalDistDF, aes(dist, fill=group)) + 
  geom_density(alpha=.5)  + theme_bw() + 
  scale_fill_manual(values = COL)
  
p10 <- ggplot(distalDistDF, aes(dist, fill=group)) + 
  geom_density(alpha=.5) + scale_x_log10() + theme_bw() + 
  scale_fill_manual(values = COL)

# sampled nonPar distal by dist
p11 <- ggplot(distaNonParlDistDF, aes(dist, fill=group)) + 
  geom_density(alpha=.5)  + theme_bw() + 
  scale_fill_manual(values = COL2)
  
p12 <- ggplot(distaNonParlDistDF, aes(dist, fill=group)) + 
  geom_density(alpha=.5) + scale_x_log10() + theme_bw() + 
  scale_fill_manual(values = COL2)

pdf(paste0(outPrefix, ".sampling_cis_and_distal.density.pdf"), w=9, h=15)
    grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, ncol=2, nrow=6)
dev.off()

pdf(paste0(outPrefix, ".sampling_cis_and_distal.pdf"), w=9, h=5.25)
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

#-----------------------------------------------------------------------
# Hi-C analysis with data from Rao et al. 2014
#-----------------------------------------------------------------------

# Adds Hi-C contact frequencies to a gene pair data set
sampClosePairs <- bplapply(sampClosePairs, addHiCfreq, tssGR, HiClist, inParallel=FALSE, ignoreSameBin=TRUE)
Sys.sleep(1) # hack to fix problems with bplapply on MOGON
sampClosePairs <- bplapply(sampClosePairs, addHiCfreq, tssGR, HiClistNorm, label="HiCnorm", inParallel=FALSE, ignoreSameBin=TRUE)
Sys.sleep(1) # hack to fix problems with bplapply on MOGON

sampClosePairsNonOVL <- bplapply(sampClosePairsNonOVL, addHiCfreq, tssGR, HiClist, inParallel=FALSE, ignoreSameBin=TRUE)
Sys.sleep(1) # hack to fix problems with bplapply on MOGON
sampClosePairsNonOVL <- bplapply(sampClosePairsNonOVL, addHiCfreq, tssGR, HiClistNorm, label="HiCnorm", inParallel=FALSE, ignoreSameBin=TRUE)
Sys.sleep(1) # hack to fix problems with bplapply on MOGON

sampClosePairsNonPar <- bplapply(sampClosePairsNonPar, addHiCfreq, tssGR, HiClist, inParallel=FALSE, ignoreSameBin=TRUE)
Sys.sleep(1) # hack to fix problems with bplapply on MOGON
sampClosePairsNonPar <- bplapply(sampClosePairsNonPar, addHiCfreq, tssGR, HiClistNorm, label="HiCnorm", inParallel=FALSE, ignoreSameBin=TRUE)
Sys.sleep(1) # hack to fix problems with bplapply on MOGON

sampDistalPairs <- bplapply(sampDistalPairs, addHiCfreq, tssGR, HiClist, inParallel=FALSE, ignoreSameBin=TRUE)
Sys.sleep(1) # hack to fix problems with bplapply on MOGON
sampDistalPairs <- bplapply(sampDistalPairs, addHiCfreq, tssGR, HiClistNorm, label="HiCnorm", inParallel=FALSE, ignoreSameBin=TRUE)
Sys.sleep(1) # hack to fix problems with bplapply on MOGON

sampDistalPairsNonPar <- bplapply(sampDistalPairsNonPar, addHiCfreq, tssGR, HiClist, inParallel=FALSE, ignoreSameBin=TRUE)
Sys.sleep(1) # hack to fix problems with bplapply on MOGON
sampDistalPairsNonPar <- bplapply(sampDistalPairsNonPar, addHiCfreq, tssGR, HiClistNorm, label="HiCnorm", inParallel=FALSE, ignoreSameBin=TRUE)
Sys.sleep(1) # hack to fix problems with bplapply on MOGON


# combine all sampling replicates to one data frame
sampClosePairsCombined <- do.call("rbind", sampClosePairs)
sampCloseNonOVLPairsCombined <- do.call("rbind", sampClosePairsNonOVL)
sampCloseNonParPairsCombined <- do.call("rbind", sampClosePairsNonPar)
sampDistalPairsCombined <- do.call("rbind", sampDistalPairs)
sampDistalNonParPairsCombined <- do.call("rbind", sampDistalPairsNonPar)

#------------------------------------------------------------------------
# close pairs
plotDFcloseHiC = data.frame(
        group = rep(c("paralogs", "sampled"), c(nrow(closePairs), nrow(sampClosePairsCombined))),
        HiCraw = c(closePairs[,"HiCfreq"], sampClosePairsCombined[,"HiCfreq"]),
        HiCnorm = c(closePairs[,"HiCnorm"], sampClosePairsCombined[,"HiCnorm"]),
        dist=abs(c(closePairs[,"dist"], sampClosePairsCombined[,"dist"]))/10^3
    )

plotDFcloseHiC$HiCrawNoZero = plotDFcloseHiC$HiCraw
plotDFcloseHiC$HiCrawNoZero[plotDFcloseHiC$HiCraw == 0] = NA
plotDFcloseHiC$HiCnormNoZero = plotDFcloseHiC$HiCnorm
plotDFcloseHiC$HiCnormNoZero[plotDFcloseHiC$HiCnorm == 0] = NA

#------------------------------------------------------------------------
# NonOVL
plotDFcloseNonOVLHiC = data.frame(
        group = rep(c("paralogs", "sampled"), c(nrow(closePairsNonOVL), nrow(sampCloseNonOVLPairsCombined))),
        HiCraw = c(closePairsNonOVL[,"HiCfreq"], sampCloseNonOVLPairsCombined[,"HiCfreq"]),
        HiCnorm = c(closePairsNonOVL[,"HiCnorm"], sampCloseNonOVLPairsCombined[,"HiCnorm"]),
        dist=abs(c(closePairsNonOVL[,"dist"], sampCloseNonOVLPairsCombined[,"dist"]))/10^3
    )

plotDFcloseNonOVLHiC$HiCrawNoZero = plotDFcloseNonOVLHiC$HiCraw
plotDFcloseNonOVLHiC$HiCrawNoZero[plotDFcloseNonOVLHiC$HiCraw == 0] = NA
plotDFcloseNonOVLHiC$HiCnormNoZero = plotDFcloseNonOVLHiC$HiCnorm
plotDFcloseNonOVLHiC$HiCnormNoZero[plotDFcloseNonOVLHiC$HiCnorm == 0] = NA

#------------------------------------------------------------------------
# NonPar
plotDFcloseNonParHiC = data.frame(
        group = rep(c("paralogs", "sampled"), c(nrow(closePairsNonOVL), nrow(sampCloseNonParPairsCombined))),
        HiCraw = c(closePairsNonOVL[,"HiCfreq"], sampCloseNonParPairsCombined[,"HiCfreq"]),
        HiCnorm = c(closePairsNonOVL[,"HiCnorm"], sampCloseNonParPairsCombined[,"HiCnorm"]),
        dist=abs(c(closePairsNonOVL[,"dist"], sampCloseNonParPairsCombined[,"dist"]))/10^3
    )

plotDFcloseNonParHiC$HiCrawNoZero = plotDFcloseNonParHiC$HiCraw
plotDFcloseNonParHiC$HiCrawNoZero[plotDFcloseNonParHiC$HiCraw == 0] = NA
plotDFcloseNonParHiC$HiCnormNoZero = plotDFcloseNonParHiC$HiCnorm
plotDFcloseNonParHiC$HiCnormNoZero[plotDFcloseNonParHiC$HiCnorm == 0] = NA
    
#------------------------------------------------------------------------
# distal pairs

plotDFdistalHiC = data.frame(
        group = c(rep("paralogs", nrow(distalPairs)), rep("sampled", nrow(sampDistalPairsCombined))),
        HiCraw = c(distalPairs[,"HiCfreq"], sampDistalPairsCombined[,"HiCfreq"]),
        HiCnorm = c(distalPairs[,"HiCnorm"], sampDistalPairsCombined[,"HiCnorm"]),
        dist=abs(c(distalPairs[,"dist"], sampDistalPairsCombined[,"dist"]))/10^3

    )
    
plotDFdistalHiC$HiCrawNoZero = plotDFdistalHiC$HiCraw
plotDFdistalHiC$HiCrawNoZero[plotDFdistalHiC$HiCraw == 0] = NA
plotDFdistalHiC$HiCnormNoZero = plotDFdistalHiC$HiCnorm
plotDFdistalHiC$HiCnormNoZero[plotDFdistalHiC$HiCnorm == 0] = NA

#------------------------------------------------------------------------
# non Par distal pairs

plotDFdistalNonParHiC = data.frame(
        group = c(rep("paralogs", nrow(distalPairs)), rep("sampled", nrow(sampDistalNonParPairsCombined))),
        HiCraw = c(distalPairs[,"HiCfreq"], sampDistalNonParPairsCombined[,"HiCfreq"]),
        HiCnorm = c(distalPairs[,"HiCnorm"], sampDistalNonParPairsCombined[,"HiCnorm"]),
        dist=abs(c(distalPairs[,"dist"], sampDistalNonParPairsCombined[,"dist"]))/10^3

    )
plotDFdistalNonParHiC$HiCrawNoZero = plotDFdistalNonParHiC$HiCraw
plotDFdistalNonParHiC$HiCrawNoZero[plotDFdistalNonParHiC$HiCraw == 0] = NA
plotDFdistalNonParHiC$HiCnormNoZero = plotDFdistalNonParHiC$HiCnorm
plotDFdistalNonParHiC$HiCnormNoZero[plotDFdistalNonParHiC$HiCnorm == 0] = NA

#------------------------------------------------------------------------
# DEBUG: plot dist difference


# sampled close by dist
close1 <- ggplot(closeDistDF, aes(dist, fill=group)) + 
  geom_density(alpha=.5) + theme_bw() + 
  scale_fill_manual(values = COL)
  
close2 <- ggplot(closeDistDF, aes(dist, fill=group)) + 
  geom_density(alpha=.5) + scale_x_log10() + theme_bw() + 
  scale_fill_manual(values = COL)
close3 <- dotplotWithDensityLogXY(plotDFcloseHiC[nrow(plotDFcloseHiC):1,], "dist", "HiCraw", "group", COL)
close4 <- dotplotWithDensityLogXY(plotDFcloseHiC[nrow(plotDFcloseHiC):1,], "dist", "HiCnorm", "group", COL)


# nonOVL sampled close by dist
nonOVL1 <- ggplot(closeNonOVLDistDF, aes(dist, fill=group)) + 
  geom_density(alpha=.5) + theme_bw() + 
  scale_fill_manual(values = COL)
  
nonOVL2 <- ggplot(closeNonOVLDistDF, aes(dist, fill=group)) + 
  geom_density(alpha=.5) + scale_x_log10() + theme_bw() + 
  scale_fill_manual(values = COL)

nonOVL3 <- dotplotWithDensityLogXY(plotDFcloseNonOVLHiC[nrow(plotDFcloseNonOVLHiC):1,], "dist", "HiCraw", "group", COL)
nonOVL4 <- dotplotWithDensityLogXY(plotDFcloseNonOVLHiC[nrow(plotDFcloseNonOVLHiC):1,], "dist", "HiCnorm", "group", COL)


# nonPar sampled close by dist
closeNonPar1 <- ggplot(closeNonParDistDF, aes(dist, fill=group)) + 
  geom_density(alpha=.5) + theme_bw() + 
  scale_fill_manual(values = COL2)
  
closeNonPar2 <- ggplot(closeNonParDistDF, aes(dist, fill=group)) + 
  geom_density(alpha=.5) + scale_x_log10() + theme_bw() + 
  scale_fill_manual(values = COL2)

closeNonPar3 <- dotplotWithDensityLogXY(plotDFcloseNonParHiC[nrow(plotDFcloseNonParHiC):1,], "dist", "HiCraw", "group", COL2)
closeNonPar4 <- dotplotWithDensityLogXY(plotDFcloseNonParHiC[nrow(plotDFcloseNonParHiC):1,], "dist", "HiCnorm", "group", COL2)


# sampled distal by dist
dist1 <- ggplot(distalDistDF, aes(dist, fill=group)) + 
  geom_density(alpha=.5)  + theme_bw() + 
  scale_fill_manual(values = COL)
  
dist2 <- ggplot(distalDistDF, aes(dist, fill=group)) + 
  geom_density(alpha=.5) + scale_x_log10() + theme_bw() + 
  scale_fill_manual(values = COL)
dist3 <- dotplotWithDensityLogXY(plotDFdistalHiC[nrow(plotDFdistalHiC):1,], "dist", "HiCraw", "group", COL)
dist4 <- dotplotWithDensityLogXY(plotDFdistalHiC[nrow(plotDFdistalHiC):1,], "dist", "HiCnorm", "group", COL)


# sampled nonPar distal by dist
distNonPar1 <- ggplot(distaNonParlDistDF, aes(dist, fill=group)) + 
  geom_density(alpha=.5)  + theme_bw() + 
  scale_fill_manual(values = COL2)
  
distNonPar2 <- ggplot(distaNonParlDistDF, aes(dist, fill=group)) + 
  geom_density(alpha=.5) + scale_x_log10() + theme_bw() + 
  scale_fill_manual(values = COL2)

distNonPar3 <- dotplotWithDensityLogXY(plotDFdistalNonParHiC[nrow(plotDFdistalNonParHiC):1,], "dist", "HiCraw", "group", COL2)
distNonPar4 <- dotplotWithDensityLogXY(plotDFdistalNonParHiC[nrow(plotDFdistalNonParHiC):1,], "dist", "HiCnorm", "group", COL2)

# all
# sampled nonPar distal by dist
all1 <- ggplot(rbind(plotDFcloseHiC, plotDFdistalHiC)[(nrow(plotDFcloseHiC)+nrow(plotDFdistalHiC)):1,], aes(dist, fill=group)) + 
  geom_density(alpha=.5)  + theme_bw() + 
  scale_fill_manual(values = COL)
  
all2 <- ggplot(rbind(plotDFcloseHiC, plotDFdistalHiC)[(nrow(plotDFcloseHiC)+nrow(plotDFdistalHiC)):1,], aes(dist, fill=group)) + 
  geom_density(alpha=.5) + scale_x_log10() + theme_bw() + 
  scale_fill_manual(values = COL)

all3 <- dotplotWithDensityLogXY(rbind(plotDFcloseHiC, plotDFdistalHiC)[(nrow(plotDFcloseHiC)+nrow(plotDFdistalHiC)):1,], "dist", "HiCraw", "group", COL)

all4 <- dotplotWithDensityLogXY(rbind(plotDFcloseHiC, plotDFdistalHiC)[(nrow(plotDFcloseHiC)+nrow(plotDFdistalHiC)):1,], "dist", "HiCnorm", "group", COL)


pdf(paste0(outPrefix, ".close_and_distal_pairs.Hi-C_vs_dist_all.ggboxplot.pdf"), w=30, h=40)
#~     grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12, ncol=2, nrow=6) 
    grid.arrange(
        close1,close2,close3,close4, 
        nonOVL1,nonOVL2,nonOVL3,nonOVL4, 
        closeNonPar1,closeNonPar2,closeNonPar3,closeNonPar4, 
        dist1,dist2,dist3,dist4, 
        distNonPar1,distNonPar2,distNonPar3,distNonPar4, 
        all1,all2,all3,all4, 
        ncol=4, nrow=6) 
dev.off()

wilcox.test(HiCnormNoZero ~ group, data=plotDFcloseHiC)
wilcox.test(HiCnormNoZero ~ group, data=plotDFcloseNonOVLHiC)
wilcox.test(HiCnormNoZero ~ group, data=plotDFcloseNonParHiC)
wilcox.test(HiCnormNoZero ~ group, data=plotDFdistalHiC)
wilcox.test(HiCnormNoZero ~ group, data=plotDFdistalNonParHiC)
wilcox.test(HiCnormNoZero ~ group, data=rbind(plotDFcloseHiC, plotDFdistalHiC))

#=======================================================================
# save workspace image
#=======================================================================
save.image(WORKIMAGE_FILE)
