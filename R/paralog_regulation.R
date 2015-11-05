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

# ISSUES TO BE FIXED BEFORE FINAL NUMBERS:
# - TAD bed file and ranges +/- 1 error (Rao org data use 0-based including coords)
# - make all gene pair data frames as.character()
# - check mapping of enhancers to genes with tssGR object and id remapping

# TO FINALIZE PIPELINE FOR SUBMISSION:
# X redesign code to run on MOGON server
# X use orghologMouse instead of orthologAll data set
# X use seed() command for reproducible randomizations
# - remove the following parts from the analysis (if not included in manuscript):
#   X TF motif analysis
#   X conserved TADs over all cell types
#   x MIC score of expression correlation
# - put sameTAD annotations in the cisPairs not in the GR and combine replicated sampled pairs early

require(biomaRt)        # to retrieve human paralogs from Ensembl
require(stringr)        # for some string functionality
require(RColorBrewer)   # for nice colors
require(GenomicRanges)  # for genomic intervals and overlap calculation
require(rtracklayer)    # for import.bed
require(plyr)           # count() function
require(gplots)         # heatmap.2 function
require(data.table)     # for data.table object
require(VennDiagram)    # for area-proportional Venn-Diagrams 
require(ggplot2)        # for nice plots
require(BiocParallel)   # for parallel computing

#-----------------------------------------------------------------------
# Load parameters from external script
#-----------------------------------------------------------------------
# read parameter script location from command-line argument
#~ args <- commandArgs(trailingOnly = TRUE)
#~ PARAM_SCRIPT=args[1]

PARAM_SCRIPT="R/paralog_regulation.param.v11.R"
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
ehGR = import.bed(EH_FILE_FANTOM5, seqinfo=seqInfoRealChrom)

# parse regulatory map from FANTOM5
fantomData = parseFANTOM5RegMap(REGMAP_FILE_FANTOM5, allGenes, ehGR, seqInfoRealChrom)

map=fantomData[["map"]]

# remap indices of used genes to tssGR indices
usedGeneIDs = names(fantomData[["tssGR"]])[map[,1]]
map[,1] = match(usedGeneIDs, names(tssGR))

# remove genes that could not be mapped back to the ENSG ID
map = map[!is.na(map[,1]),]
map$dist = getMapDist(map, tssGR, ehGR)

# make GR for associations
mapGR = getMapAsGR(map, tssGR, ehGR, strand.as.direction=TRUE)

# map gene symbols to linked enhancer IDs 
gene2ehID = getGenetoEhIDmapping(names(tssGR)[map[,1]], map[,2])

# add number of linked enhancer to tssGR
tssGR$linked_enhancer =  sapply(gene2ehID[names(tssGR)], length)

#-----------------------------------------------------------------------
# Parse Hi-C domains
#-----------------------------------------------------------------------

# parse TAD data sets as list of GRanges
RaoTADs = bplapply(RaoDomainFiles, parseDomainsRao, disjoin=FALSE, seqinfo=seqInfo)
Sys.sleep(3)
DixonTADs <- bplapply(DixonDomainFiles, import.bed, seqinfo=seqInfo)
Sys.sleep(3)

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

#~ allBoundaries = lapply(allTADs, getNonOverlappingBoundaries, min.gapwidth=1)

# write all TADs to BED files:
for (tadName in names(allTADs)){
    TAD = allTADs[[tadName]]
    export(sort(TAD), paste0(outPrefix, ".TAD_data.", tadName, ".bed"))
}

orgStr2Name = c(mmusculus="mouse", cfamiliaris="dog")

# parse mouse and dog TADs
speciesSeqInfo = list("mmusculus"=seqInfoMouse, "cfamiliaris"=seqInfoDog)


speciesTADs = lapply(1:2, function(i)
    parseRudanTADs(RudanFile, sheet=i, disjoin=FALSE, seqinfo=speciesSeqInfo[[i]])
    )
names(speciesTADs) = names(speciesSeqInfo)

speciesHiC = list(
    "mmusculus" = parseRudanHiC("data/Rudan2015/GSE65126_HiC_mouse_liver_merged_50000.txt", seqInfoMouse, resolution=50*10^3), 
    "cfamiliaris"= parseRudanHiC("data/Rudan2015/GSE65126_HiC_dog_liver_merged_50000.txt", seqInfoDog, resolution=50*10^3)
)

#-----------------------------------------------------------------------
# parse baseline gene expression data set:
#-----------------------------------------------------------------------

# expDFlist is already loaded in the script "R/data.expression.R"
nExp = length(expDFlist)

# plot overlap of ENSG IDs in expression data and gene set:
for (expName in names(expDFlist)) {

    expDF = expDFlist[[expName]]

    setList = list(names(tssGR), row.names(expDF))
    names(setList) = c("ENSEMBL", paste(expName, "expression"))
    vennPlotFile = paste0(outPrefix, ".", expName, ".geneIDs_ensembl_vs_expression.vennDiagram.png")
    # plot venn diagram
    v = venn.diagram(setList, filename=vennPlotFile, imagetype="png", fill=COL2, 
            main=paste("Ensembl Gene ID overlap in", expName), main.cex=1.3, main.fontfamily=3,
            cex=1.3, cat.cex=1.3, fontfamily=3, cat.fontfamily=1.3, cat.default.pos="text", 
        )
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
#             +-- cisPairs           
#             +-- distalCisPairs           


# paralogPairs are loaded from data.ensembl.R

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

# add HGNC symbols
allCisPairs = addHGNC(allCisPairs, tssGR)

# add same Strand info:
allCisPairs = addSameStrand(allCisPairs, tssGR)

# add number of common enhancers
allCisPairs = addCommonEnhancer(allCisPairs, gene2ehID)

# add position of enhancers:
allCisPairs = addRelativeEnhancerPosition(allCisPairs, tssGR, gene2ehID, ehGR)

# add pairwise correlations of gene expression over all tissues
for (expName in names(expDFlist)) {
    message(paste("INFO: annotate pairs with expression correlation form:", expName))
    expDF = expDFlist[[expName]]
    
    allCisPairs = addCor(allCisPairs, expDF, colName=paste0(expName, "_expCor"))
     
}

# add information of one-two-one orthologs in other species
allCisPairs = addOrthologAnnotation(allCisPairs, orthologsSpeciesList[["mmusculus"]], "mmusculus", tssGRmouse, speciesTADs[["mmusculus"]], speciesHiC[["mmusculus"]][[1]], speciesHiC[["mmusculus"]][[2]] )
allCisPairs = addOrthologAnnotation(allCisPairs, orthologsSpeciesList[["cfamiliaris"]], "cfamiliaris", tssGRdog, speciesTADs[["cfamiliaris"]], speciesHiC[["cfamiliaris"]][[1]], speciesHiC[["cfamiliaris"]][[2]])

# add promoter-promoter contacts from caputre Hi-C
allCisPairs[,"captureC_raw"] = getPairwiseMatrixScore(allCisPairs, captureHiC[["raw"]], tssGR, replaceZeroByNA=TRUE)
allCisPairs[,"captureC_ObsExp"] = getPairwiseMatrixScore(allCisPairs, captureHiC[["obsExp"]], tssGR, replaceZeroByNA=TRUE)

#-----------------------------------------------------------------------
# Filter for close and distal pairs
#-----------------------------------------------------------------------

# get close cis pairs
cisPairs = allCisPairs[abs(allCisPairs$dist) <= MAX_DIST,] 

# get distal pairs
distalCisPairs = allCisPairs[abs(allCisPairs$dist) > MAX_DIST,] 

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

# save before sampling
save.image(paste0(WORKIMAGE_FILE, ".only_annotation.Rdata"))

#=======================================================================
# 3.) Sample random control/background data sets
#=======================================================================

#-----------------------------------------------------------------------
# Sample pairs with equal probability from all genes
#-----------------------------------------------------------------------
randPairs = bplapply(1:N_RAND, function(x){getRandomPairs(nrow(paralogPairsUniqGuniqP), names(tssGR))})
Sys.sleep(3)

randPairs = lapply(randPairs, addSameStrand, tssGR)

# filter for random pairs in Cis
randPairsInCis = lapply(randPairs, getCisPairs, tssGR)

randPairsInCis <- bplapply(randPairsInCis, addOrthologAnnotation, orthologsSpeciesList[["mmusculus"]], "mmusculus",  tssGRmouse, speciesTADs[["mmusculus"]], speciesHiC[["mmusculus"]][[1]], speciesHiC[["mmusculus"]][[2]], inParallel=FALSE)
Sys.sleep(3)
randPairsInCis <- bplapply(randPairsInCis, addOrthologAnnotation, orthologsSpeciesList[["cfamiliaris"]], "cfamiliaris", tssGRdog, speciesTADs[["cfamiliaris"]], speciesHiC[["cfamiliaris"]][[1]], speciesHiC[["cfamiliaris"]][[2]], inParallel=FALSE)
Sys.sleep(3)

#-----------------------------------------------------------------------
# Sample cis pairs according to enhancer number in paralogs and linear distance in paralog gene pairs
#-----------------------------------------------------------------------

# get all possible gene pairs within MAX_DIST bp
allGenePairs = getAllGenePairs(tssGR, maxDist=MAX_DIST)

# get sample weights according to enahncer number and distance
cisWeights = getSampleWeightsByDistAndEnhancers(allGenePairs, tssGR, cisPairs, adjust=DENSITY_BW_ADJUST)

# sample according to enahncer number and distance
randCisPairs = bplapply(1:N_RAND, function(x){ 
    sampleFromAllPairsByWeight(n=nrow(cisPairs), hitDF=allGenePairs, tssGR, weight=cisWeights)
    })
Sys.sleep(3)

#-----------------------------------------------------------------------
# Sample distal cis pairs the same way
#-----------------------------------------------------------------------

# Now sample from all possible gene pairs within DISTAL_MIN_DIST - DISTAL_MAX_DIST bp
allDistalGenePairs = getAllGenePairs(tssGR, maxDist=DISTAL_MAX_DIST, minDist=DISTAL_MIN_DIST)


# get observed probabilty density of distances and linked enhancer
#~ distalWeight = getSampleWeightsByDistAndEnhancers(allDistalGenePairs, tssGR, distalCisPairs, adjust=DENSITY_BW_ADJUST)

distalWeight = getSampleWeightsByDist(hitDF=allDistalGenePairs, sourcePairs=distalCisPairs, adjust=DENSITY_BW_ADJUST)

# sample according to enahncer number and distance
#~ randDistalCisPairs = replicate(N_RAND, 
#~     sampleFromAllPairsByWeight(n=nrow(distalCisPairs), hitDF=allDistalGenePairs, tssGR, weight=distalWeight)
#~     , simplify=FALSE)

randDistalCisPairs = bplapply(1:N_RAND, function(x){ 
    sampleFromAllPairsByWeight(n=nrow(distalCisPairs), hitDF=allDistalGenePairs, tssGR, weight=distalWeight)
    })
Sys.sleep(3)


#-----------------------------------------------------------------------
# annotate sampled gene pairs
#-----------------------------------------------------------------------
randCisPairs = lapply(randCisPairs, addHGNC, tssGR)
randCisPairs = lapply(randCisPairs, addSameStrand, tssGR)
randCisPairs = bplapply(randCisPairs, addCommonEnhancer, gene2ehID)
Sys.sleep(3)
randCisPairs = bplapply(randCisPairs, addRelativeEnhancerPosition, tssGR, gene2ehID, ehGR)
Sys.sleep(3)

for (expName in names(expDFlist)) {
    expDF = expDFlist[[expName]]
    
    message(paste("INFO: Annotate sampled pairs with expression form:", expName))
    randCisPairs = bplapply(randCisPairs, addCor, expDF, colName=paste0(expName, "_expCor"))
    Sys.sleep(3)

}

# save temp
save.image(paste0(WORKIMAGE_FILE, ".sampling_and_after_expression_annotation.Rdata"))

randCisPairs = bplapply(randCisPairs, addOrthologAnnotation, orthologsSpeciesList[["mmusculus"]], "mmusculus", tssGRmouse, speciesTADs[["mmusculus"]], speciesHiC[["mmusculus"]][[1]], speciesHiC[["mmusculus"]][[2]], inParallel=FALSE)
Sys.sleep(3)

randCisPairs = bplapply(randCisPairs, addOrthologAnnotation, orthologsSpeciesList[["cfamiliaris"]], "cfamiliaris", tssGRdog, speciesTADs[["cfamiliaris"]], speciesHiC[["cfamiliaris"]][[1]], speciesHiC[["cfamiliaris"]][[2]], inParallel=FALSE)
Sys.sleep(3)


# add promoter-promoter contacts from caputre Hi-C
randCisPairs <- bplapply(1:N_RAND, function(i){
    randCisPairs[[i]][,"captureC_raw"] <- getPairwiseMatrixScore(randCisPairs[[i]], captureHiC[["raw"]], tssGR, replaceZeroByNA=TRUE)
    randCisPairs[[i]][,"captureC_ObsExp"] <- getPairwiseMatrixScore(randCisPairs[[i]], captureHiC[["obsExp"]], tssGR, replaceZeroByNA=TRUE)
    return(randCisPairs[[i]])
})
Sys.sleep(3)

# same for sampled distal pairs
randDistalCisPairs <- bplapply(1:N_RAND, function(i){
    randDistalCisPairs[[i]][,"captureC_raw"] <- getPairwiseMatrixScore(randDistalCisPairs[[i]], captureHiC[["raw"]], tssGR, replaceZeroByNA=TRUE)
    randDistalCisPairs[[i]][,"captureC_ObsExp"] <- getPairwiseMatrixScore(randDistalCisPairs[[i]], captureHiC[["obsExp"]], tssGR, replaceZeroByNA=TRUE)
    return(randDistalCisPairs[[i]])
})
Sys.sleep(3)


randDistalCisPairs = lapply(randDistalCisPairs, addHGNC, tssGR)
randDistalCisPairs = lapply(randDistalCisPairs, addSameStrand, tssGR)

randDistalCisPairs = bplapply(randDistalCisPairs, addOrthologAnnotation, orthologsSpeciesList[["mmusculus"]], "mmusculus", tssGRmouse, speciesTADs[["mmusculus"]], speciesHiC[["mmusculus"]][[1]], speciesHiC[["mmusculus"]][[2]], inParallel=FALSE)
Sys.sleep(3)

randDistalCisPairs = bplapply(randDistalCisPairs, addOrthologAnnotation, orthologsSpeciesList[["cfamiliaris"]], "cfamiliaris", tssGRdog, speciesTADs[["cfamiliaris"]], speciesHiC[["cfamiliaris"]][[1]], speciesHiC[["cfamiliaris"]][[2]], inParallel=FALSE)
Sys.sleep(3)


#-----------------------------------------------------------------------
# save a work image after sampling and annotation before Hi-C.
#-----------------------------------------------------------------------
save.image(paste0(WORKIMAGE_FILE, ".sampling_and_annotation_beforeHiC.Rdata"))

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


# combine all sampling replicates to one data frame
randDistalCisPairsCombined = do.call("rbind", randDistalCisPairs)
randCisPairsCombined = do.call("rbind", randCisPairs)

if ( !USE_LOCAL_HIC_CONTACTS) {

    # Adds Hi-C contact frequencies to a gene pari data set
    distalCisPairs = addHiCfreq(distalCisPairs, tssGR, HiClist)
    randDistalCisPairsCombined = addHiCfreq(randDistalCisPairsCombined, tssGR, HiClist)
    #~ randCisPairs = lapply(randCisPairs, addHiCfreq, tssGR, HiClist)
    
    distalCisPairs = addHiCfreq(distalCisPairs, tssGR, HiClistNorm, label="HiCnorm")
    randDistalCisPairsCombined = addHiCfreq(randDistalCisPairsCombined, tssGR, HiClistNorm, label="HiCnorm")

    # save or load downloaded data 
    save(distalCisPairs, randDistalCisPairsCombined, file=paste0(outDataPrefix, ".Hi-C.", CELL, ".", HIC_RESOLUTION, ".distal_and_sampled_pairs.RData"))

}else{
    load(paste0(outDataPrefix, ".Hi-C.", CELL, ".", HIC_RESOLUTION, ".distal_and_sampled_pairs.RData"))
}

#-----------------------------------------------------------------------
# save a work image after sampling and annotation.
#-----------------------------------------------------------------------
save.image(paste0(WORKIMAGE_FILE, ".sampling_and_annotation.Rdata"))

#=======================================================================
# 4.) Run analysis
#=======================================================================

#-----------------------------------------------------------------------
# verify same distribution of linked enhancer in paralog genes and randomly sampled genes
#-----------------------------------------------------------------------
pdf(paste0(outPrefix, ".sampling_cis_dist.pdf"))
    par(cex=1, lwd=1.5, mfrow=c(2,1))
    # verify same distribution of distances
    paraDist = cisPairs$dist / 10^3
    randDist = unlist(lapply(randCisPairs, function(d){d[,"dist"]})) / 10^3
    hist(abs(paraDist), 50, col=COL[1],
    main="Paralog pairs", xlab="Distance (kb)")
    hist(abs(randDist), 50, col=COL[2],
    main="Sampled pairs", xlab="Distance (kb)")    
dev.off()

# get distributions of linked enhancers and distances
paraLinkedEnhancer = tssGR[c(cisPairs[,1], cisPairs[,2])]$linked_enhancer
randLinkedEnhancer = unlist(lapply(randCisPairs, function(gP) tssGR[c(gP[,1], gP[,2])]$linked_enhancer))
paraDist = cisPairs$dist / 10^3
randDist = unlist(lapply(randCisPairs, function(d){d[,"dist"]})) / 10^3
paraDistalDist = distalCisPairs$dist / 10^3
randDistalDist = unlist(lapply(randDistalCisPairs, function(d){d[,"dist"]})) / 10^3

pdf(paste0(outPrefix, ".sampling_cis_and_distal.pdf"))
    par(cex=1, lwd=1.5, mfrow=c(3,3))

    qqplot(paraLinkedEnhancer, randLinkedEnhancer, xlab="Enhancers in paralog genes", ylab="Enhancers in sampled genes", main="QQ-Plot of linked enhancers\n in close cis pairs")
    abline(0,1, col="red")

    hist(paraLinkedEnhancer[paraLinkedEnhancer<=50], 50, col=COL[1],
    main="Close cis paralogs", xlab="Enhancers")
    hist(randLinkedEnhancer[randLinkedEnhancer<=50], 50, col=COL[2],
    main="Close cis sampled genes", xlab="Enhancers")    

    # verify same distribution of distances    
    qqplot(abs(paraDist), abs(randDist), xlab="Disance between paralogs", ylab="Distance between sampled gene pairs", main="QQ-Plot of distances\n in in close cis pairs")
    abline(0,1, col="red")

    hist(abs(paraDist), 50, col=COL[1],
    main="Close cis paralogs", xlab="Distance (kb)")
    hist(abs(randDist), 50, col=COL[2],
    main="Close cis sampled genes", xlab="Distance (kb)")    

    # same for distal pairs
    qqplot(abs(paraDistalDist), abs(randDistalDist), log="xy",
    xlab="Disance between paralogs", ylab="Distance between sampled gene pairs", main="QQ-Plot of distances\n in distal cis pairs")
    abline(0,1, col="red")

    hist(abs(paraDistalDist), 50, col=COL[1],
    main="Distal cis paralogs", xlab="Distance (kb)")
    hist(abs(randDistalDist), 50, col=COL[2],
    main="Distal cis sampled genes", xlab="Distance (kb)")    

dev.off()

pdf(paste0(outPrefix, ".sampling_cis_and_distal.qqplot.pdf"))
    par(cex=1, lwd=1.5, mfrow=c(2,2))
    
    closeLim = c(0, MAX_DIST)/10^3
    distalLim = c(DISTAL_MIN_DIST, DISTAL_MAX_DIST)/10^3
    
    qqplot(abs(paraDist), abs(randDist), xlim=closeLim, ylim=closeLim, xlab="Disance between paralogs", ylab="Distance between sampled gene pairs", main="QQ-Plot of distances\n in in close cis pairs")
    abline(0,1, col="red")

    qqplot(log10(abs(paraDist)+1), log10(abs(randDist)+1), xlim=log10(closeLim+1), ylim=log10(closeLim+1), xlab="log10 Disance between paralogs", ylab=" log10 Distance between sampled gene pairs", main="QQ-Plot of distances\n in in close cis pairs")
    abline(0,1, col="red")

    qqplot(abs(paraDistalDist), abs(randDistalDist),
    xlab="Disance between paralogs", ylab="Distance between sampled gene pairs", main="QQ-Plot of distances\n in distal cis pairs")
    abline(0,1, col="red")

    qqplot(log10(abs(paraDistalDist)+1), log10(abs(randDistalDist)+1), xlim=log10(distalLim+1), ylim=log10(distalLim+1),
    xlab="Disance between paralogs", ylab="Distance between sampled gene pairs", main="QQ-Plot of distances\n in distal cis pairs")
    abline(0,1, col="red")
    
dev.off()

#-----------------------------------------------------------------------
# plot number of genes in each group
#-----------------------------------------------------------------------
pdf(paste0(outPrefix, ".number_of_genes.pdf"), width=3.5)

    nrGenes = sapply(list(paralogs, nonParalogs), length)
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

    paraPercent = 100 * nrow(allCisPairs) / nrow(paralogPairsUniqGuniqP)

    
    randPercent = sapply(randPairsInCis, nrow) * 100 / nrow(randPairs[[1]])

    height=c(paraPercent, mean(randPercent))

    par(cex=1.5, lwd=2)
    bp = my.barplot(height, 
        names=paste0(c("Paralogs\n (n=", "Random pairs\n (n="), c(nrow(paralogPairsUniqGuniqP), paste0(N_RAND, "x",nrow(paralogPairsUniqGuniqP))), c(")", ")")), 
        addValues=TRUE, col=COL_RAND,
        main="Gene pairs on\n same chromosome", 
        ylab="%")
    
    error.bar(bp,height, c(NA, sd(randPercent)))
        
dev.off()

#---------------------------------------------------------------
# fraction of paralog pairs with same strand on same chrom
#---------------------------------------------------------------
paraPercent = sapply(list(paralogPairsUniqGuniqP, distalCisPairs, cisPairs), function(gP) percentTrue(gP[,"sameStrand"]))
randPercent = lapply(list(randPairs, randDistalCisPairs,randCisPairs), function(gpl) {
        sapply(gpl, function(gP) percentTrue(gP[,"sameStrand"]))
        })

height=rbind(paraPercent, sapply(randPercent, mean))
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
wstest = wilcox.test(abs(cisPairs$dist[cisPairs$sameStrand]), abs(cisPairs$dist[!cisPairs$sameStrand]))
pdf(paste0(outPrefix, ".paralogPairs_sameStrand_dist.boxplot.pdf"), width=3.5)
    par(cex=1.5, lwd=2)
    bp = my.boxplot(list(abs(cisPairs$dist[cisPairs$sameStrand])/10^3 ,abs(cisPairs$dist[!cisPairs$sameStrand])/10^3),
        returnBp=TRUE, border=COL[1], names=c("Same strand", "Opposite strand"),
        ylab="Distance [kb]")
    
    add_pval_all_allPairs(1:length(bp$n), 1.1*max(bp$out), pval_mat=cbind(NA, signif(wstest$p.value, 3)), xpd=TRUE)

dev.off()

#-----------------------------------------------------------------------
# Paralog pairs across chromosomes
#-----------------------------------------------------------------------

# get number of pairs across each pair of chromosomes
interChromPairsMatrix = interChromPairMatrix(paralogPairsUniqGuniqP, tssGR)
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
Sys.sleep(3)

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
paraDist = allCisPairs[abs(allCisPairs[,"dist"]) <= MAX_DIST ,"dist"]

# get cispairs and distance form uniformaly random genes
uniformRandCisPairs = lapply(randPairs, getCisPairs, tssGR)
uniformRandCisPairs = lapply(uniformRandCisPairs, addPairDist, tssGR)
randDist = unlist(lapply(uniformRandCisPairs, function(d){d[,"dist"]}))
randDist = randDist[abs(randDist) <= MAX_DIST]

sampledDist = unlist(lapply(randCisPairs, function(d){d[,"dist"]}))
    
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
allRandDists = unlist(lapply(randCisPairs, function(d){d[,"dist"]}))


for (D in DIST_TH) {
    pdf(paste0(outPrefix, ".paralogPairs_cis_within_", D/10^3, "kb.hist.pdf"))
        par(cex=1.3, lwd=2, mfrow=c(2,1))
        hist(abs(cisPairs[abs(cisPairs$dist) <= D, "dist"])/10^3, 50, col=COL[1], main=paste("Linear distance distribution between paralogs\n on same chromosome within", D/10^3, "kb"), xlab="Distance (kb)")

        hist(abs(allRandDists[abs(allRandDists) <= D])/10^3, 50, col=COL[2], main=paste("Linear distance distribution between sampled gene pairs\n on same chromosome within", D/10^3, "kb"), xlab="Distance (kb)")
    dev.off()
}

pdf(paste0(outPrefix, ".paralogPairs_on_same_chrom_distance.boxplot.pdf"))
    par(cex=1.3, lwd=2)
    dists = lapply(c(list(abs(cisPairs$dist)), lapply(randCisPairs, function(df) abs(df[,"dist"]))), function(x){x/10^6}) 
    boxplot(dists, border=c(COL[1], rep(COL[2], N_RAND)),
        names=c("Paralog\n Pairs", 1:N_RAND),
        ylab="Distance (Mb)",main="Linear distance distribution\n between gene pairs on same chromosome"
    )
    mtext("Random pairs", side=1, line=2, at=mean(1:N_RAND +1), cex=1.3)
dev.off()

pdf(paste0(outPrefix, ".paralogPairs_on_same_chrom_distance.log_boxplot.pdf"))
    par(cex=1.3, lwd=2)
    dists = lapply(c(list(log10(abs(cisPairs$dist)+1)), lapply(randCisPairs, function(df) log10(abs(df[,"dist"])+1))), function(x){x/10^6}) 
    boxplot(dists, border=c(COL[1], rep(COL[2], N_RAND)),
        names=c("Paralog\n Pairs", 1:N_RAND),
        ylab="Distance (log_10 Mb)",main="Linear distance distribution\n between gene pairs on same chromosome"
    )
    mtext("Random pairs", side=1, line=2, at=mean(1:N_RAND +1), cex=1.3)
dev.off()

#-----------------------------------------------------------------------
# Plot size distribution of domains
#-----------------------------------------------------------------------
sizeList = lapply(allTADs, function(gr) width(gr)/10^3)
numberList = lapply(allTADs, length)

pdf(paste0(outPrefix, ".Compare_domain_size.boxplot.pdf"))
    par(cex=1.3, lwd=1.5)
    my.boxplot(sizeList, names=paste0(gsub("_", " ", names(sizeList)), "n=", numberList),
        main="Size distribution of TADs",
        ylab="TAD size [kb]", col=COL_DOMAIN)
dev.off()

pdf(paste0(outPrefix, ".Compare_domain_size.log10.boxplot.pdf"))
    par(cex=1.3, lwd=1.5)
    my.boxplot(lapply(sizeList, function(x){log10(x*10^3)}), names=paste0(gsub("_", " ", names(sizeList)), "n=", numberList),
        main="Size distribution of TADs",
        ylab="TAD size (log_10 bp)", col=COL_DOMAIN)
dev.off()

#-----------------------------------------------------------------------
# Co-occurance in same TAD
#-----------------------------------------------------------------------

# make GRanges objects for cis paralog pairs and random paris on same chromosome
cisParaGR = getPairAsGR(cisPairs, tssGR)
cisRandGR = bplapply(randCisPairs, getPairAsGR, tssGR)
Sys.sleep(3)


#~ # co-occurance within the same domain
#~ cisParaGR = addWithinSubject(cisParaGR, hESC_TADs, "hESC_TADs")
#~ cisRandGR = lapply(cisRandGR, addWithinSubject, hESC_TADs, "hESC_TADs")
#~ 
#~ # same with Hi-C domains from Rao et al 2014
#~ cisParaGR = addWithinSubject(cisParaGR, raoDomDisjoin, "Rao_Domains")
#~ cisRandGR = lapply(cisRandGR, addWithinSubject, raoDomDisjoin, "Rao_Domains")

for(tadName in names(allTADs)){

    TAD = allTADs[[tadName]]
    message(paste("INFO: Compute overlap with TADs from:", tadName))

    # co-occurance within the same domain
    cisParaGR = addWithinSubject(cisParaGR, TAD, tadName)
    cisRandGR = bplapply(cisRandGR, addWithinSubject, TAD, tadName)
    Sys.sleep(3)

}


# Plot fraction of paralog pairs within TAD
for ( D in DIST_TH){
    
    pvalues = c()
    heightMat = c()
    sdMat = c()

    # repeat analysis for all TADs from Dixon et al and Rao et al    
    for( LABEL in names(allTADs) ){
        
        # get boolean vectors
        paraInTAD = mcols(cisParaGR)[width(cisParaGR)<=D, LABEL]
        randInTAD = lapply(cisRandGR, function(gr) mcols(gr)[width(gr)<=D, LABEL])

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

        pvalues = c(pvalues, fs.test$p.value)
        heightMat = cbind(heightMat, heights)
        sdMat = c(sdMat, sdRand)

        pdf(paste0(outPrefix, ".paralogPairs_cis_within_", D/10^3, "_kb_within_", LABEL, ".barplot.pdf"), width=3.5)
    
            par(cex=1.5, lwd=2)
            yMax = 1.3*max(heights) + sdRand
            
            bp = my.barplot(heights, addValues=TRUE, yMax=yMax,
                names=c("Paralogs", "Sampled pairs"),
                main=paste("Paralog pairs\n(within", D/10^3, "kb) in\n", LABEL, "TAD"), ylab="Gene pairs with TSSs in same TAD [%]", col=COL)

            error.bar(bp,heights, c(NA,  sdRand), lwd=2)
            
            # pvalue matrix for only two group comparision
            pval_mat = matrix(c(NA, signif(fs.test$p.value, 3)), 1)
            add_pval_all_allPairs(bp, ymax=1.1*max(heights)+sdRand, pval_mat)
    
        dev.off()
        
        #---------------------------------------------------------------
        # compare same strand
        #---------------------------------------------------------------
        # create contingency table and run Fisher test
        paraStrand = mcols(cisParaGR)[width(cisParaGR)<=D, "sameStrand"]
        randStrand = lapply(cisRandGR, function(gr) mcols(gr)[width(gr)<=D, "sameStrand"])

        contab = rbind(
            para=table(paraStrand),  
            rand=table(unlist(randStrand))
        )
        fs.test = fisher.test(contab)
        fs.test.inTAD = fisher.test(
           rbind(
            para=table(paraStrand & paraInTAD),  
            rand=table(unlist(randStrand) & unlist(randInTAD))) 
        )
        # get percent values for barplot
        paraPercent = percentTrue(paraStrand)
        randPercent = sapply(randStrand, percentTrue) 
        paraPercentBoth = percentTrue(paraStrand & paraInTAD)
        randPercentBoth = sapply(1:N_RAND, function(i) percentTrue(randInTAD[[i]] & randStrand[[i]]) ) 
        
        heights = matrix(c(paraPercent, mean(randPercent, na.rm=TRUE), paraPercentBoth, mean(randPercentBoth, na.rm=TRUE)), 2)
        sdRand = c(sd(randPercent, na.rm=TRUE), sd(randPercentBoth, na.rm=TRUE))
        
        # TODO make plot nicer
        pdf(paste0(outPrefix, ".paralogPairs_cis_within_", D/10^3, "_kb.SameStrand_within_", LABEL, ".barplot.pdf")) #, width=3.5
    
            par(cex=1.5, lwd=2)
            yMax = 1.4*max(heights) + max(sdRand, na.rm=TRUE)
            
            bp = my.barplot(heights, addValues=TRUE, yMax=yMax,
                names=c("Same strand", "Same strand\n and same TAD"),
                main=paste("Paralogs (within", D/10^3, "kb)\n with same strand\n in", LABEL, "TAD"), ylab="Gene pairs [%]", col=COL, srt=0, adj=NULL)

            error.bar(bp[2,], heights[2,], sdRand, lwd=2)
            
            add_pval_two_pairs(bp, values=1.1*heights, pvalues=c(fs.test$p.value, fs.test.inTAD$p.value))
            
            legend("topright", c("Paralogs", "Sampled"), fill=COL)

        dev.off()
        
        
    }
    
    # Plot for TADs from all cell-types together
    pdf(paste0(outPrefix, ".paralogPairs_cis_within_", D/10^3, "_kb_within_allTADs.barplot.pdf"), width=14)

        par(cex=1.5, lwd=2)
        bp = my.barplot(heightMat, addValues=TRUE, beside=TRUE, names=names(allTADs), col=COL, yMax=1.3*max(heightMat), digits=2, ylab="Gene pairs in same TAD [%]")
        error.bar(bp, heightMat, rbind(rep(NA, length(allTADs)), sdMat), length=0.05)

        add_pval_two_pairs(bp, heightMat, pvalues, offset=.15*max(heightMat)) #, digits=1

        legend("topleft", c("Paralogs", "Sampled pairs"), fill=COL, bty="n")
    dev.off()
    
}

#-----------------------------------------------------------------------
# Number of enhancers linked to each gene in the groups
#-----------------------------------------------------------------------

# show bias of number of enhancers linked to each gene:
paraEhPercent = percentCounts(tssGR[paralogs]$linked_enhancer)
nonParaEhPercent = percentCounts(tssGR[nonParalogs]$linked_enhancer)
mergedCounts = merge(paraEhPercent, nonParaEhPercent, all=TRUE, by="count")

# plot distribution of number of enhancers
pdf(paste0(outPrefix, ".number_of_enhancers.barplot.pdf"))
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
    legend("topright", legend=c("Paralog genes", "Non-Paralogs"), fill=COL2)
dev.off()

# plot fraction of pairs without and with at least one shared enhancers



# get enhancer-promotor distances 
tab = data.table(dist=abs(map$dist), gene=map$Gene)
enhancerDist = tab[,list(
    max=max(dist),
    mean=mean(dist),
    min=min(dist)
        ),by=gene]

parlogsDists = enhancerDist[enhancerDist$gene %in% paralogsHGNC,]
nonParlogsDists = enhancerDist[enhancerDist$gene %in% nonParalogsHGNC,]

# nearest enhancer
min.test = wilcox.test(parlogsDists$min, nonParlogsDists$min)

pdf(paste0(outPrefix, ".enhancer_dist_min.pdf"), 4.5, 7)
    par(cex=1.2)
    ymax = max(c(parlogsDists$min/10^3, nonParlogsDists$min/10^3))
    boxplot(parlogsDists$min/10^3, nonParlogsDists$min/10^3, 
        log="y", border=COL2, lwd=2, ylim=c(1, ymax+2*ymax),
        ylab="Distance to nearest enhancer (kb)")
    xlabels=paste(
        c("Paralogs", "Non-paralogs"), 
        "\n n =", 
        c(nrow(parlogsDists), nrow(nonParlogsDists)), 
        "\nMean =", 
        signif(c(mean(parlogsDists$min/10^3), mean(nonParlogsDists$min/10^3)),3)
        )
    axis(1, at = c(1, 2), labels=xlabels, line=1, tick=FALSE)
    pval_mat = matrix(c(NA, signif(min.test$p.value, 3)), 1)
    add_pval_all_allPairs(c(1,2), ymax=ymax+.5*ymax, pval_mat,offset_fac=.5)
dev.off()

#-----------------------------------------------------------------------
# Shared enhancer between paralog genes
#-----------------------------------------------------------------------

# iterate over distance threshold
for (D in DIST_TH){

    # subset of paralog pairs with distance < 500kb
    closeCisPairs = cisPairs[abs(cisPairs$dist) < D,]
    closeCisRand = lapply(randCisPairs, function(gp){gp[abs(gp$dist) < D,]})
    
    # build a data set with the frequency of occurrence numbers
    mergedCounts = merge(
        percentCounts(closeCisPairs$commonEnhancer), 
        percentCounts(rbind.fill(closeCisRand)[,"commonEnhancer"]), all=TRUE, by="count")    
    names(mergedCounts) = c("count", "paralogs", "random")
    
    # plot distribution of number of enhancers
    pdf(paste0(outPrefix, ".paralogPairs_cis_within_", D/10^3,"kb_common_enhancer.barplot.pdf"))
        maxCount = 5
        countTab = mergedCounts[seq(maxCount+1),]
        countTab[maxCount+1,] = countTab[maxCount+1,] +  colSums(mergedCounts[seq(maxCount+2, nrow(mergedCounts)),], na.rm=TRUE)
        height=t(countTab[,seq(2, ncol(countTab))])
        names=countTab[,1]
        names[maxCount+1] = paste0(maxCount, "+")
        par(cex=1.5)
        bp = barplot(height, names.arg=names,
            beside=TRUE, col=COL, ylim=c(0,100),
            main=paste("Shared enhancers among paralog pairs\n(within", D/10^3, "kb)"),
            xlab="Number of shared enhancers", ylab="%")
        text(bp, height, signif(height, 2), pos=3)
        legend("topright", legend=c("Paralog pairs", "Sampled pairs"), fill=COL)
    dev.off()

    # create contingency table and run Fisher test
    paraHasShared = closeCisPairs$commonEnhancer > 0
    randHasShared = lapply(closeCisRand, function(df){df$commonEnhancer > 0})
    contab = rbind(
        para=table(paraHasShared),  
        rand=table(unlist(randHasShared))
    )
    fs.test = fisher.test(contab)

    # plot fraction of pairs with at least one shared enhancer
    paraPercent = percentTrue(paraHasShared)
    randPercent = sapply(randHasShared, percentTrue) 
    heights = c(paraPercent, mean(randPercent, na.rm=TRUE))
    sdRand = sd(sapply(randHasShared, percentTrue), na.rm=TRUE)

    pdf(paste0(outPrefix, ".paralogPairs_cis_within_", D/10^3,"kb_has_shared_enhancer.barplot.pdf"), 3.5, 7)

        par(cex=1.5, lwd=2)
        ymax = max(heights) + sdRand + .2*max(heights)
        bp = my.barplot(heights, yMax=1.4*max(heights), addValues=TRUE,
            names=c("Paralog pairs", "Sampled pairs"),
            main=paste("Pairs (within", D/10^3, "kb) \n with\n shared enhancers"),
            ylab="Gene pairs with shared enhancer [%]", col=COL)
        error.bar(bp,heights, c(NA,  sdRand), lwd=2)
        
        # pvalue matrix for only two group comparision
        pval_mat = matrix(c(NA, signif(fs.test$p.value, 3)), 1)
        add_pval_all_allPairs(bp, ymax=1.1*max(heights)+sdRand, pval_mat)

    dev.off()
    
}

#-----------------------------------------------------------------------
# Hi-C analysis with data from Rao et al. 2014
#-----------------------------------------------------------------------


# Make data.frame for plotting with ggplot2
breaks = seq(DISTAL_MIN_DIST, DISTAL_MAX_DIST, length.out=10) / 10^3

plotDF = data.frame(
        group = c(rep("paralogs", nrow(distalCisPairs)), rep("sampled", nrow(randDistalCisPairsCombined))),
        HiCraw = c(distalCisPairs[,"HiCfreq"], randDistalCisPairsCombined[,"HiCfreq"]),
        HiCnorm = c(distalCisPairs[,"HiCnorm"], randDistalCisPairsCombined[,"HiCnorm"]),
        captureC_raw = c(distalCisPairs[,"captureC_raw"], randDistalCisPairsCombined[,"captureC_raw"]),
        captureC_ObsExp = c(distalCisPairs[,"captureC_ObsExp"], randDistalCisPairsCombined[,"captureC_ObsExp"]),
        dist=abs(c(distalCisPairs[,"dist"], randDistalCisPairsCombined[,"dist"]))/10^3,
        distBin=as.factor(breaks[.bincode(abs(c(distalCisPairs[,"dist"], randDistalCisPairsCombined[,"dist"]))/10^3, breaks)])
    )
    
plotDF$HiCrawNoZero = plotDF$HiCraw
plotDF$HiCrawNoZero[plotDF$HiCraw == 0] = NA
plotDF$HiCnormNoZero = plotDF$HiCnorm
plotDF$HiCnormNoZero[plotDF$HiCnorm == 0] = NA

p = ggplot(plotDF, aes(x=distBin, y=HiCraw)) + scale_y_log10() +
    geom_boxplot(aes(x=distBin, colour = group), lwd=1.5) + 
    scale_color_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + 
    guides(fill=guide_legend(title="")) +
    labs(y="Hi-C counts", x="Linear distance bin [kb]")

# save ggplot as pdf
ggsave(p, file=paste0(outPrefix, ".distal_pairs.Hi-C_raw_contacts_by_dist.boxplot.pdf"))

p = ggplot(plotDF, aes(x=distBin, y=HiCnorm)) + scale_y_log10() +
    geom_boxplot(aes(x=distBin, colour = group), lwd=1.5) + 
    scale_color_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + 
    guides(fill=guide_legend(title="")) +
    labs(y="Normalized Hi-C counts", x="Linear distance bin [kb]")

# save ggplot as pdf
ggsave(p, file=paste0(outPrefix, ".distal_pairs.Hi-C_normalized_contacts_by_dist.boxplot.pdf"))

# Plot capture C data between promoters in separate distance bins    
p = ggplot(plotDF, aes(x=distBin, y=captureC_raw)) + scale_y_log10() +
    geom_boxplot(aes(x=distBin, colour = group), lwd=1.5) + 
    scale_color_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + 
    guides(fill=guide_legend(title="")) +
    labs(y="Capture-Hi-C counts", x="Linear distance bin [kb]")

# save ggplot as pdf
ggsave(p, file=paste0(outPrefix, ".distal_pairs.captureC_raw_by_dist.boxplot.pdf"))

p = ggplot(plotDF, aes(x=distBin, y=captureC_ObsExp)) +
    geom_boxplot(aes(x=distBin, colour = group), lwd=1.5) + 
    scale_color_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + 
    guides(fill=guide_legend(title="")) +
    labs(y="Capture Hi-C (Observed/Expected)", x="Linear distance bin [kb]")

# save ggplot as pdf
ggsave(p, file=paste0(outPrefix, ".distal_pairs.captureC_ObsExp_by_dist.boxplot.pdf"))

#------------------------------------------------------------------------
# now for all pairs (not by dist)
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
    labs(y="Normalized Hi-C counts", x="", title=paste0("P = ", signif(ws.test$p.value, 3))) +
    scale_x_discrete(labels=xlabels )

# save ggplot as pdf
ggsave(p, file=paste0(outPrefix, ".distal_pairs.Hi-C_normalized_contacts.ggboxplot.pdf"), width=3.5)

#------------------------------------------------------------------------
# do the same without zero Hi-C counts
#------------------------------------------------------------------------
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
    labs(y="Normalized Hi-C counts", x="", title=paste0("P = ", signif(ws.test$p.value, 3))) +
    scale_x_discrete(labels=xlabels )

# save ggplot as pdf
ggsave(p, file=paste0(outPrefix, ".distal_pairs.Hi-C_normalized_contacts_noZero.ggboxplot.pdf"), width=3.5)

#------------------------------------------------------------------------
# for capture Hi-C
#------------------------------------------------------------------------
xlabels=paste0(c("Paralogs", "Sampled"), 
    "\n n=", applyToSubset(plotDF, function(v) sum(!is.na(v)), "captureC_raw", "group"), 
    "\n med=", signif(applyToSubset(plotDF, median, "captureC_raw", "group", na.rm=TRUE), 3),
    "\n avg=", signif(applyToSubset(plotDF, mean, "captureC_raw", "group", na.rm=TRUE), 3)
    )
ws.test = wilcox.test(plotDF[,"captureC_raw"] ~ plotDF[,"group"])

p = ggplot(plotDF, aes(x=group, y=captureC_raw, color=group)) + scale_y_log10() +
    geom_boxplot(lwd=1.5) + 
    scale_color_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), legend.position="none") + 
    guides(fill=guide_legend(title="")) +
    labs(y="Capture-Hi-C counts", x="", title=paste0("P = ", signif(ws.test$p.value, 3))) +
    scale_x_discrete(labels=xlabels )

# save ggplot as pdf
ggsave(p, file=paste0(outPrefix, ".distal_pairs.captureC_raw.ggboxplot.pdf"), width=3.5)
#------------------------------------------------------------------------
xlabels=paste0(c("Paralogs", "Sampled"), 
    "\n n=", applyToSubset(plotDF, function(v) sum(!is.na(v)), "captureC_ObsExp", "group"), 
    "\n med=", signif(applyToSubset(plotDF, median, "captureC_ObsExp", "group", na.rm=TRUE), 3),
    "\n avg=", signif(applyToSubset(plotDF, mean, "captureC_ObsExp", "group", na.rm=TRUE), 3)
    )
ws.test = wilcox.test(plotDF[,"captureC_ObsExp"] ~ plotDF[,"group"])

p = ggplot(plotDF, aes(x=group, y=captureC_ObsExp, color=group)) + scale_y_log10() +
    geom_boxplot(lwd=1.5) + 
    scale_color_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), legend.position="none") + 
    guides(fill=guide_legend(title="")) +
    labs(y="Capture Hi-C (Observed/Expected)", x="", title=paste0("P = ", signif(ws.test$p.value, 3))) +
    scale_x_discrete(labels=xlabels )

# save ggplot as pdf
ggsave(p, file=paste0(outPrefix, ".distal_pairs.captureC_ObsExp.ggboxplot.pdf"), width=3.5)

#------------------------------------------------------------------------
# plot interactions for close cis pairs (within MAX_DIST)
#------------------------------------------------------------------------
plotDFclose = data.frame(
        group = c(rep("paralogs", nrow(cisPairs)), rep("sampled", nrow(randCisPairsCombined))),
        captureC_raw = c(cisPairs[,"captureC_raw"], randCisPairsCombined[,"captureC_raw"]),
        captureC_ObsExp = c(cisPairs[,"captureC_ObsExp"], randCisPairsCombined[,"captureC_ObsExp"])
    )

xlabels=paste0(c("Paralogs", "Sampled"), 
    "\n n=", applyToSubset(plotDFclose, function(v) sum(!is.na(v)), "captureC_raw", "group"), 
    "\n med=", signif(applyToSubset(plotDFclose, median, "captureC_raw", "group", na.rm=TRUE), 3),
    "\n avg=", signif(applyToSubset(plotDFclose, mean, "captureC_raw", "group", na.rm=TRUE), 3)
    )
ws.test = wilcox.test(plotDFclose[,"captureC_raw"] ~ plotDFclose[,"group"])

p = ggplot(plotDFclose, aes(x=group, y=captureC_raw, color=group)) + scale_y_log10() +
    geom_boxplot(lwd=1.5) + 
    scale_color_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), legend.position="none") + 
    guides(fill=guide_legend(title="")) +
    labs(y="Capture-Hi-C counts", x="", title=paste0("P = ", signif(ws.test$p.value, 3))) +
    scale_x_discrete(labels=xlabels )

# save ggplot as pdf
ggsave(p, file=paste0(outPrefix, ".close_pairs.captureC_raw.ggboxplot.pdf"), width=3.5)
#------------------------------------------------------------------------
xlabels=paste0(c("Paralogs", "Sampled"), 
    "\n n=", applyToSubset(plotDFclose, function(v) sum(!is.na(v)), "captureC_ObsExp", "group"), 
    "\n med=", signif(applyToSubset(plotDFclose, median, "captureC_ObsExp", "group", na.rm=TRUE), 3),
    "\n avg=", signif(applyToSubset(plotDFclose, mean, "captureC_ObsExp", "group", na.rm=TRUE), 3)
    )
ws.test = wilcox.test(plotDFclose[,"captureC_ObsExp"] ~ plotDFclose[,"group"])

p = ggplot(plotDFclose, aes(x=group, y=captureC_ObsExp, color=group)) + scale_y_log10() +
    geom_boxplot(lwd=1.5) + 
    scale_color_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=20), legend.position="none") + 
    guides(fill=guide_legend(title="")) +
    labs(y="Capture Hi-C (Observed/Expected)", x="", title=paste0("P = ", signif(ws.test$p.value, 3))) +
    scale_x_discrete(labels=xlabels )

# save ggplot as pdf
ggsave(p, file=paste0(outPrefix, ".close_pairs.captureC_ObsExp.ggboxplot.pdf"), width=3.5)


#-----------------------------------------------------------------------
# Expression analysis
#-----------------------------------------------------------------------

# combine all sampled data sets:
randCisPairsAll = do.call("rbind", randCisPairs)

paraExpCorSummaryList=list()

# iterate over all expression data sets
for (expName in names(expDFlist)){

    expDF = expDFlist[[expName]]

    # plot distribution of correlation coefficient
    paraExpCor = cisPairs[,paste0(expName, "_expCor")]
    sampledExpCor = randCisPairsAll[,paste0(expName, "_expCor")]
    
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
    commonEhGroups = getGroupsByCounts(cisPairs$commonEnhancer, maxCount)
    
    # make ggplot 
    p = ggplot(cisPairs, aes(x=commonEhGroups, y=paraExpCor)) + geom_boxplot(fill=NA, aes(colour = commonEhGroups), lwd=1.5) + geom_jitter(alpha=.3, size=2, aes(colour = commonEhGroups)) + scale_colour_manual(values= colorRampPalette(c("gray", COL_EH))(maxCount+1), guide=FALSE) + theme_bw() + theme(text = element_text(size=20)) + labs(y="Pearson correlation coefficient", x = "Number of shared enhancers")
    
    # save ggplot as pdf
    ggsave(p, file=paste0(outPrefix,".expCor_by_shared_enhancers.", expName, ".boxplot.pdf"), width=7, height=7)

    #-------------------------------------------------------------------
    # Compare expression against localization in common TAD
    #-------------------------------------------------------------------
    for (tadName in names(allTADs)){
        inTADgroup <- factor(mcols(cisParaGR)[, tadName], levels=c(FALSE, TRUE), labels=c("Not in same TAD", "In same TAD"))
        
        # test difference with wilcoxon test
        ws.test = wilcox.test(paraExpCor ~  inTADgroup)
    
        p = ggplot(cisPairs, aes(x=inTADgroup, y=paraExpCor)) + geom_boxplot(fill=NA, aes(colour = inTADgroup), lwd=1.5) + geom_jitter(alpha=.3, size=2, aes(colour = inTADgroup)) + scale_colour_manual(values= colorRampPalette(c("gray", COL_TAD[1]))(2), guide=FALSE) + theme_bw() + theme(text = element_text(size=20), axis.text.x = element_text(angle = 45, hjust = 1)) + labs(y="Pearson correlation coefficient", x = "", title=paste0("p = ",signif(ws.test$p.value,2)))
        # save ggplot as pdf
        ggsave(p, file=paste0(outPrefix,".expCor_", expName, ".by_sahred_TAD_", tadName, ".boxplot.pdf"), width=3.5, height=7)
    }

    # plot for all cell types: 
    expPlotDF = data.frame(
        expCor = rep(paraExpCor, length(allTADs)),
        withinTAD = unlist(lapply(names(allTADs), function(tadName) factor(mcols(cisParaGR)[, tadName], levels=c(TRUE, FALSE), labels=c("In same TAD", "Not in same TAD")))),
        cells = rep(factor(gsub("_", " ", names(allTADs)), gsub("_", " ", names(allTADs))) , each=length(paraExpCor))
    )
    p = ggplot(expPlotDF, aes(x=withinTAD, y=expCor, colour=withinTAD)) + geom_violin(aes(colour = withinTAD), lwd=1, alpha=.25) + geom_boxplot(aes(color=withinTAD), fill=NA, width=.25, lwd=1)  + 
    facet_wrap(~ cells) + scale_colour_manual(values= colorRampPalette(c(COL_TAD[1], "gray"))(2), guide=FALSE) + theme_bw() + theme(text = element_text(size=20), axis.text.x = element_text(angle = 45, hjust = 1)) + labs(y="Pearson correlation coefficient", x = "") +  theme(legend.position="bottom")

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

#-----------------------------------------------------------------------
# LRRC8 Gene family
#-----------------------------------------------------------------------

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
# 5.) Compare one2one orthologs of human paralogs in mouse and dog
#=======================================================================
#allCisPairs
#randPairsInCis


# iterate over all species:
for (orgStr in c("mmusculus", "cfamiliaris")){
    
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
            ylab="Pairs on some chrom [%]")        
        error.bar(bp,height, c(NA, sd(randPercent)))
    dev.off()
    
    # compare linear distance between orthologs of human paralogs and random genes
    orthoDistStr = paste0(orgStr, "_dist")
    orthoDist = abs(allCisPairs[abs(allCisPairs[, orthoDistStr]) <= MAX_DIST, orthoDistStr]) /10^3
    randDist = abs(unlist(
        lapply(randPairsInCis, function(d){
            d[abs(d[, orthoDistStr]) <= MAX_DIST, orthoDistStr]
            })
        )) / 10^3
    
    pdf(paste0(outPrefix, ".paralogPairs_orthologs_dist_", orgName, ".hist.pdf"))
        par(lwd=2, cex=1.5, mfrow=c(2,1))

        hist(orthoDist, 50, col=COL_ORTHO[1],
            main="Disance between orthologs of paralogs", xlab="Distance (kb)")
        hist(randDist, 50, col=COL_RAND[2],
            main="Distance between orthologs of random gene pairs", xlab="Distance (kb)")    
    dev.off()
    
    # compare distance to dist-wighted sampled pairs    
    orthoDistStr = paste0(orgStr, "_dist")
    
    # take only those pairs into acount that have one-to-one orthologs on same chrom within 1 MB.
    caseSubset = !is.na(cisPairs[, orthoDistStr]) & abs(cisPairs[, orthoDistStr]) <= MAX_DIST
    randSubset = lapply(randCisPairs, function(gP) !is.na(gP[, orthoDistStr]) & abs(gP[, orthoDistStr]) <= MAX_DIST)
    
    # get distances of four sets
    paraDist = abs(cisPairs[caseSubset,"dist"]) /10^3
    randDist = abs(unlist(lapply(1:N_RAND, function(i) randCisPairs[[i]][randSubset[[i]],"dist"]))) / 10^3

    orthoDist = abs(cisPairs[caseSubset, orthoDistStr]) /10^3
    randOrthoDist = abs(unlist(
            lapply(1:N_RAND, function(i) randCisPairs[[i]][randSubset[[i]], orthoDistStr]) 
        )) / 10^3
        
    pdf(paste0(outPrefix, ".paralogPairs_orthologs_dist_sampled_", orgName, ".hist.pdf"))
        par(lwd=2, cex=1.5, mfrow=c(2,2))

        hist(paraDist, 50, col=COL[1],
            main="Disance between\n paralogs", xlab="Distance (kb)")
        hist(randDist, 50, col=COL[2],
            main="Distance between\n sampled gene pairs", xlab="Distance (kb)")    
        hist(orthoDist, 50, col=COL_ORTHO[1],
            main="Disance between\n orthologs of paralogs", xlab="Distance (kb)")
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
    orthoInTAD = cisPairs[cisPairs[,paste0(orgStr, "_one2one")], paste0(orgStr, "_TAD")]
    randInTAD = lapply(randCisPairs, function(gP) gP[gP[,paste0(orgStr, "_one2one")], paste0(orgStr, "_TAD")])

    # create contingency table and run Fisher test
    contab = rbind(
        ortho=table(orthoInTAD),  
        rand=table(unlist(randInTAD))
    )
    fs.test = fisher.test(contab)
        
    pdf(paste0(outPrefix, ".paralogPairs_orthologs_inSameTAD", orgName, ".barplot.pdf"), width=3.5)

        paraPercent = percentTrue(orthoInTAD)
        randPercent = sapply(randInTAD, percentTrue)        
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

        paraInTAD <- mcols(cisParaGR)[, tadName]
        randInTAD = lapply(cisRandGR, function(gr) mcols(gr)[, tadName])

        humanAndOrthoInTAD = cisPairs[cisPairs[,paste0(orgStr, "_one2one")] & paraInTAD, paste0(orgStr, "_TAD")]
        humanAndRandInTAD = lapply(1:length(randCisPairs), function(i){ 
            gP = randCisPairs[[i]]
            gP[gP[,paste0(orgStr, "_one2one")] & randInTAD[[i]], paste0(orgStr, "_TAD")]
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
    # Compare Hi-C counts of orthologs
    
    # combine all random pairs to one DF
    randCisPairsAll = do.call("rbind", randCisPairs)    
    plotDF = data.frame(
        group = c(rep("paralogs", nrow(cisPairs)), rep("sampled", nrow(randCisPairsAll))),
        HiCraw = c(cisPairs[,paste0(orgStr, "_HiC")], randCisPairsAll[,paste0(orgStr, "_HiC")]),
        HiCnorm = c(cisPairs[,paste0(orgStr, "_HiCnorm")], randCisPairsAll[,paste0(orgStr, "_HiCnorm")])
    )
    
    #------------------------------------------------------------------------
    xlabels=paste0(c("Paralogs", "Sampled"), 
        "\n n=", applyToSubset(plotDF, function(v) sum(!is.na(v)), "HiCraw", "group"), 
        "\n med=", signif(applyToSubset(plotDF, median, "HiCraw", "group", na.rm=TRUE), 3),
        "\n avg=", signif(applyToSubset(plotDF, mean, "HiCraw", "group", na.rm=TRUE), 3)
        )
    ws.test = wilcox.test(HiCnorm ~ group, data=plotDF)
    p = ggplot(plotDF, aes(x=group, y=HiCraw)) + scale_y_log10() +
        #geom_jitter(alpha=.3, size=2, aes(colour = group), guide=FALSE) + 
        geom_boxplot(fill=NA, aes(colour = group), lwd=1.5) + 
        scale_color_manual(values=COL_ORTHO, guide=FALSE) +
        theme_bw() + theme(text = element_text(size=20)) + 
        labs(y=paste("Hi-C counts in", orgName), x="", title=paste0("P = ", signif(ws.test$p.value, 3))) + scale_x_discrete(labels=xlabels)
    
    # save ggplot as pdf
    ggsave(p, file=paste0(outPrefix, ".paralogPairs_orthologs_HiCraw", orgName, ".boxplot.pdf"), width=3.5, height=7)
    
    # Normalized Hi-C counts
    #------------------------------------------------------------------------
    xlabels=paste0(c("Paralogs", "Sampled"), 
        "\n n=", applyToSubset(plotDF, function(v) sum(!is.na(v)), "HiCnorm", "group"), 
        "\n med=", signif(applyToSubset(plotDF, median, "HiCnorm", "group", na.rm=TRUE), 3),
        "\n avg=", signif(applyToSubset(plotDF, mean, "HiCnorm", "group", na.rm=TRUE), 3)
        )
    ws.test = wilcox.test(HiCnorm ~ group, data=plotDF)
    p = ggplot(plotDF, aes(x=group, y=HiCnorm)) + scale_y_log10() + 
        #geom_jitter(alpha=.3, size=2, aes(colour = group), guide=FALSE) + 
        geom_boxplot(fill=NA, aes(colour = group), lwd=1.5) + 
        scale_color_manual(values=COL_ORTHO, guide=FALSE) +
        theme_bw() + theme(text = element_text(size=20)) + 
        labs(y=paste("log2(observed / expected) Hi-C counts in", orgName), x="", title=paste0("P = ", signif(ws.test$p.value, 3))) + scale_x_discrete(labels=xlabels)
    
    # save ggplot as pdf
    ggsave(p, file=paste0(outPrefix, ".paralogPairs_orthologs_HiCnorm", orgName, ".boxplot.pdf"), width=3.5, height=7)

    #-------------------------------------------------------------------
    # do the same for distal cis pairs
    #-------------------------------------------------------------------
    # combine all random pairs to one DF
    randDistalCisPairsAll = do.call("rbind", randDistalCisPairs)    
    plotDF = data.frame(
        group = c(rep("paralogs", nrow(distalCisPairs)), rep("sampled", nrow(randDistalCisPairsAll))),
        HiCraw = c(distalCisPairs[,paste0(orgStr, "_HiC")], randDistalCisPairsAll[,paste0(orgStr, "_HiC")]),
        HiCnorm = c(distalCisPairs[,paste0(orgStr, "_HiCnorm")], randDistalCisPairsAll[,paste0(orgStr, "_HiCnorm")])
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
    ws.test = wilcox.test(HiCraw ~ group, data=plotDF)
    
    p = ggplot(plotDF, aes(x=group, y=HiCraw)) + scale_y_log10() + 
        #geom_jitter(alpha=.3, size=2, aes(colour = group), guide=FALSE) + 
        geom_boxplot(fill=NA, aes(colour = group), lwd=1.5) + 
        scale_color_manual(values=COL_ORTHO, guide=FALSE) +
        theme_bw() + theme(text = element_text(size=20)) + 
        labs(y=paste("Hi-C counts in ", orgName), x="", title=paste0("P = ", signif(ws.test$p.value, 3))) + scale_x_discrete(labels=xlabels)
    
    # save ggplot as pdf
    ggsave(p, file=paste0(outPrefix, ".paralogPairs_orthologs.distalPairs_HiCraw.", orgName, ".boxplot.pdf"), width=3.5, height=7)
    
    #------------------------------------------------------------------------
    # Normalized Hi-C counts
    xlabels=paste0(c("Paralogs", "Sampled"), 
        "\n n=", applyToSubset(plotDF, function(v) sum(!is.na(v)), "HiCnorm", "group"), 
        "\n med=", signif(applyToSubset(plotDF, median, "HiCnorm", "group", na.rm=TRUE), 3),
        "\n avg=", signif(applyToSubset(plotDF, mean, "HiCnorm", "group", na.rm=TRUE), 3)
        )    
    ws.test = wilcox.test(HiCnorm ~ group, data=plotDF)
    p = ggplot(plotDF, aes(x=group, y=HiCnorm)) + scale_y_log10() + 
        #geom_jitter(alpha=.3, size=2, aes(colour = group), guide=FALSE) + 
        geom_boxplot(fill=NA, aes(colour = group), lwd=1.5) + 
        scale_color_manual(values=COL_ORTHO, guide=FALSE) +
        theme_bw() + theme(text = element_text(size=20)) + 
        labs(y=paste("Normalized Hi-C counts in ", orgName), x="", title=paste0("P = ", signif(ws.test$p.value, 3))) + scale_x_discrete(labels=xlabels)
    
    # save ggplot as pdf
    ggsave(p, file=paste0(outPrefix, ".paralogPairs_orthologs.distalPairs_HiCnorm.", orgName, ".boxplot.pdf"), width=3.5, height=7)
    
    #------------------------------------------------------------------------
    # No zero raw counts
    xlabels=paste0(c("Paralogs", "Sampled"), 
        "\n n=", applyToSubset(plotDF, function(v) sum(!is.na(v)), "HiCrawNoZero", "group"), 
        "\n med=", signif(applyToSubset(plotDF, median, "HiCrawNoZero", "group", na.rm=TRUE), 3),
        "\n avg=", signif(applyToSubset(plotDF, mean, "HiCrawNoZero", "group", na.rm=TRUE), 3)
        )
    ws.test = wilcox.test(HiCrawNoZero ~ group, data=plotDF)
    
    p = ggplot(plotDF, aes(x=group, y=HiCrawNoZero)) + scale_y_log10() + 
        #geom_jitter(alpha=.3, size=2, aes(colour = group), guide=FALSE) + 
        geom_boxplot(fill=NA, aes(colour = group), lwd=1.5) + 
        scale_color_manual(values=COL_ORTHO, guide=FALSE) +
        theme_bw() + theme(text = element_text(size=20)) + 
        labs(y=paste("Hi-C counts in ", orgName), x="", title=paste0("P = ", signif(ws.test$p.value, 3))) + scale_x_discrete(labels=xlabels)
    
    # save ggplot as pdf
    ggsave(p, file=paste0(outPrefix, ".paralogPairs_orthologs.distalPairs_HiCrawNoZero.", orgName, ".boxplot.pdf"), width=3.5, height=7)    
    #------------------------------------------------------------------------
    # Normalized Hi-C counts without Zero counts
    xlabels=paste0(c("Paralogs", "Sampled"), 
        "\n n=", applyToSubset(plotDF, function(v) sum(!is.na(v)), "HiCnormNoZero", "group"), 
        "\n med=", signif(applyToSubset(plotDF, median, "HiCnormNoZero", "group", na.rm=TRUE), 3),
        "\n avg=", signif(applyToSubset(plotDF, mean, "HiCnormNoZero", "group", na.rm=TRUE), 3)
        )    
    ws.test = wilcox.test(HiCnormNoZero ~ group, data=plotDF)
    p = ggplot(plotDF, aes(x=group, y=HiCnormNoZero)) + scale_y_log10() + 
        #geom_jitter(alpha=.3, size=2, aes(colour = group), guide=FALSE) + 
        geom_boxplot(fill=NA, aes(colour = group), lwd=1.5) + 
        scale_color_manual(values=COL_ORTHO, guide=FALSE) +
        theme_bw() + theme(text = element_text(size=20)) + 
        labs(y=paste("Normalized Hi-C counts in ", orgName), x="", title=paste0("P = ", signif(ws.test$p.value, 3))) + scale_x_discrete(labels=xlabels)
    
    # save ggplot as pdf
    ggsave(p, file=paste0(outPrefix, ".paralogPairs_orthologs.distalPairs_HiCnormNoZero.", orgName, ".boxplot.pdf"), width=3.5, height=7)
    
}

#=======================================================================
# 6.) Mouse and dog paralog analysis
#=======================================================================
runBasicParalogAnalysis(paste0(outPrefix, ".Mouse"), paralogPairsMouse, "mmusculus_paralog_ds", tssGRmouse, speciesTADs[["mmusculus"]], tissueName="Mouse", HiClist=speciesHiC[["mmusculus"]][[1]], HiClistNorm=speciesHiC[["mmusculus"]][[2]])
message("Finish Mouse paralog analysis!")

runBasicParalogAnalysis(paste0(outPrefix, ".Dog"), paralogPairsDog, "cfamiliaris_paralog_ds", tssGRdog, speciesTADs[["cfamiliaris"]], tissueName="Dog", HiClist=speciesHiC[["cfamiliaris"]][[1]], HiClistNorm=speciesHiC[["cfamiliaris"]][[2]])
message("Finish Dog paralog analysis!")

#=======================================================================
# Compare sahred enhancer positions
#=======================================================================
nEhPos <- 3
ehColNames = c("eh_up", "eh_cent", "eh_down")
ehPosFactor <- factor(c("upstream", "center", "downstream"), levels=c("upstream", "center", "downstream"), labels=c("upstream", "center", "downstream"))

ehPosDF = data.frame(
        "group" = c(rep("paralogs", nParaPairs*nEhPos), rep("sampled", nRandPairs*nEhPos)),
        "enahncer_position"= unlist(list(rep(ehPosFactor, each=nParaPairs), rep(ehPosFactor, each=nRandPairs))),
        "eh_counts"=c(unlist(cisPairs[,ehColNames]), unlist(randCisPairsCombined[,ehColNames])),
        "sameIMR90_TAD"=factor(c(rep(mcols(cisParaGR)[,"Rao_IMR90"], nEhPos), rep(mcols(cisRandGRall)[,"Rao_IMR90"], nEhPos)), levels=c(TRUE, FALSE), labels=c("same TAD (Rao IMR90)", "not same TAD"))
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
        theme_bw() + facet_grid(group~sameIMR90_TAD, scales="free_y") + 
        labs(y="Summed counts", x="Enhancer position", title= "Shared enhancer position relative to genes pairs") + 
        scale_fill_manual(values=colorRampPalette(brewer.pal(9, "Oranges"))(length(unique(ehPosDF$eh_count))), guide = guide_legend("Enhancer\ncounts")) 
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

# combine all sampled data sets:
cisRandGRall <- do.call("c", cisRandGR)
cisRandAll <- mcols(cisRandGRall)

nSpecies=length(orgStr2Name) 
nTAD <- length(allTADs)
nParaPairs <- nrow(cisPairs)
nRandPairs <- sum(sapply(cisRandGR, length))


# check if "one2one ortholog" is a good approximation for duplication age
seqSimDF <- data.frame(
    "species"= rep(factor(orgStr2Name, orgStr2Name), each=nParaPairs),
    "one2one" = c(cisPairs[,"mmusculus_one2one"], cisPairs[,"cfamiliaris_one2one"]),
    "Ds"= rep(cisPairs[,"hsapiens_paralog_ds"],nSpecies),
    "Dn"= rep(cisPairs[,"hsapiens_paralog_dn"],nSpecies)
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
    dist = rep(c(cisPairs[,"dist"], randCisPairsCombined[,"dist"]), nExp),
    commonEnhancer = rep(c(cisPairs[,"commonEnhancer"], randCisPairsCombined[,"commonEnhancer"]), nExp),
    DupAge_beforeMouse = factor(rep(c(cisPairs[,"mmusculus_one2one"], randCisPairsCombined[,"mmusculus_one2one"]), nExp), c(FALSE, TRUE), c("Young", "Old")),
    DupAge_beforeDog = factor(rep(c(cisPairs[,"cfamiliaris_one2one"], randCisPairsCombined[,"cfamiliaris_one2one"]), nExp), c(FALSE, TRUE), c("Young", "Old"))
)

pdf(paste0(outPrefix,"old_vs_jung.beforeMouse.Expression_correlation.all_data.boxplot.pdf"))
    ggplot(paraExpCorDF, na.rm=TRUE, aes(x=genepairs, y = expCor, fill=DupAge_beforeMouse)) + theme_bw() + labs(y="Pearson correlation coefficient", title= "Co-expression correlation over tissues")  + facet_grid(.~dataset) + theme(legend.position = "bottom") + theme(text = element_text(size=15)) + scale_color_manual(values=COL) + scale_fill_manual(values=COL_AGE) + geom_boxplot(lwd=1)
dev.off()

pdf(paste0(outPrefix,"old_vs_jung.beforeMouse.Expression_correlation.all_data.violinplot.pdf"))
    ggplot(paraExpCorDF, na.rm=TRUE, aes(x=genepairs, y = expCor, fill=DupAge_beforeMouse)) + theme_bw() + labs(y="Pearson correlation coefficient", title= "Co-expression correlation over tissues")  + facet_grid(.~dataset) + theme(legend.position = "bottom") + theme(text = element_text(size=15)) + scale_color_manual(values=COL) + scale_fill_manual(values=COL_AGE) + geom_violin(adjust = .2)
dev.off()

pdf(paste0(outPrefix,"old_vs_jung.beforeDog.Expression_correlation.all_data.boxplot.pdf"))
    ggplot(paraExpCorDF, na.rm=TRUE, aes(x=genepairs, y = expCor, fill=DupAge_beforeDog)) + theme_bw() + labs(y="Pearson correlation coefficient", title= "Co-expression correlation over tissues")  + facet_grid(.~dataset) + theme(legend.position = "bottom") + theme(text = element_text(size=15)) + scale_color_manual(values=COL) + scale_fill_manual(values=COL_AGE) + geom_boxplot(lwd=1)
dev.off()

pdf(paste0(outPrefix,"old_vs_jung.beforeMouse.sharedEnhancer.pdf"))
#~     ggplot(paraExpCorDF[1:(nrow(paraExpCorDF)/nExp),], aes(commonEnhancer, group=DupAge_beforeMouse, fill=DupAge_beforeMouse)) + 
#~     geom_histogram(aes(y = ..density..), breaks=0:8, position="dodge") + xlim(0,8) +
#~     theme_bw() + 
#~     labs(y="Number of pairs", title= "Number of shared enhancer")  + theme(legend.position = "bottom") + scale_color_manual(values=COL) + scale_fill_manual(values=COL_AGE) + 
#~     facet_grid(genepairs~.)
    ggplot(paraExpCorDF[1:(nrow(paraExpCorDF)/nExp),], aes(commonEnhancer, group=DupAge_beforeMouse, fill=DupAge_beforeMouse)) + 
    geom_bar(aes(y = ..density..), badwidth=1, width=.5, position="dodge") + xlim(0,10) +
    theme_bw() + 
    labs(y="Density", title= "Number of shared enhancer")  + theme(legend.position = "bottom") + scale_color_manual(values=COL) + scale_fill_manual(values=COL_AGE) + 
    facet_grid(genepairs~.)
dev.off()


pdf(paste0(outPrefix,"old_vs_jung.dist.pdf"))
    ggplot(paraExpCorDF[1:(nrow(paraExpCorDF)/nExp),], aes(abs(dist)/10^3, group=DupAge_beforeMouse, colour=DupAge_beforeMouse, fill=DupAge_beforeMouse)) + 
#~     scale_x_log10() +
#~     geom_histogram() +
    geom_density(alpha = 0.1, lwd=1.2) +
    theme_bw() + 
    labs(y="Density", x="Distance [kb]")  + theme(legend.position = "bottom") + scale_color_manual(values=COL_AGE) + scale_fill_manual(values=COL_AGE) +
    facet_grid(genepairs~.)
dev.off()


tadDF <- data.frame(
    "sameTAD" = c(
            unlist(rbind(mcols(cisParaGR)[,names(allTADs)])),
            unlist(rbind(mcols(cisRandGRall)[,names(allTADs)]))
                ),
    "genepairs" = c(rep("Paralog", nParaPairs*nTAD), rep("Sampled", nRandPairs*nTAD)),    
    "cellType" = c(rep(names(allTADs), each=nParaPairs), rep(names(allTADs), each=nRandPairs)),    
#~     DupAge_beforeMouse = factor(rep(c(cisPairs[,"mmusculus_one2one"], randCisPairsCombined[,"mmusculus_one2one"]), nTAD), c(FALSE, TRUE), c("Young", "Old")),
    DupAge_beforeMouse = factor(c(rep(cisPairs[,"mmusculus_one2one"], nTAD), rep(randCisPairsCombined[,"mmusculus_one2one"], nTAD)), c(FALSE, TRUE), c("Young", "Old")),
    DupAge_beforeDog = factor(c(rep(cisPairs[,"cfamiliaris_one2one"], nTAD), rep(randCisPairsCombined[,"cfamiliaris_one2one"], nTAD)), c(FALSE, TRUE), c("Young", "Old"))
)

# get the fraction for each combination
fracTadDF <- ddply(tadDF, .(genepairs, cellType, DupAge_beforeMouse), summarize, frac=percentTrue(sameTAD))

pdf(paste0(outPrefix,"old_vs_jung.sameTAD.pdf"))

    ggplot(fracTadDF, aes(y=frac, x=cellType, fill=DupAge_beforeMouse)) + 
    geom_bar(stat="identity", position="dodge") +
    theme_bw() + 
    labs(y="Percent of pairs in same TAD", x="")  + theme(legend.position = "bottom") +  scale_fill_manual(values=COL_AGE) +
    geom_text(aes(label=signif(frac,3)), position=position_dodge(width=1), vjust=-0.25, size=3) + 
    facet_grid(genepairs~.) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

#-----------------------------------------------------------------------
nParaDist <- nrow(distalCisPairs)
nRandDist <- nrow(randDistalCisPairsCombined)

hicDF = data.frame(
        group = c(rep("paralogs", nParaDist), rep("sampled", nRandDist)),
        DupAge_beforeMouse = factor(c(distalCisPairs[,"mmusculus_one2one"], randDistalCisPairsCombined[,"mmusculus_one2one"]), c(FALSE, TRUE), c("Young", "Old")),
        DupAge_beforeDog = factor(c(distalCisPairs[,"cfamiliaris_one2one"], randDistalCisPairsCombined[,"cfamiliaris_one2one"]), c(FALSE, TRUE), c("Young", "Old")),
        HiCraw = c(distalCisPairs[,"HiCfreq"], randDistalCisPairsCombined[,"HiCfreq"]),
        HiCnorm = c(distalCisPairs[,"HiCnorm"], randDistalCisPairsCombined[,"HiCnorm"]),
        captureC_raw = c(distalCisPairs[,"captureC_raw"], randDistalCisPairsCombined[,"captureC_raw"]),
        captureC_ObsExp = c(distalCisPairs[,"captureC_ObsExp"], randDistalCisPairsCombined[,"captureC_ObsExp"]),
        dist=abs(c(distalCisPairs[,"dist"], randDistalCisPairsCombined[,"dist"]))/10^3,
        distBin=as.factor(breaks[.bincode(abs(c(distalCisPairs[,"dist"], randDistalCisPairsCombined[,"dist"]))/10^3, breaks)])
)
hicDF$HiCrawNoZero = hicDF$HiCraw
hicDF$HiCrawNoZero[hicDF$HiCraw == 0] = NA
hicDF$HiCnormNoZero = hicDF$HiCnorm
hicDF$HiCnormNoZero[hicDF$HiCnorm == 0] = NA

pdf(paste0(outPrefix,"old_vs_jung.Hi-C.pdf"), w=3.5, h=7)
    ggplot(hicDF, aes(x=group, y=HiCraw, colour = DupAge_beforeMouse))  +
    geom_boxplot(lwd=1.5) + scale_y_log10() + 
    scale_color_manual(values=COL_AGE, name="") +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + 
    guides(fill=guide_legend(title="")) +
    labs(y="Hi-C counts", x="")
dev.off()

pdf(paste0(outPrefix,"old_vs_jung.Hi-C_norm.pdf"), w=3.5, h=7)
    ggplot(hicDF, aes(x=group, y=HiCnorm, colour = DupAge_beforeMouse))  +
    geom_boxplot(lwd=1.5) + scale_y_log10() + 
    scale_color_manual(values=COL_AGE, name="") +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + 
    guides(fill=guide_legend(title="")) +
    labs(y="Normalized Hi-C counts", x="")
dev.off()

pdf(paste0(outPrefix,"old_vs_jung.CaptureC_raw.pdf"), w=3.5, h=7)
    ggplot(hicDF, aes(x=group, y=captureC_raw, colour = DupAge_beforeMouse))  +
    geom_boxplot(lwd=1.5) + scale_y_log10() + 
    scale_color_manual(values=COL_AGE, name="") +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + 
    guides(fill=guide_legend(title="")) +
    labs(y="Capture-C counts", x="")
dev.off()

pdf(paste0(outPrefix,"old_vs_jung.captureC_ObsExp.pdf"), w=3.5, h=7)
    ggplot(hicDF, aes(x=group, y=captureC_ObsExp, colour = DupAge_beforeMouse))  +
    geom_boxplot(lwd=1.5) + scale_y_log10() + 
    scale_color_manual(values=COL_AGE, name="") +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + 
    guides(fill=guide_legend(title="")) +
    labs(y="Capture-C (Obs/Exp)", x="")
dev.off()

#=======================================================================
# Best candidate picking...
#=======================================================================
#~ candidates = with(cisPairs, commonEnhancer>3 & sameStrand & abs(ENCODE_cell_lines_expCor) > .6 ) & cisParaGR$stable_TADs 
#~ table(candidates) 

#cisPairs[candidates,]

# write the data of Fetuin-A and Fetuin-B (ENSG00000145192, ENSG00000090512)
#grep("ENSG00000145192", cisPairs[,2])

write.table(t(as.matrix(mcols(cisParaGR[515,]))), file=paste0(outPrefix, ".Fetuin_AB.annotation.txt"), sep="\t", quote=FALSE, col.names=FALSE)

#=======================================================================
# Olfactory receptors:
#=======================================================================
orPairIDX <- which(cisPairs[,1] %in% ORids & cisPairs[,2] %in% ORids)
orDistalIDX <- which(distalCisPairs[,1] %in% ORids & distalCisPairs[,2] %in% ORids)


# get subset of pairs taht are OR
orPairs <- cisPairs[orPairIDX,]
orPairGR <- cisParaGR[orPairIDX]

orDistalPairs <- distalCisPairs[orDistalIDX,]

# Make data.frame for plotting with ggplot2
breaks = seq(DISTAL_MIN_DIST, DISTAL_MAX_DIST, length.out=10) / 10^3
plotDF = data.frame(
        group = c(rep("OR", nrow(orDistalPairs)), rep("paralogs", nrow(distalCisPairs)), rep("sampled", nrow(randDistalCisPairsCombined))),
        HiCraw = c(orDistalPairs[,"HiCfreq"], distalCisPairs[,"HiCfreq"], randDistalCisPairsCombined[,"HiCfreq"]),
        HiCnorm = c(orDistalPairs[,"HiCnorm"], distalCisPairs[,"HiCnorm"], randDistalCisPairsCombined[,"HiCnorm"]),
        captureC_raw = c(orDistalPairs[,"captureC_raw"], distalCisPairs[,"captureC_raw"], randDistalCisPairsCombined[,"captureC_raw"]),
        captureC_ObsExp = c(orDistalPairs[,"captureC_ObsExp"],  distalCisPairs[,"captureC_ObsExp"], randDistalCisPairsCombined[,"captureC_ObsExp"]),
        dist=abs(c(orDistalPairs[,"dist"], distalCisPairs[,"dist"], randDistalCisPairsCombined[,"dist"]))/10^3,
        distBin=as.factor(breaks[.bincode(abs(c(orDistalPairs[,"dist"], distalCisPairs[,"dist"], randDistalCisPairsCombined[,"dist"]))/10^3, breaks)])
    )


plotDF$HiCrawNoZero = plotDF$HiCraw
plotDF$HiCrawNoZero[plotDF$HiCraw == 0] = NA
plotDF$HiCnormNoZero = plotDF$HiCnorm
plotDF$HiCnormNoZero[plotDF$HiCnorm == 0] = NA

groupCol = c(COL_FAMILY[6], COL)
p =  ggplot(plotDF, aes(x=group, y=HiCnorm, colour = group)) + scale_y_log10() +
    geom_boxplot(lwd=1.5) + 
    scale_color_manual(values=groupCol, guide=FALSE) +
    theme_bw() + theme(text = element_text(size=20), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + 
    guides(fill=guide_legend(title="")) +
    labs(y="Normalized Hi-C counts", x="")

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
# save workspace image
#=======================================================================
save.image(WORKIMAGE_FILE)
