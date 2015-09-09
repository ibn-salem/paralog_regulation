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
# # take duplication age (dN) into account
# - for faster pipeline Run: outsource parsing of (and scanning) of TF motifs
# - Recheck sampling of gene pairs by distance.
#       - try to sample separatly for dist and enhancer
#       - plot distal pair distance with sampled pairs with log-log qqplot
#       - Maybe Sample according to distance of allCisPairs and divide close and distal afterwards

# DISCUSSION: (after discussion with Miguel on 16.06.15
# - leave out expression data for the first publication
# - Question of interest: Paralogs fitting to TAD structure of genome?
# - RECHECK carefully mouse and dog Hi-C data (and comparision to sampled genes)
# - change scale of distance box plot of orthologs
# - Synteny breaks of between mouse and human around TAD boundaries.
# X Use >= 1MB as distal pair cut-off
# - linear distance correlation for orthologs of sampled genes (should be less correlation?)
# X Put numbers to Hi-C boxplot in ortholog analysis

# FURTHER INTERESTING STUFF:
# - GTEx Consortium RNA-seq data for expression analysis
# - use other functional gene pairs (KEGG, GO (level?), PPI)
# - Colocalization of paralog pairs in other organisms (% species with shared chrom)
# - Use mouse Hi-C data from the Rao et al. 2014 paper

# EXAMPLE:
# - Check HoxA and HoxD locus in detail, as well as, IGf2/H19 locus (Kurukut et al. PNAS 2006)
# - Use the PRC1 complex as example: 
#       - check for CBX2,4,8 on chr17 and CBX6,7 on chr22
#       - PHC1,2,3 are on diff. chrom

# ADDITIONAL ANALYSIS:
# - repeat all analysis with all pairs, only two pairs
# - CHECK: ratio of synonymous mutations used for pair choosing?
# - check robustness to 1MB distance cutoff
# - build sampled background separately for enhancer, and TAD/Hi-C 
# - Check stable TADs, check definition. Why not significant?
# - Significance test on expression correlation
# - check linear distance conservation with correlation p-value
# - Compute fraction of one-to-one orthologs within the same TAD from 
#   only those human paralogs that are in the same TAD
# - include size of mouse and dog TADs in the size boxplot of all TADs

# TO FINALIZE PIPELINE FOR SUBMISSION:
# - redesign code to run on MOGON server
# - make all gene pair data frames as.character()
# - use orghologMouse instead of orthologAll data set
# - use seed() command for reproducible randomizations

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

# load some custom functions
source("R/functions.plot.R")
source("R/functions.regMap.R")
source("R/functions.genePairs.R")
source("R/functions.genePairs.randomization.R")
source("R/functions.genePairs.paralog_analysis.R")

# load ensemble data sets of genes and paralog pairs
source("R/data.ensembl.R")
source("R/data.expression.R")  # load expression data from EBI expression atlas
source("R/data.captureHiC.R")  # load capture Hi-C data between promoters from Mifsud et al. 2015

#-----------------------------------------------------------------------
# Define some parameters
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
    Dixon_hESC="data/Dixon2012/hESC.hg19.bed",
    Dixon_IMR90="data/Dixon2012/IMR90.hg19.bed"
)

# Mouse and dog TADs from Rudan et al. 2015
RudanFile = "data/Rudan2015/mmc2.xlsx"

# Hi-C data
CELL="IMR90"
HIC_RESOLUTION=50*10^3 # 50kb
HIC_DATA_DIR="data/Rao2014"

# two colors for paralog and non-paralog genes
COL=brewer.pal(9, "Set1")[c(1,9)]   # for paralog vs. sampled genes
COL_RAND=c(COL[1], brewer.pal(8, "Dark2")[8])   # random genes
COL2=brewer.pal(9, "Set1")[c(1,2)]  # for paralog vs. non-paralog
COL_DOMAIN=brewer.pal(9, "Set3")
COL_FAMILY=brewer.pal(8, "Dark2")
COL_EH = brewer.pal(9, "Set1")[5]
COL_TAD = brewer.pal(4, "Set1")[3]
COL_ORTHO = c(brewer.pal(12, "Paired")[5], brewer.pal(8, "Pastel2")[8])

# use local data or downlaod data from ensemble
USE_LOCAL = TRUE
USE_LOCAL_HIC_CONTACTS = FALSE
N_RAND=100              # number of random permutations of whole data sets
DENSITY_BW_ADJUST=0.1   # parameter to adjust bandwidth in density estimation for sampling

MAX_DIST=10^6
#~ DISTAL_MAX_DIST=10^7
#~ DISTAL_MIN_DIST=10^6
DISTAL_MAX_DIST=10^9
DISTAL_MIN_DIST=10^6

DIST_TH=c(10^6, 10^5)

PROMOTER_UPSTREAM=800
PROMOTER_DOWNSTREAM=200

VERSION="v08"
outDataPrefix = "results/paralog_regulation/EnsemblGRCh37_paralog_genes"
outPrefix = paste0("results/paralog_regulation/", VERSION, "_maxDist_", MAX_DIST, "_nrand_", N_RAND, "/EnsemblGRCh37_paralog_genes")
# create directory, if not exist
dir.create(dirname(outPrefix), recursive=TRUE, showWarnings = FALSE)

WORKIMAGE_FILE = paste0(outPrefix, "workspace.Rdata")
#load(WORKIMAGE_FILE)
#source("R/functions.genePairs.paralog_analysis.R")

#-----------------------------------------------------------------------
# make GenomicRange objects for TSS of Genes:
#-----------------------------------------------------------------------
tssGR = getTssGRfromENSEMBLGenes(genes, seqInfoRealChrom, colNames=c("hgnc_symbol"))
tssGR$gene_size = genes[names(tssGR), "end_position"] - genes[names(tssGR), "start_position"]

promoterGR = promoters(tssGR, upstream=PROMOTER_UPSTREAM, downstream=PROMOTER_DOWNSTREAM)
export(promoterGR, con=paste0(outDataPrefix, ".promoters.bed"), format="BED")

## TODO:
## Here the matrix-scan analysis should be executed. 
## See matrix-scan_pipeline.sh script
##

# get motif match table:
matrixScanFile = paste0(outDataPrefix, ".promoters.bed.names.fa.jaspar.matrix-scan.uniquePos")
motifTable = parseMatrixScan(matrixScanFile, names(promoterGR))
motifTableBinary = motifTable
motifTableBinary[motifTable>0] = 1

# import motif hit table:

write.table(motifTable, file=paste0(matrixScanFile, ".counts.tab"), sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)

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
RaoTADs = lapply(RaoDomainFiles, parseDomainsRao, disjoin=FALSE, seqinfo=seqInfo)
allTADs = c(
    RaoTADs,
    lapply(DixonDomainFiles, import.bed, seqinfo=seqInfo),
    "stable_TADs"=getConservedTADs(RaoTADs, maxgap=10^4, n=3)
    )
    
#~ allBoundaries = lapply(allTADs, getNonOverlappingBoundaries, min.gapwidth=1)

# write all TADs to BED files:
for (tadName in names(allTADs)){
    TAD = allTADs[[tadName]]
    export(sort(TAD), paste0(outPrefix, ".TAD_data.", tadName, ".bed"))
}

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

require(VennDiagram)

# plot overlap of ENSG IDs in expression data and gene set:
for (expName in names(expDFlist)) {

    expDF = expDFlist[[expName]]

    setList = list(names(tssGR), row.names(expDF))
    names(setList) = c("ENSEMBL", paste(expName, "expression"))
    vennPlotFile = paste0(outPrefix, ".", expName, ".geneIDs_ensembl_vs_expression.vennDiagram.pdf")
    # plot venn diagram
    v = venn.diagram(setList, filename=vennPlotFile, fill=COL2, 
            main=paste("Ensembl Gene ID overlap in", expName), main.cex=1.3, main.fontfamily=3,
            cex=1.3, cat.cex=1.3, fontfamily=3, cat.fontfamily=1.3, cat.default.pos="text", 
        )
}

#=======================================================================
# Annotate and filter paralog gene pairs
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


# add pairwise correlations of gene expression over all tissues
for (expName in names(expDFlist)) {
    message(paste("INFO: annotate pairs with expression correlation form:", expName))
    expDF = expDFlist[[expName]]
    
    allCisPairs = addCor(allCisPairs, expDF, colName=paste0(expName, "_expCor"))
    
    # add maximal information coefficient (MIC)
    allCisPairs = addMIC(allCisPairs, expDF, colName=paste0(expName, "_MIC"))
 
}

# add correlation over motif matches
allCisPairs = addCor(allCisPairs, motifTable, colName="motif_cor")
allCisPairs = addMIC(allCisPairs, motifTable, colName="motif_MIC")
allCisPairs = addCor(allCisPairs, motifTableBinary, colName="motif_binary_cor")
allCisPairs = addMIC(allCisPairs, motifTableBinary, colName="motif_binary_MIC")

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


#=======================================================================
# Sample random control/background data sets
#=======================================================================

#-----------------------------------------------------------------------
# Sample pairs with equal probability from all genes
#-----------------------------------------------------------------------
randPairs = replicate(10*N_RAND, getRandomPairs(nrow(paralogPairsUniqGuniqP), names(tssGR)), simplify=FALSE)

randPairs = lapply(randPairs, addSameStrand, tssGR)

# filter for random pairs in Cis
randPairsInCis = lapply(randPairs, getCisPairs, tssGR)

randPairsInCis = lapply(randPairsInCis, addOrthologAnnotation, orthologsSpeciesList[["mmusculus"]], "mmusculus", tssGRmouse, speciesTADs[["mmusculus"]], speciesHiC[["mmusculus"]][[1]], speciesHiC[["mmusculus"]][[2]])
randPairsInCis = lapply(randPairsInCis, addOrthologAnnotation, orthologsSpeciesList[["cfamiliaris"]], "cfamiliaris", tssGRdog, speciesTADs[["cfamiliaris"]], speciesHiC[["cfamiliaris"]][[1]], speciesHiC[["cfamiliaris"]][[2]])

#-----------------------------------------------------------------------
# Sample cis pairs according to enhancer number in paralogs and linear distance in paralog gene pairs
#-----------------------------------------------------------------------

# get all possible gene pairs within MAX_DIST bp
allGenePairs = getAllGenePairs(tssGR, maxDist=MAX_DIST)

# get sample weights according to enahncer number and distance
cisWeights = getSampleWeightsByDistAndEnhancers(allGenePairs, tssGR, cisPairs, adjust=DENSITY_BW_ADJUST)

# sample according to enahncer number and distance
randCisPairs = replicate(N_RAND, 
    sampleFromAllPairsByWeight(n=nrow(cisPairs), hitDF=allGenePairs, tssGR, weight=cisWeights)
    , simplify=FALSE)

    
#-----------------------------------------------------------------------
# Sample distal cis pairs the same way
#-----------------------------------------------------------------------

# Now sample from all possible gene pairs within DISTAL_MIN_DIST - DISTAL_MAX_DIST bp
allDistalGenePairs = getAllGenePairs(tssGR, maxDist=DISTAL_MAX_DIST, minDist=DISTAL_MIN_DIST)


# get observed probabilty density of distances and linked enhancer
#~ distalWeight = getSampleWeightsByDistAndEnhancers(allDistalGenePairs, tssGR, distalCisPairs, adjust=DENSITY_BW_ADJUST)

distalWeight = getSampleWeightsByDist(hitDF=allDistalGenePairs, sourcePairs=distalCisPairs, adjust=DENSITY_BW_ADJUST)

# sample according to enahncer number and distance
randDistalCisPairs = replicate(N_RAND, 
    sampleFromAllPairsByWeight(n=nrow(distalCisPairs), hitDF=allDistalGenePairs, tssGR, weight=distalWeight)
    , simplify=FALSE)
    

#-----------------------------------------------------------------------
# annotate sampled gene pairs
#-----------------------------------------------------------------------
randCisPairs = lapply(randCisPairs, addHGNC, tssGR)
randCisPairs = lapply(randCisPairs, addSameStrand, tssGR)
randCisPairs = lapply(randCisPairs, addCommonEnhancer, gene2ehID)

for (expName in names(expDFlist)) {
    expDF = expDFlist[[expName]]
    
    message(paste("INFO: Annotate sampled pairs with expression form:", expName))
    randCisPairs = lapply(randCisPairs, addCor, expDF, colName=paste0(expName, "_expCor"))
    randCisPairs = lapply(randCisPairs, addMIC, expDF, colName=paste0(expName, "_MIC"))

}

# add correlation over motif matches
randCisPairs = lapply(randCisPairs, addCor, motifTable, colName="motif_cor")
randCisPairs = lapply(randCisPairs, addMIC, motifTable, colName="motif_MIC")
randCisPairs = lapply(randCisPairs, addCor, motifTableBinary, colName="motif_binary_cor")
randCisPairs = lapply(randCisPairs, addMIC, motifTableBinary, colName="motif_binary_MIC")

randCisPairs = lapply(randCisPairs, addOrthologAnnotation, orthologsSpeciesList[["mmusculus"]], "mmusculus", tssGRmouse, speciesTADs[["mmusculus"]], speciesHiC[["mmusculus"]][[1]], speciesHiC[["mmusculus"]][[2]])
randCisPairs = lapply(randCisPairs, addOrthologAnnotation, orthologsSpeciesList[["cfamiliaris"]], "cfamiliaris", tssGRdog, speciesTADs[["cfamiliaris"]], speciesHiC[["cfamiliaris"]][[1]], speciesHiC[["cfamiliaris"]][[2]])

# add promoter-promoter contacts from caputre Hi-C
for (i in 1:N_RAND){
    randCisPairs[[i]][,"captureC_raw"] <- getPairwiseMatrixScore(randCisPairs[[i]], captureHiC[["raw"]], tssGR, replaceZeroByNA=TRUE)
    randCisPairs[[i]][,"captureC_ObsExp"] <- getPairwiseMatrixScore(randCisPairs[[i]], captureHiC[["obsExp"]], tssGR, replaceZeroByNA=TRUE)
    
    # same for sampled distal pairs
    randDistalCisPairs[[i]][,"captureC_raw"] <- getPairwiseMatrixScore(randDistalCisPairs[[i]], captureHiC[["raw"]], tssGR, replaceZeroByNA=TRUE)
    randDistalCisPairs[[i]][,"captureC_ObsExp"] <- getPairwiseMatrixScore(randDistalCisPairs[[i]], captureHiC[["obsExp"]], tssGR, replaceZeroByNA=TRUE)
}

randDistalCisPairs = lapply(randDistalCisPairs, addHGNC, tssGR)
randDistalCisPairs = lapply(randDistalCisPairs, addSameStrand, tssGR)

randDistalCisPairs = lapply(randDistalCisPairs, addOrthologAnnotation, orthologsSpeciesList[["mmusculus"]], "mmusculus", tssGRmouse, speciesTADs[["mmusculus"]], speciesHiC[["mmusculus"]][[1]], speciesHiC[["mmusculus"]][[2]])
randDistalCisPairs = lapply(randDistalCisPairs, addOrthologAnnotation, orthologsSpeciesList[["cfamiliaris"]], "cfamiliaris", tssGRdog, speciesTADs[["cfamiliaris"]], speciesHiC[["cfamiliaris"]][[1]], speciesHiC[["cfamiliaris"]][[2]])


#-----------------------------------------------------------------------
# save a work image after sampling and annotation.
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
    
    writeHiCinteractions(HiClist[1], paste0(outDataPrefix, ".Hi-C.", CELL, ".", HIC_RESOLUTION, ".map.track.txt"))

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
# 5.) Run analysis
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


pdf(paste0(outPrefix, ".sampling_cis_and_distal.pdf"))
    par(cex=1, lwd=1.5, mfrow=c(3,3))

    paraLinkedEnhancer = tssGR[c(cisPairs[,1], cisPairs[,2])]$linked_enhancer
    randLinkedEnhancer = unlist(lapply(randCisPairs, function(gP) tssGR[c(gP[,1], gP[,2])]$linked_enhancer))
    
    qqplot(paraLinkedEnhancer, randLinkedEnhancer, xlab="Enhancers in paralog genes", ylab="Enhancers in sampled genes", main="QQ-Plot of linked enhancers\n in close cis pairs")
    abline(0,1, col="red")

    hist(paraLinkedEnhancer[paraLinkedEnhancer<=50], 50, col=COL[1],
    main="Close cis paralogs", xlab="Enhancers")
    hist(randLinkedEnhancer[randLinkedEnhancer<=50], 50, col=COL[2],
    main="Close cis sampled genes", xlab="Enhancers")    

    # verify same distribution of distances
    paraDist = cisPairs$dist / 10^3
    randDist = unlist(lapply(randCisPairs, function(d){d[,"dist"]})) / 10^3
    
    qqplot(abs(paraDist), abs(randDist), xlab="Disance between paralogs", ylab="Distance between sampled gene pairs", main="QQ-Plot of distances\n in in close cis pairs")
    abline(0,1, col="red")

    hist(abs(paraDist), 50, col=COL[1],
    main="Close cis paralogs", xlab="Distance (kb)")
    hist(abs(randDist), 50, col=COL[2],
    main="Close cis sampled genes", xlab="Distance (kb)")    

    # same for distal pairs
    paraDistalDist = distalCisPairs$dist / 10^3
    randDistalDist = unlist(lapply(randDistalCisPairs, function(d){d[,"dist"]})) / 10^3
    
    qqplot(abs(paraDistalDist), abs(randDistalDist), log="xy",
    xlab="Disance between paralogs", ylab="Distance between sampled gene pairs", main="QQ-Plot of distances\n in distal cis pairs")
    abline(0,1, col="red")

    hist(abs(paraDistalDist), 50, col=COL[1],
    main="Distal cis paralogs", xlab="Distance (kb)")
    hist(abs(randDistalDist), 50, col=COL[2],
    main="Distal cis sampled genes", xlab="Distance (kb)")    

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
randChromPairMatrix = lapply(randPairs, interChromPairMatrix, tssGR)

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

pdf(paste0(outPrefix, ".random_genes_distance.hist.pdf"))
    par(cex=1.5, lwd=2, mfrow=c(2,1))
    hist(abs(paraDist)/10^3, 50, col=COL_RAND[1],
    main="Distance between paralog genes", xlab="Distance (kb)")
    hist(abs(randDist)/10^3, 50, col=COL_RAND[2],
    main="Distance between random genes", xlab="Distance (kb)")    
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
pdf(paste0(outPrefix, ".hESC_TAD_size.hist.pdf"))
    par(cex=1.5, lwd=2)
    hist(width(allTADs[["Dixon_hESC"]])/10^3, col=COL_DOMAIN[1],
        main="Size distribution of hESC TADs\n (Dixon et al. 2012)",
        xlab="Domain size (kb)")
dev.off()
pdf(paste0(outPrefix, ".Rao_GM12878_TAD_size.hist.pdf"))
    par(cex=1.5, lwd=2)
    hist(width(allTADs[["Rao_GM12878"]])/10^3, col=COL_DOMAIN[2],
        main="Size distribution of GM12878 TADs\n (Rao et al. 2014)",
        xlab="TAD size (kb)")
dev.off()

sizeList = lapply(allTADs, function(gr) width(gr)/10^3)

pdf(paste0(outPrefix, ".Compare_domain_size.boxplot.pdf"))
    par(cex=1.3, lwd=1.5)
    my.boxplot(sizeList, names=gsub("_", " ", names(sizeList)),
        main="Size distribution of TADs",
        ylab="TAD size [kb]", col=COL_DOMAIN)
dev.off()

#-----------------------------------------------------------------------
# Co-occurance in same interaction domain
#-----------------------------------------------------------------------

# make GRanges objects for cis paralog pairs and random paris on same chromosome
cisParaGR = getPairAsGR(cisPairs, tssGR)
cisRandGR = lapply(randCisPairs, getPairAsGR, tssGR)


#~ # co-occurance within the same domain
#~ cisParaGR = addWithinSubject(cisParaGR, hESC_TADs, "hESC_TADs")
#~ cisRandGR = lapply(cisRandGR, addWithinSubject, hESC_TADs, "hESC_TADs")
#~ 
#~ # same with Hi-C domains from Rao et al 2014
#~ cisParaGR = addWithinSubject(cisParaGR, raoDomDisjoin, "Rao_Domains")
#~ cisRandGR = lapply(cisRandGR, addWithinSubject, raoDomDisjoin, "Rao_Domains")

for(tadName in names(allTADs)){

    TAD = allTADs[[tadName]]

    # co-occurance within the same domain
    cisParaGR = addWithinSubject(cisParaGR, TAD, tadName)
    cisRandGR = lapply(cisRandGR, addWithinSubject, TAD, tadName)
    
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

        add_pval_two_pairs(bp, heightMat, pvalues, offset=.15*max(heightMat), digits=1)

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

# for capture Hi-C
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
paraExpMICSummaryList=list()

# iterate over all expression data sets
for (expName in names(expDFlist)){

    expDF = expDFlist[[expName]]

    # plot distribution of correlation coefficient
    paraExpCor = cisPairs[,paste0(expName, "_expCor")]
    sampledExpCor = randCisPairsAll[,paste0(expName, "_expCor")]

    paraExpMIC = cisPairs[,paste0(expName, "_MIC")]
    sampledExpMIC = randCisPairsAll[,paste0(expName, "_MIC")]
    
    pdf(paste0(outPrefix,".sampled_genes.Expression_correlation.", expName, ".hist.pdf"))
        par(cex=1.5, lwd=2, mfrow=c(2,1))
        breaks = seq(-1,1,length.out=51)
        xlab = paste0("Gene expression correlation over n=", ncol(expDF), " tissues [Pearson's R]")
        hist(paraExpCor, breaks, col=COL[1],
        main=paste0("Paralog genes co-expression (",expName,")"), xlab=xlab)
        hist(sampledExpCor, breaks, col=COL[2],
        main=paste0("Sampled genes co-expression (",expName,")"), xlab=xlab)    
    dev.off()

    pdf(paste0(outPrefix,".sampled_genes.Expression_MIC.", expName, ".hist.pdf"))
        par(cex=1.5, lwd=2, mfrow=c(2,1))
        breaks = seq(0,1,length.out=51)
        xlab = paste0("Maximal information coefficient (MIC) over n=", ncol(expDF), " tissues")
        hist(paraExpMIC, breaks, col=COL[1],
        main=paste0("Paralog genes (",expName,")"), xlab=xlab)
        hist(sampledExpMIC, breaks, col=COL[2],
        main=paste0("Sampled genes (",expName,")"), xlab=xlab)    
    dev.off()
    
    # add vectors to summary plot list
    paraExpCorSummaryList = c(paraExpCorSummaryList, list(paraExpCor, sampledExpCor))
    paraExpMICSummaryList = c(paraExpMICSummaryList, list(paraExpMIC, sampledExpMIC))
    
    # plot distribution of mutual information (MI)
    #paraExpMI = cisPairs[,paste0(expName, "_expMI")]
    #sampledExpMI = randCisPairsAll[,paste0(expName, "_expMI")]
    
#~     pdf(paste0(outPrefix,".sampled_genes.Expression_MI.", expName, ".hist.pdf"))
#~         par(cex=1.5, lwd=2, mfrow=c(2,1))
#~         breaks = seq(0,1,length.out=51)
#~         xlab = paste0("Mutual information of gene expression over n=", ncol(expDF), " tissues [bits]")
#~         hist(paraExpMI, breaks, col=COL[1],
#~         main=paste0("Paralog genes co-expression (",expName,")"), xlab=xlab)
#~         hist(sampledExpMI, breaks, col=COL[2],
#~         main=paste0("Sampled genes co-expression (",expName,")"), xlab=xlab)    
#~     dev.off()
#~     
#~     # add vectors to summary plot list
#~     paraExpMISummaryList = c(paraExpMISummaryList, list(paraExpMI, sampledExpMI))

    #-------------------------------------------------------------------
    # Filter for exclusively expressed genes with an non-coexistence pattern
    #-------------------------------------------------------------------
    MIN_MICR2 = .25
    MAX_R = 0
    
    r = cisPairs[,paste0(expName, "_expCor")]
    mic = cisPairs[,paste0(expName, "_MIC")]
    micr2 =  mic - r**2
    
    plotPairs = cisPairs[which(micr2 >= MIN_MICR2 & r <= MAX_R), ]

    plotAllExp(plotPairs, expDF, tssGR,
        paste0(outPrefix,".expression.noncoexistence.", expName, ".minMICR2_", MIN_MICR2, "_maxR_", MAX_R, ".hist.pdf"))


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
    
        p = ggplot(cisPairs, aes(x=inTADgroup, y=paraExpCor)) + geom_boxplot(fill=NA, aes(colour = inTADgroup), lwd=1.5) + geom_jitter(alpha=.3, size=2, aes(colour = inTADgroup)) + scale_colour_manual(values= colorRampPalette(c("gray", COL_TAD))(2), guide=FALSE) + theme_bw() + theme(text = element_text(size=20), axis.text.x = element_text(angle = 45, hjust = 1)) + labs(y="Pearson correlation coefficient", x = "", title=paste0("p = ",signif(ws.test$p.value,2)))
        # save ggplot as pdf
        ggsave(p, file=paste0(outPrefix,".expCor_", expName, ".by_sahred_TAD_", tadName, ".boxplot.pdf"), width=3.5, height=7)
    }

    # plot for all cell types: 
    expPlotDF = data.frame(
        expCor = rep(paraExpCor, length(allTADs)),
        expMIC = rep(paraExpMIC, length(allTADs)),
        withinTAD = unlist(lapply(names(allTADs), function(tadName) factor(mcols(cisParaGR)[, tadName], levels=c(TRUE, FALSE), labels=c("In same TAD", "Not in same TAD")))),
        cells = rep(factor(gsub("_", " ", names(allTADs)), gsub("_", " ", names(allTADs))) , each=length(paraExpCor))
    )
    p = ggplot(expPlotDF, aes(x=withinTAD, y=expCor, colour=withinTAD)) + geom_violin(aes(colour = withinTAD), lwd=1, alpha=.25) + geom_boxplot(aes(color=withinTAD), fill=NA, width=.25, lwd=1)  + 
    facet_wrap(~ cells) + scale_colour_manual(values= colorRampPalette(c(COL_TAD, "gray"))(2), guide=FALSE) + theme_bw() + theme(text = element_text(size=20), axis.text.x = element_text(angle = 45, hjust = 1)) + labs(y="Pearson correlation coefficient", x = "") +  theme(legend.position="bottom")

    # save ggplot as pdf
    ggsave(p, file=paste0(outPrefix,".expCor_", expName, ".by_sahred_TAD.allTADs.boxplot.pdf"), w=7, h=7)
    
    p = ggplot(expPlotDF, aes(x=withinTAD, y=expMIC, colour=withinTAD)) + geom_violin(aes(colour = withinTAD), lwd=1, alpha=.25) + geom_boxplot(aes(color=withinTAD), fill=NA, width=.25, lwd=1)  + 
    facet_wrap(~ cells) + scale_colour_manual(values= colorRampPalette(c(COL_TAD, "gray"))(2), guide=FALSE) + theme_bw() + theme(text = element_text(size=20), axis.text.x = element_text(angle = 45, hjust = 1)) + labs(y="Maximal information coefficient (MIC)", x = "") +  theme(legend.position="bottom")

    # save ggplot as pdf
    ggsave(p, file=paste0(outPrefix,".expMIC_", expName, ".by_sahred_TAD.allTADs.boxplot.pdf"), w=7, h=7)
    
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

# same for mutual information
pdf(paste0(outPrefix,".sampled_genes.Expression_MIC.all_data.boxplot.pdf"))
    par(cex=1.3, lwd=2, mar=c(8, 4, 1, 2))
    boxplot(paraExpMICSummaryList, border=COL, ylab="Maximal information coefficient (MIC)", names=NA, xaxt="n")
    
    # add x-axis labels
    labPos = colMeans(matrix(seq(2*nExp), nrow=2))
    axis(1, at=labPos, labels=FALSE)
    text(x=labPos, y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]),labels=names(expDFlist), srt=45, adj=1, xpd=TRUE)
    
dev.off()

# make data.frame with additional columns indicating the groups
paraExpCorSummaryDF <- data.frame(
                expCor = unlist(paraExpCorSummaryList), 
                expMIC = unlist(paraExpMICSummaryList), 
                expMIC_R2 = unlist(paraExpMICSummaryList) - unlist(paraExpCorSummaryList)**2, 
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

# same for MIC
pdf(paste0(outPrefix,".sampled_genes.Expression_MIC.all_data.violinplot.pdf"))
    ggplot(paraExpCorSummaryDF, na.rm=TRUE, aes(x=genepairs, y = expMIC, fill=genepairs)) + geom_violin(adjust = .2) + theme_bw() + labs(y="MIC", x = "", title= "Co-expression over tissues") + scale_x_discrete(labels="")  + facet_grid(.~dataset) + theme(legend.position = "bottom") + scale_fill_manual(values=COL) + geom_boxplot(aes(fill=NULL), width=.2, colour="black", outlier.shape = NA)+guides(fill=guide_legend(title=""))
dev.off()

# same for MIC-R2
pdf(paste0(outPrefix,".sampled_genes.Expression_MIC-R2.all_data.violinplot.pdf"))
    ggplot(paraExpCorSummaryDF, na.rm=TRUE, aes(x=genepairs, y = expMIC_R2, fill=genepairs)) + geom_violin(adjust = .2) + theme_bw() + labs(y="MIC - R^2", x = "", title= "Co-expression over tissues") + scale_x_discrete(labels="")  + facet_grid(.~dataset) + theme(legend.position = "bottom") + scale_fill_manual(values=COL) + geom_boxplot(aes(fill=NULL), width=.2, colour="black", outlier.shape = NA)+guides(fill=guide_legend(title=""))
dev.off()

#-----------------------------------------------------------------------
# TF motif analysis in promoters
#-----------------------------------------------------------------------
pdf(paste0(outPrefix,".motif_in_promoter.example.heatmap.pdf"))

    simple.heatmap(as.matrix(motifTable[39:80, 1:50]), labRow=tssGR[row.names(motifTable)[39:80]]$hgnc_symbol, xlab="TF motifs", ylab="Promoters", margins=c(2,2))

dev.off()

# plot distribution of correlation coefficient
paraMotifCor = cisPairs[,"motif_cor"]
sampledMotifCor = randCisPairsAll[,"motif_cor"]

paraMotifMIC = cisPairs[,"motif_MIC"]
sampledMotifMIC = randCisPairsAll[,"motif_MIC"]

pdf(paste0(outPrefix,".sampled_genes.motifs_correlation.density.pdf"))
    par(cex=1.5, lwd=3)
    yMax = 1.3 * max(density(paraMotifCor)$y)
    xlab = paste0("TF motifs correlation over n=", ncol(motifTable), " TF motifs [Pearson's R]")
    plot(density(paraMotifCor), xlim=c(-1,1), ylim=c(0,yMax), col=COL[1],
    main="TF Motif Correlation", xlab=xlab)
    lines(density(sampledMotifCor), col=COL[2])
    legend("topleft", c("Paralog genes", "Sampled genes"), col=COL, lty=1)    
dev.off()

pdf(paste0(outPrefix,".sampled_genes.motifs_correlation.hist.pdf"))
    par(cex=1.5, lwd=2, mfrow=c(2,1))
    breaks = seq(-1,1,length.out=51)
    xlab = paste0("TF motifs correlation over n=", ncol(motifTable), " TF motifs [Pearson's R]")
    hist(paraMotifCor, breaks, col=COL[1],
    main="Paralog genes", xlab=xlab)
    hist(sampledMotifCor, breaks, col=COL[2],
    main="Sampled genes", xlab=xlab)    
dev.off()

pdf(paste0(outPrefix,".sampled_genes.motifs_MIC.hist.pdf"))
    par(cex=1.5, lwd=2, mfrow=c(2,1))
    breaks = seq(0,1,length.out=51)
    xlab = paste0("TF Motif associations over n=", ncol(motifTable), " TF motifs [MIC]")
    hist(paraMotifMIC, breaks, col=COL[1],
    main="Paralog genes", xlab=xlab)
    hist(sampledMotifMIC, breaks, col=COL[2],
    main="Sampled genes", xlab=xlab)    
dev.off()

# plot distribution of correlation coefficient 
paraMotifCor = cisPairs[,"motif_binary_cor"]
sampledMotifCor = randCisPairsAll[,"motif_binary_cor"]

paraMotifMIC = cisPairs[,"motif_binary_MIC"]
sampledMotifMIC = randCisPairsAll[,"motif_binary_MIC"]

pdf(paste0(outPrefix,".sampled_genes.motifs_binary_correlation.density.pdf"))
    par(cex=1.5, lwd=3)
    yMax = 1.3 * max(density(paraMotifCor)$y)
    xlab = paste0("TF motifs correlation over n=", ncol(motifTable), " TF motifs [Pearson's R]")
    plot(density(paraMotifCor), xlim=c(-1,1), ylim=c(0,yMax), col=COL[1],
    main="TF Motif Correlation", xlab=xlab)
    lines(density(sampledMotifCor), col=COL[2])
    legend("topleft", c("Paralog genes", "Sampled genes"), col=COL, lty=1)    
dev.off()

pdf(paste0(outPrefix,".sampled_genes.motifs_binary_MIC.density.pdf"))
    par(cex=1.5, lwd=3)
    yMax = 1.3 * max(density(paraMotifMIC)$y)
    xlab = paste0("TF motifs correlation over n=", ncol(motifTableBinary), " TF motifs [MIC]")
    plot(density(paraMotifMIC), xlim=c(0,1), ylim=c(0,yMax), col=COL[1],
    main="TF Motif Correlation", xlab=xlab)
    lines(density(sampledMotifMIC), col=COL[2])
    legend("topleft", c("Paralog genes", "Sampled genes"), col=COL, lty=1)    
dev.off()

pdf(paste0(outPrefix,".sampled_genes.motifs_binary_correlation.hist.pdf"))
    par(cex=1.5, lwd=2, mfrow=c(2,1))
    breaks = seq(-1,1,length.out=51)
    xlab = paste0("TF motifs correlation over n=", ncol(motifTable), " TF motifs [Pearson's R]")
    hist(paraMotifCor, breaks, col=COL[1],
    main="Paralog genes", xlab=xlab)
    hist(sampledMotifCor, breaks, col=COL[2],
    main="Sampled genes", xlab=xlab)    
dev.off()

pdf(paste0(outPrefix,".sampled_genes.motifs_binary_MIC.hist.pdf"))
    par(cex=1.5, lwd=2, mfrow=c(2,1))
    breaks = seq(0,1,length.out=51)
    xlab = paste0("TF Motif associations over n=", ncol(motifTable), " TF motifs [MIC]")
    hist(paraMotifMIC, breaks, col=COL[1],
    main="Paralog genes", xlab=xlab)
    hist(sampledMotifMIC, breaks, col=COL[2],
    main="Sampled genes", xlab=xlab)    
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
# Compare get one2one orthologs of human paralogs in mouse and dog
#=======================================================================
#allCisPairs
#randPairsInCis

orgStr2Name = c(mmusculus="mouse", cfamiliaris="dog")

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
    
#~ 
#~     # correlate linear distances of random genes in human with mouse orthologs
#~     pdf(paste0(outPrefix, ".paralogPairs_orthologs_rand_dist_", orgName, ".dotplot.pdf"))
#~ 
#~         humanDist = abs(unlist(
#~             lapply(randPairsInCis, function(d){
#~                 subSet = (abs(d[, orthoDistStr]) <= MAX_DIST) & (abs(d[,"dist"])<= MAX_DIST)
#~                 d[, orthoDistStr]
#~                 })
#~             )) / 10^3
#~ 
#~         maxDistSubset = abs(allCisPairs[, "dist"]) <= MAX_DIST & abs(allCisPairs[, paste0(orgStr, "_dist")]) <= MAX_DIST 
#~ 
#~         humanDist = abs(allCisPairs[maxDistSubset, "dist"]) #/10^3
#~         orthoDist = abs(allCisPairs[maxDistSubset, paste0(orgStr, "_dist")]) #/10^3
#~ 
#~         r = cor(humanDist[!is.na(orthoDist)], orthoDist[!is.na(orthoDist)])
#~         p = cor.test(humanDist[!is.na(orthoDist)], orthoDist[!is.na(orthoDist)])$p.value
#~         
#~         par(lwd=2, cex=1.5)
#~         plot(log10(humanDist),  log10(orthoDist),
#~             xlim=c(3, log10(MAX_DIST)), ylim=c(3, log10(MAX_DIST)),
#~             main=paste("Distance in human and", orgName),
#~             xlab="log_10 distance in human", ylab=paste("log_10 distance in", orgName) ) #, col=rgb(0,0,0,.5)
#~         legend("topleft", paste("R =", signif(r, 3)))
#~          
#~     dev.off()

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
#~     for (tadName in names(allTADs)){
#~ 
#~         paraInTAD <- mcols(cisParaGR)[, tadName]
#~         randInTAD = lapply(cisRandGR, function(gr) mcols(gr)[, tadName])
#~ 
#~         humanAndOrthoInTAD = cisPairs[cisPairs[,paste0(orgStr, "_one2one")] & paraInTAD, paste0(orgStr, "_TAD")]
#~         humanAndRandInTAD = lapply(1:length(randCisPairs), function(i){ 
#~             gP = randCisPairs[[i]]
#~             gP[gP[,paste0(orgStr, "_one2one")] & randInTAD[[i]], paste0(orgStr, "_TAD")]
#~             })
#~ 
#~         # create contingency table and run Fisher test
#~         contab = rbind(
#~             para=table(humanAndOrthoInTAD),  
#~             rand=table(unlist(humanAndRandInTAD))
#~         )
#~         fs.test = fisher.test(contab)
#~     }
        
    
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
# Mouse paralog analysis
#=======================================================================
runBasicParalogAnalysis(paste0(outPrefix, ".Mouse"), paralogPairsMouse, "mmusculus_paralog_dn", tssGRmouse, speciesTADs[["mmusculus"]], tissueName="Mouse", HiClist=speciesHiC[["mmusculus"]][[1]], HiClistNorm=speciesHiC[["mmusculus"]][[2]])
message("Finish Mouse paralog analysis!")

#=======================================================================
# Dog and dog paralog analysis
#=======================================================================
runBasicParalogAnalysis(paste0(outPrefix, ".Dog"), paralogPairsDog, "cfamiliaris_paralog_dn", tssGRdog, speciesTADs[["cfamiliaris"]], tissueName="Dog", HiClist=speciesHiC[["cfamiliaris"]][[1]], HiClistNorm=speciesHiC[["cfamiliaris"]][[2]])
message("Finish Dog paralog analysis!")


#=======================================================================
# Best candidate picking...
#=======================================================================
candidates = with(cisPairs, commonEnhancer>3 & sameStrand & abs(ENCODE_cell_lines_expCor) > .6 ) & cisParaGR$stable_TADs 

table(candidates) 

#cisPairs[candidates,]

#=======================================================================
# Seq similarity of downstream paralog for pairs with only one ortholog in mouse
#=======================================================================

#mouseOrthologsInHumanALL
#~ singleMouseOrthlogs = mouseOrthologsInHumanALL[mouseOrthologsInHumanALL$hsapiens_homolog_orthology_type == "ortholog_one2many",]

#~ cisPairs[,1] = as.character(cisPairs[,1])
#~ cisPairs[,2] = as.character(cisPairs[,2])



commonOrthologs = getUniqueCommonOrthologs(cisPairs, humanToOrtho, orthoGeneCol="mmusculus_homolog_ensembl_gene")

#-----------------------------------------------------------------------
# get percent identety for both genes
#-----------------------------------------------------------------------
getPercentID <- function(g, o, humanToOrtho, orthoGeneCol= "mmusculus_homolog_ensembl_gene", simCol="mmusculus_homolog_perc_id_r1"){
    
    if( is.na(o) ){ 
        return(NA)
    }
    uniqOrtholog = g == humanToOrtho[,1] & o == humanToOrtho[,orthoGeneCol]
    
    sim = humanToOrtho[uniqOrtholog, simCol][1]

    return(sim)
}


cisPairs$g1_orthoSim = unlist(mapply(getPercentID, cisPairs[,1], commonOrthologs, MoreArgs=list(humanToOrtho=humanToOrtho, orthoGeneCol= "mmusculus_homolog_ensembl_gene", simCol="mmusculus_homolog_perc_id_r1")))
cisPairs$g2_orthoSim = unlist(mapply(getPercentID, cisPairs[,2], commonOrthologs, MoreArgs=list(humanToOrtho=humanToOrtho, orthoGeneCol= "mmusculus_homolog_ensembl_gene", simCol="mmusculus_homolog_perc_id_r1")))


cisPairs$firstUpstream = getFirstUpstream(cisPairs, tssGR)

# n1 means that the upsteam paralog is more similar to the ortholog
cisPairs$upstreamSim = ifelse(cisPairs$firstUpstream, cisPairs$g1_orthoSim ,cisPairs$g2_orthoSim) > ifelse(cisPairs$firstUpstream, cisPairs$g2_orthoSim ,cisPairs$g1_orthoSim)

cisPairs$downstreamSim = ifelse(cisPairs$firstUpstream, cisPairs$g1_orthoSim ,cisPairs$g2_orthoSim) < ifelse(cisPairs$firstUpstream, cisPairs$g2_orthoSim ,cisPairs$g1_orthoSim) 

cisPairs$eqSim = cisPairs$g1_orthoSim == cisPairs$g2_orthoSim

sapply(list(cisPairs$eqSim, cisPairs$upstreamSim, cisPairs$downstreamSim), sum, na.rm=TRUE)

addSingleOrthologs <- function(genePairs, singleOrthologs){
    
    gene2Ortho = singleOrthologs
    rownames(gene2Ortho) = singleOrthologs[,2]

    g1 = as.character(genePairs[,1])
    g2 = as.character(genePairs[,2])

    hasSingleOrtho = (g1 %in% singleOrthologs[,2]) & (g2 %in% singleOrthologs[,2])
    
    genePairs$g1_ortho_perc_id = singleOrthologs[singleOrthologs[,2]==g1]
    
}

#=======================================================================
# save workspace image
#=======================================================================
save.image(WORKIMAGE_FILE)
