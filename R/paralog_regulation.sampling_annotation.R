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


require(biomaRt)        # to retrieve human paralogs from Ensembl
require(stringr)        # for some string functionality
require(GenomicRanges)  # for genomic intervals and overlap calculation
require(rtracklayer)    # for import.bed
require(plyr)           # count() function
require(data.table)     # for data.table object
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
source("R/data.gene_age.R")  # load duplication age for paralogs as computed by Pablo Mier Munoz

#=======================================================================
# 1.) Parse data
#=======================================================================
if (!LOAD_INPUT_DATA){
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
    # Parse TADs from Rao et al 2014 and Dixon et al 2012
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
        stableTADs,
        DixonTADs
    )
    
    message("INFO: Finshed parsing of TADs.")
    
    # write all TADs to BED files:
    for (tadName in names(allTADs)){
        TAD = allTADs[[tadName]]
        export(sort(TAD), paste0(outPrefix, ".TAD_data.", tadName, ".bed"))
    }
    
    #-----------------------------------------------------------------------
    # Parse sub compartments in GM12878 Rao et al 2014
    #-----------------------------------------------------------------------
    # parse subcompartments by restricting ranges to the  chromosome sizes
    subCompGR <- IRanges::trim(import.bed(RaoSubcompartmentFile, seqinfo=seqInfo))
    # add name column for subcompartment and compartment type
    subCompGR$subcomp <- ifelse(subCompGR$name != "NA", subCompGR$name, NA)
    subCompGR$comp <- substr(subCompGR$subcomp, 1,1)
    
    # merge all adjacent A and B subcompartment to get full compartments
    compA <- reduce(subCompGR[subCompGR$comp == "A" & !is.na(subCompGR$comp)])
    compA$comp <- "A"
    compB <- reduce(subCompGR[subCompGR$comp == "B" & !is.na(subCompGR$comp)])
    compB$comp <- "B"
    compGR <- sort(c(compA, compB))
    
    # annotate tssGR with compartment and sub-compartment
    #~ tssCompHit <- findOverlaps(tssGR, subCompGR)
    #~ mcols(tssGR)[queryHits(tssCompHit), "comp"] <- subCompGR[subjectHits(tssCompHit)]$comp
    #~ mcols(tssGR)[queryHits(tssCompHit), "subcomp"] <- subCompGR[subjectHits(tssCompHit)]$subcomp
    
    #-----------------------------------------------------------------------
    # Load Hi-C data from Rao et al. 2014
    #-----------------------------------------------------------------------
    if ( !USE_LOCAL) {
    
        # parse normalized Hi-C map from Rao et al. 
        HiClist = parseRaoHiC(CELL, HIC_RESOLUTION, HIC_DATA_DIR, seqInfo)

        HiClistRaw = parseRaoHiC(CELL, HIC_RESOLUTION, HIC_DATA_DIR, seqInfo, normStr=NULL)
        
        # parse expected contacts by distance for each chromosome
        expectedHiCList <- parseRaoExpected(CELL, HIC_RESOLUTION, HIC_DATA_DIR)
    
        
        # save data for faster loading next time
        save(HiClist, HiClistRaw, expectedHiCList, file=paste0(outDataPrefix, ".Hi-C.", CELL, ".", HIC_RESOLUTION, ".HiClist_expectedHiCList.RData"))
        
        # write Hi-C interaction in track.txt format for visuallization in Epi-Genome browser
    #~     writeHiCinteractions(HiClist[1], paste0(outDataPrefix, ".Hi-C.", CELL, ".", HIC_RESOLUTION, ".map.track.txt"))
    
    }else{
        load(paste0(outDataPrefix, ".Hi-C.", CELL, ".", HIC_RESOLUTION, ".HiClist_expectedHiCList.RData"))
    }
    
    #-----------------------------------------------------------------------
    # Parse Hi-C data for other species
    #-----------------------------------------------------------------------
    orgStr2Name = c(mmusculus="mouse", cfamiliaris="dog")
    
    # parse mouse and dog TADs
    speciesSeqInfo = list("mmusculus"=seqInfoMouse, "cfamiliaris"=seqInfoDog)
    
    
    speciesTADs = lapply(1:2, function(i) parseRudanTADs(RudanFile, sheet=i, disjoin=FALSE, seqinfo=speciesSeqInfo[[i]]))
    names(speciesTADs) = names(speciesSeqInfo)
    
    speciesHiC = list(
        "mmusculus" = parseRudanHiC("data/Rudan2015/GSE65126_HiC_mouse_liver_merged_50000.txt", seqInfoMouse, resolution=50*10^3), 
        "cfamiliaris"= parseRudanHiC("data/Rudan2015/GSE65126_HiC_dog_liver_merged_50000.txt", seqInfoDog, resolution=50*10^3)
    )
    
    #-----------------------------------------------------------------------
    # Parse HIPPIE 
    #-----------------------------------------------------------------------
    hippieDF <- read.table("data/HIPPIE/hippie_current.txt", header=FALSE, sep="\t")
    hippie <- data.frame(
        symbol1 = str_split_fixed(as.character(hippieDF[,1]), "_", 2)[,1],
        symbol2 = str_split_fixed(as.character(hippieDF[,3]), "_", 2)[,1],
        score = hippieDF[,5],
        stringsAsFactors=FALSE)

    hippie$ID <- getPairIDsorted(entrezToENSG[as.character(hippieDF[,2])], entrezToENSG[as.character(hippieDF[,4])])
            
    # save all loaded data as image file
    save.image(paste0(WORKIMAGE_FILE, ".loaded_data.Rdata"))
}else{
    load(paste0(WORKIMAGE_FILE, ".loaded_data.Rdata"))
}

#=======================================================================
# 2.) Filter and annotate paralog gene pairs
#=======================================================================
if (!LOAD_PAIRS) {
    
    # FILTERING ORDER:
    # paralogPairs                      w dups        n=
    # +-- paralogPairsUniqP             w/o dups      n=
    # +-- paralogPairsUniqG             w dups        n=
    #     +-- paralogPairsUniq
    #         +-- allCisPairs           
    #             +-- closePairs           
    #             +-- distalPairs           
    
    
    # the data.frame "paralogPairs" is loaded from data.ensembl.R
    
    # remove double entries of the form A-B and B-A
    paralogPairsUniqP = uniquePair(paralogPairs)
    
    # filter out overlapping gene pairs
    nonOVL <- nonOverlappingGenePairs(paralogPairsUniqP, genesGR)
    paralogPairsUniqPnonOVL <- paralogPairsUniqP[nonOVL,]
    
    # write all paralog pairs to output file:
    write.table(paralogPairsUniqPnonOVL, file=paste0(outPrefix, ".paralog_pairs.paralogPairsUniqPnonOVL.txt"),
        sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    
    #~ 
    #~ # get for each gene only one unique pair, the one with highest similarity
    #~ # this is computed by an maximum weight matching
    #~ paralogPairsWithDS = paralogPairs[!is.na(paralogPairs[,"hsapiens_paralog_ds"]),]
    #~ paralogPairsUniqG = uniquePairPerGeneBySim(paralogPairsWithDS, -1*paralogPairsWithDS[,"hsapiens_paralog_ds"])
    #~     
    #~ # get only a unique pair order (one of A-B, B-A) form the unique pairs
    #~ paralogPairsUniq = uniquePair(paralogPairsUniqG)
    #~ 
    #~ # subset of paralog pairs that are located on the same chromosome
    #~ allCisPairs = getCisPairs(paralogPairsUniq, tssGR)
    
    # get for each gene only one unique pair, the one with highest similarity
    # this is computed by an maximum weight matching
    paralogPairsUniqPnonOVLWithDS = paralogPairsUniqPnonOVL[!is.na(paralogPairsUniqPnonOVL[,"hsapiens_paralog_ds"]),]
    
    paralogPairsUniqPnonOVLUniqG = uniquePairPerGeneBySim(paralogPairsUniqPnonOVLWithDS, -1*paralogPairsUniqPnonOVLWithDS[,"hsapiens_paralog_ds"])
    
    # get only a unique pair order (one of A-B, B-A) form the unique pairs
    paralogPairsUniq = uniquePair(paralogPairsUniqPnonOVLUniqG)
    
    # subset of paralog pairs that are located on the same chromosome
    allCisPairs = getCisPairs(paralogPairsUniq, tssGR)
    
    #-----------------------------------------------------------------------
    # Annotate:
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
    
    # add number of common enhancers
    allCisPairs = addCommonEnhancer(allCisPairs, gene2ehID)
    
    # add position of enhancers:
    allCisPairs = addRelativeEnhancerPosition(allCisPairs, tssGR, gene2ehID, ehGR)
    
    
    # add expression in IMR90 cells
    # expDFlist is already loaded in the script "R/data.expression.R"
    allCisPairs = addPairExp(allCisPairs, expDFlist[["ENCODE_cell_lines"]], expCol="IMR_90", label="exp_IMR90")
    
    nExp = length(expDFlist)
    
    # add pairwise correlations of gene expression over all tissues
    for (expName in names(expDFlist)) {
        
        message(paste("INFO: annotate pairs with expression correlation form:", expName))
        expDF = expDFlist[[expName]]
        
        allCisPairs = addCor(allCisPairs, expDF, colName=paste0(expName, "_expCor"))
    }
    
    # add same TAD annotation
    
    # make GRanges objects for cis paralog pairs and random paris on same chromosome
    allCisPairsGR = getPairAsGR(allCisPairs, tssGR)
    
    for(tadName in names(allTADs)){
        message(paste("INFO: Compute overlap with TADs from:", tadName))
        # co-occurance within the same domain
        allCisPairsGR = addWithinSubject(allCisPairsGR, allTADs[[tadName]], tadName)
    }
  
    # assign annotation in GRanges object to gene pair data.frames
    allCisPairs[,names(allTADs)] <- data.frame( mcols(allCisPairsGR)[, names(allTADs)] )

    for(tadName in names(allTADs)){
        message(paste("INFO: Compute common subset of overallping TADs from:", tadName))
        allCisPairs = addSubTADmode(allCisPairs, allTADs[[tadName]], tssGR, paste0(tadName, "_subTAD"))
    }
    
    # Adds Hi-C contact frequencies to a gene pair data set
    allCisPairs = addHiCfreq(allCisPairs, tssGR, HiClistRaw, label="HiCRaw", ignoreSameBin=TRUE)
    Sys.sleep(3) # hack to fix problems with bplapply on MOGON
    allCisPairs = addHiCfreq(allCisPairs, tssGR, HiClist, label="HiC", ignoreSameBin=TRUE)
    Sys.sleep(3) # hack to fix problems with bplapply on MOGON
    allCisPairs = addHiCobsExp(allCisPairs, tssGR, expectedHiCList, HIC_RESOLUTION, HiClabel="HiC", label="HiCobsExp")
    
    
    # add promoter-promoter contacts from caputre Hi-C
    allCisPairs[,"captureC_raw"] = getPairwiseMatrixScoreByName(allCisPairs, captureHiC[["raw"]], replaceZeroByNA=TRUE)
    allCisPairs[,"captureC_ObsExp"] = getPairwiseMatrixScoreByName(allCisPairs, captureHiC[["obsExp"]], replaceZeroByNA=TRUE)
    
    # add common compartment and subcompartment
    paralogPairsUniq <- addCommonCompartment(paralogPairsUniq, tssGR, compGR, subCompGR)
    allCisPairs <- addCommonCompartment(allCisPairs, tssGR, compGR, subCompGR)
    
    # add information of one-two-one orthologs in other species
    for (orgStr in names(orgStr2Name)){
    
        allCisPairs = addOrthologAnnotation(allCisPairs, orthologsSpeciesList[[orgStr]], orgStr, speciesTssGR[[orgStr]], speciesTADs[[orgStr]], speciesHiC[[orgStr]][[1]], speciesHiC[[orgStr]][[2]] )
        Sys.sleep(3) # hack to fix problems with bplapply on MOGON
    
        # add age related boolean flags
        allCisPairs <- addAgeFlags(allCisPairs, orthologsSpeciesList[[orgStr]], orgStr)
    
    }
    
    # add duplication age
    allCisPairs <- addAge(allCisPairs, pair.ac, uniPSwissToEnsgDF)
    
    # add HIPPIE
    allCisPairs <- addHIPPIE(allCisPairs, hippie)

    # save allCisPairs with all annotations
    write.table(allCisPairs, file=paste0(outPrefix, ".paralog_pairs.allCisPairs.annotated.txt"),
        sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    
    #~ allCisPairs <- read.table(paste0(outPrefix, ".paralog_pairs.allCisPairs.annotated.txt"), header=TRUE, sep="\t")
    #~ allCisPairs[,1] <- as.character(allCisPairs[,1])
    #~ allCisPairs[,2] <- as.character(allCisPairs[,2])
    
    #-----------------------------------------------------------------------
    # Filter for close and distal pairs
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
    # - sampCisPairs        := sampled by distance of allCisPairs
    # - sampEhStrPairs      := sampled by distance, numb. of enhancer, and strand of allCisPairs
    message("INFO: Start sampling of gene pairs...")
    
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
    # Sample cis pairs according to linear distance in paralog gene pairs
    #-----------------------------------------------------------------------
    # get all possible gene pairs within MAX_DIST bp
    allCloseGenePairsOVL <- getAllGenePairs(tssGR, maxDist=MAX_DIST, minDist=1)
    
    # filter out overlapping pairs
    closeNonOVL <- nonOverlappingGenePairs(allCloseGenePairsOVL, genesGR, useIDs=TRUE)
    allCloseGenePairs <- allCloseGenePairsOVL[closeNonOVL,]
    
    # get sample weights according to distance
    closeDistBreaks <- seq(log10(1), log10(MAX_DIST), length.out=N_SAMPLING_BIN+1)
    closeWeights <- weightsByBin(log10(abs(closePairs$dist)), log10(abs(allCloseGenePairs$dist)), breaks=closeDistBreaks)
    
    # sample close pairs according to distance weight
    sampClosePairs <- bplapply(1:N_RAND, function(x){ 
        sampleFromAllPairsByWeight(n=nrow(closePairs), hitDF=allCloseGenePairs, tssGR, weight=closeWeights)
        })
    Sys.sleep(3) # hack to fix problems with bplapply on MOGON
    
    #-----------------------------------------------------------------------
    # Sample cis pairs according to enhancer number in paralogs and linear distance in paralog gene pairs
    #-----------------------------------------------------------------------
    # get sample weights according to enahncer number and distance
    closeAndEhWeights <- weightsByBinDistAndEnhancers(distWeight=closeWeights, sourcePairs=closePairs, hitDF=allCloseGenePairs, tssGR)
    
    # sample according to enahncer number and distance
    sampEhClosePairs <- bplapply(1:N_RAND, function(x){ 
        sampleFromAllPairsByWeight(n=nrow(closePairs), hitDF=allCloseGenePairs, tssGR, weight=closeAndEhWeights)
        })
    Sys.sleep(3) # hack to fix problems with bplapply on MOGON
    
    #-----------------------------------------------------------------------
    # Sample all cis pairs by distance only
    #-----------------------------------------------------------------------
    # Now sample from all possible gene pairs within 1 - DISTAL_MAX_DIST bp
    allCisGenePairsOVL <- getAllGenePairs(tssGR, maxDist=DISTAL_MAX_DIST, minDist=1)
    
    # filter out overlapping pairs
    allNonOVL <- nonOverlappingGenePairs(allCisGenePairsOVL, genesGR, useIDs=TRUE)
    allCisGenePairs <- allCisGenePairsOVL[allNonOVL,]
    
    # get sampling weights for distal pairs according to distance
    allBreaks <- seq(log10(1), log10(DISTAL_MAX_DIST), length.out=N_SAMPLING_BIN+1)
    allWeight <- weightsByBin(log10(abs(allCisPairs$dist)), log10(abs(allCisGenePairs$dist)), breaks=allBreaks)
    
    # sample according to distance
    sampCisPairs <- bplapply(1:N_RAND, function(x){ 
        sampleFromAllPairsByWeight(n=nrow(allCisPairs), hitDF=allCisGenePairs, tssGR, weight=allWeight)
        })
    Sys.sleep(3) # hack to fix problems with bplapply on MOGON

    #-----------------------------------------------------------------------
    # Sample all cis pairs by distance using 100 bins 
    #-----------------------------------------------------------------------

    # get sampling weights for distal pairs according to distance
    all100Breaks <- seq(log10(1), log10(DISTAL_MAX_DIST), length.out=91)
    all100Weight <- weightsByBin(log10(abs(allCisPairs$dist)), log10(abs(allCisGenePairs$dist)), breaks=all100Breaks)
    
    # sample according to distance
    sampCis100Pairs <- bplapply(1:N_RAND, function(x){ 
        sampleFromAllPairsByWeight(n=nrow(allCisPairs), hitDF=allCisGenePairs, tssGR, weight=all100Weight)
        })
    Sys.sleep(3) # hack to fix problems with bplapply on MOGON

#~     #-----------------------------------------------------------------------
#~     # Sample from all cis pairs according to distance, enhancer, and strand
#~     #-----------------------------------------------------------------------
#~     # get sample weights according to enahncer number and distance
#~     distEhStrWeights <- weightsByBinDistAndEnhancers(distWeight=closeWeights, sourcePairs=closePairs, hitDF=allCloseGenePairs, tssGR)
#~     
#~     # sample according to enahncer number and distance
#~     sampEhClosePairs <- bplapply(1:N_RAND, function(x){ 
#~         sampleFromAllPairsByWeight(n=nrow(closePairs), hitDF=allCloseGenePairs, tssGR, weight=closeAndEhWeights)
#~         })
#~     Sys.sleep(3) # hack to fix problems with bplapply on MOGON
#~     


    #-----------------------------------------------------------------------
    # Sample distal cis pairs by distance only
    #-----------------------------------------------------------------------
    # Now sample from all possible gene pairs within DISTAL_MIN_DIST - DISTAL_MAX_DIST bp
    allDistalGenePairs <- allCisGenePairs[abs(allCisGenePairs$dist) > DISTAL_MIN_DIST,]
    
    # get sampling weights for distal pairs according to distance
    distalBreaks <- seq(log10(DISTAL_MIN_DIST), log10(DISTAL_MAX_DIST), length.out=N_SAMPLING_BIN+1)
    distalWeight <- weightsByBin(log10(abs(distalPairs$dist)), log10(abs(allDistalGenePairs$dist)), breaks=distalBreaks)
    
    # sample according to distance
    sampDistalPairs <- bplapply(1:N_RAND, function(x){ 
        sampleFromAllPairsByWeight(n=nrow(distalPairs), hitDF=allDistalGenePairs, tssGR, weight=distalWeight)
        })
    Sys.sleep(3) # hack to fix problems with bplapply on MOGON

    
    #-----------------------------------------------------------------------
    # save sampled gene pairs before annotation.
    #-----------------------------------------------------------------------
    save(randPairs, randPairsInCis, sampClosePairs, sampEhClosePairs, sampDistalPairs, sampCisPairs, sampCis100Pairs, file=paste0(WORKIMAGE_FILE, ".sampled_pairs_before_annotation.Rdata"))
    #load(paste0(WORKIMAGE_FILE, ".sampled_pairs_before_annotation.Rdata"))
    message("INFO: Finished sampling of gene pairs.")
    
    #-----------------------------------------------------------------------
    # annotate sampled gene pairs
    #-----------------------------------------------------------------------
    
    # annotate gene symbols and strand
    sampClosePairs <- lapply(sampClosePairs, addHGNC, tssGR)
    sampEhClosePairs <- lapply(sampEhClosePairs, addHGNC, tssGR)
    sampDistalPairs <- lapply(sampDistalPairs, addHGNC, tssGR)
    sampCisPairs <- lapply(sampCisPairs, addHGNC, tssGR)
    sampCis100Pairs <- lapply(sampCis100Pairs, addHGNC, tssGR)
    
    sampClosePairs <- lapply(sampClosePairs, addSameStrand, tssGR)
    sampEhClosePairs <- lapply(sampEhClosePairs, addSameStrand, tssGR)
    sampDistalPairs <- lapply(sampDistalPairs, addSameStrand, tssGR)
    sampCisPairs <- lapply(sampCisPairs, addSameStrand, tssGR)
    sampCis100Pairs <- lapply(sampCis100Pairs, addSameStrand, tssGR)
    
    # annotate common enhancers
    sampEhClosePairs <- bplapply(sampEhClosePairs, addCommonEnhancer, gene2ehID)
    Sys.sleep(3) # hack to fix problems with bplapply on MOGON
    sampEhClosePairs <- bplapply(sampEhClosePairs, addRelativeEnhancerPosition, tssGR, gene2ehID, ehGR)
    Sys.sleep(3) # hack to fix problems with bplapply on MOGON
    
    # add expression in IMR90 cells
    sampClosePairs <- lapply(sampClosePairs, addPairExp, expDFlist[["ENCODE_cell_lines"]], expCol="IMR_90", label="exp_IMR90")
    sampEhClosePairs <- lapply(sampEhClosePairs, addPairExp, expDFlist[["ENCODE_cell_lines"]], expCol="IMR_90", label="exp_IMR90")
    sampDistalPairs <- lapply(sampDistalPairs, addPairExp, expDFlist[["ENCODE_cell_lines"]], expCol="IMR_90", label="exp_IMR90")
    sampCisPairs <- lapply(sampCisPairs, addPairExp, expDFlist[["ENCODE_cell_lines"]], expCol="IMR_90", label="exp_IMR90")
    sampCis100Pairs <- lapply(sampCis100Pairs, addPairExp, expDFlist[["ENCODE_cell_lines"]], expCol="IMR_90", label="exp_IMR90")
    
    # add expression correlation
    for (expName in names(expDFlist)) {
        
        message(paste("INFO: Annotate sampled pairs with expression form:", expName))
        expDF <- expDFlist[[expName]]
    
        sampClosePairs <- bplapply(sampClosePairs, addCor, expDF, colName=paste0(expName, "_expCor"))
        Sys.sleep(3) # hack to fix problems with bplapply on MOGON
    
        sampEhClosePairs <- bplapply(sampEhClosePairs, addCor, expDF, colName=paste0(expName, "_expCor"))
        Sys.sleep(3) # hack to fix problems with bplapply on MOGON
    
        sampDistalPairs <- bplapply(sampDistalPairs, addCor, expDF, colName=paste0(expName, "_expCor"))
        Sys.sleep(3) # hack to fix problems with bplapply on MOGON
    
        sampCisPairs <- bplapply(sampCisPairs, addCor, expDF, colName=paste0(expName, "_expCor"))
        Sys.sleep(3) # hack to fix problems with bplapply on MOGON

        sampCis100Pairs <- bplapply(sampCis100Pairs, addCor, expDF, colName=paste0(expName, "_expCor"))
        Sys.sleep(3) # hack to fix problems with bplapply on MOGON
        
    }
    
    # save temp
    #~ save.image(paste0(WORKIMAGE_FILE, ".sampling_and_after_expression_annotation.Rdata"))
    
    # annotate with orthologs in other species:
    for (orgStr in names(orgStr2Name)){
    
        # add ortholog information to random sampled pairs on same chromosome
        randPairsInCis <- bplapply(randPairsInCis, addOrthologAnnotation, orthologsSpeciesList[[orgStr]], orgStr, speciesTssGR[[orgStr]], speciesTADs[[orgStr]], speciesHiC[[orgStr]][[1]], speciesHiC[[orgStr]][[2]], inParallel=FALSE)
    
        sampClosePairs = bplapply(sampClosePairs, addOrthologAnnotation, orthologsSpeciesList[[orgStr]], orgStr, speciesTssGR[[orgStr]], speciesTADs[[orgStr]], speciesHiC[[orgStr]][[1]], speciesHiC[[orgStr]][[2]], inParallel=FALSE)
        Sys.sleep(3) # hack to fix problems with bplapply on MOGON
    
        sampEhClosePairs = bplapply(sampEhClosePairs, addOrthologAnnotation, orthologsSpeciesList[[orgStr]], orgStr, speciesTssGR[[orgStr]], speciesTADs[[orgStr]], inParallel=FALSE)
        Sys.sleep(3) # hack to fix problems with bplapply on MOGON
    
        sampDistalPairs = bplapply(sampDistalPairs, addOrthologAnnotation, orthologsSpeciesList[[orgStr]], orgStr, speciesTssGR[[orgStr]], speciesTADs[[orgStr]], speciesHiC[[orgStr]][[1]], speciesHiC[[orgStr]][[2]], inParallel=FALSE)
        Sys.sleep(3) # hack to fix problems with bplapply on MOGON

        sampCisPairs = bplapply(sampCisPairs, addOrthologAnnotation, orthologsSpeciesList[[orgStr]], orgStr, speciesTssGR[[orgStr]], speciesTADs[[orgStr]], speciesHiC[[orgStr]][[1]], speciesHiC[[orgStr]][[2]], inParallel=FALSE)
        Sys.sleep(3) # hack to fix problems with bplapply on MOGON
        
        sampCis100Pairs = bplapply(sampCis100Pairs, addOrthologAnnotation, orthologsSpeciesList[[orgStr]], orgStr, speciesTssGR[[orgStr]], speciesTADs[[orgStr]], speciesHiC[[orgStr]][[1]], speciesHiC[[orgStr]][[2]], inParallel=FALSE)
        Sys.sleep(3) # hack to fix problems with bplapply on MOGON
        
        # add age related boolean flags
        sampClosePairs <- bplapply(sampClosePairs, addAgeFlags, orthologsSpeciesList[[orgStr]], orgStr )
        Sys.sleep(3) # hack to fix problems with bplapply on MOGON
    
        sampEhClosePairs <- bplapply(sampEhClosePairs, addAgeFlags, orthologsSpeciesList[[orgStr]], orgStr )
        Sys.sleep(3) # hack to fix problems with bplapply on MOGON
    
        sampDistalPairs <- bplapply(sampDistalPairs, addAgeFlags, orthologsSpeciesList[[orgStr]], orgStr )
        Sys.sleep(3) # hack to fix problems with bplapply on MOGON
        
        sampCisPairs <- bplapply(sampCisPairs, addAgeFlags, orthologsSpeciesList[[orgStr]], orgStr )
        Sys.sleep(3) # hack to fix problems with bplapply on MOGON
        
        sampCis100Pairs <- bplapply(sampCis100Pairs, addAgeFlags, orthologsSpeciesList[[orgStr]], orgStr )
        Sys.sleep(3) # hack to fix problems with bplapply on MOGON
        
    }    
    
    # add annotation of beeing located in the same TAD
    
    # make GRanges objects for cis paralog pairs and random paris on same chromosome
    sampClosePairsGR = bplapply(sampClosePairs, getPairAsGR, tssGR)
    Sys.sleep(3) # hack to fix problems with bplapply on MOGON
    sampEhClosePairsGR = bplapply(sampEhClosePairs, getPairAsGR, tssGR)
    Sys.sleep(3) # hack to fix problems with bplapply on MOGON
    sampCisPairsGR = bplapply(sampCisPairs, getPairAsGR, tssGR)
    Sys.sleep(3) # hack to fix problems with bplapply on MOGON
    sampCis100PairsGR = bplapply(sampCis100Pairs, getPairAsGR, tssGR)
    Sys.sleep(3) # hack to fix problems with bplapply on MOGON
    
    for(tadName in names(allTADs)){
    
        message(paste("INFO: Compute overlap with TADs from:", tadName))
        # co-occurance within the same domain
        sampClosePairsGR = bplapply(sampClosePairsGR, addWithinSubject, allTADs[[tadName]], tadName)
        Sys.sleep(3) # hack to fix problems with bplapply on MOGON
        
        sampClosePairs <- bplapply(sampClosePairs, addSubTADmode, allTADs[[tadName]], tssGR, paste0(tadName, "_subTAD"))
        Sys.sleep(3) # hack to fix problems with bplapply on MOGON
        
        sampEhClosePairsGR = bplapply(sampEhClosePairsGR, addWithinSubject, allTADs[[tadName]], tadName)
        Sys.sleep(3) # hack to fix problems with bplapply on MOGON
        sampEhClosePairs <- bplapply(sampEhClosePairs, addSubTADmode, allTADs[[tadName]], tssGR, paste0(tadName, "_subTAD"))
        Sys.sleep(3) # hack to fix problems with bplapply on MOGON
    
        sampCisPairsGR = bplapply(sampCisPairsGR, addWithinSubject, allTADs[[tadName]], tadName)
        Sys.sleep(3) # hack to fix problems with bplapply on MOGON
        sampCisPairs <- bplapply(sampCisPairs, addSubTADmode, allTADs[[tadName]], tssGR, paste0(tadName, "_subTAD"))
        Sys.sleep(3) # hack to fix problems with bplapply on MOGON
        
        sampCis100PairsGR = bplapply(sampCis100PairsGR, addWithinSubject, allTADs[[tadName]], tadName)
        Sys.sleep(3) # hack to fix problems with bplapply on MOGON
        sampCis100Pairs <- bplapply(sampCis100Pairs, addSubTADmode, allTADs[[tadName]], tssGR, paste0(tadName, "_subTAD"))
        Sys.sleep(3) # hack to fix problems with bplapply on MOGON

    }
    
    # assign annotation in GRanges object to gene pair data.frames
    for (i in seq(N_RAND)){
        sampClosePairs[[i]][,names(allTADs)] <- data.frame( mcols(sampClosePairsGR[[i]])[,names(allTADs)] )
        sampEhClosePairs[[i]][,names(allTADs)] <- data.frame( mcols(sampEhClosePairsGR[[i]])[,names(allTADs)] )
        sampCisPairs[[i]][,names(allTADs)] <- data.frame( mcols(sampCisPairsGR[[i]])[,names(allTADs)] )
        sampCis100Pairs[[i]][,names(allTADs)] <- data.frame( mcols(sampCis100PairsGR[[i]])[,names(allTADs)] )
    }
    
    # DEBUG SAVE
#~     save(sampClosePairs, sampEhClosePairs, sampDistalPairs, sampCisPairs, sampCis100Pairs, file=paste0(WORKIMAGE_FILE, ".sampled_pairs_before_annotation.Rdata"))
    #load(paste0(WORKIMAGE_FILE, ".sampled_pairs_before_Hi-C_annotation.Rdata"))
    
    # Adds Hi-C contact frequencies to a gene pair data set
    for (i in 1:N_RAND) {
        message(paste("Work on sample:", i))

        sampClosePairs[[i]] <- addHiCfreq(sampClosePairs[[i]], tssGR, HiClistRaw, label="HiCRaw", ignoreSameBin=TRUE)
        Sys.sleep(3) # hack to fix problems with bplapply on MOGON
        sampClosePairs[[i]] <- addHiCfreq(sampClosePairs[[i]], tssGR, HiClist, label="HiC", ignoreSameBin=TRUE)
        Sys.sleep(3) # hack to fix problems with bplapply on MOGON
        sampClosePairs[[i]] <- addHiCobsExp(sampClosePairs[[i]], tssGR, expectedHiCList, HIC_RESOLUTION, HiClabel="HiC", label="HiCobsExp")

        sampDistalPairs[[i]] <- addHiCfreq(sampDistalPairs[[i]], tssGR, HiClistRaw, label="HiCRaw", ignoreSameBin=TRUE)
        Sys.sleep(3) # hack to fix problems with bplapply on MOGON
        sampDistalPairs[[i]] <- addHiCfreq(sampDistalPairs[[i]], tssGR, HiClist, label="HiC", ignoreSameBin=TRUE)
        Sys.sleep(3) # hack to fix problems with bplapply on MOGON
        sampDistalPairs[[i]] <- addHiCobsExp(sampDistalPairs[[i]], tssGR, expectedHiCList, HIC_RESOLUTION, HiClabel="HiC", label="HiCobsExp")

        sampCisPairs[[i]] <- addHiCfreq(sampCisPairs[[i]], tssGR, HiClistRaw, label="HiCRaw", ignoreSameBin=TRUE)
        Sys.sleep(3) # hack to fix problems with bplapply on MOGON
        sampCisPairs[[i]] <- addHiCfreq(sampCisPairs[[i]], tssGR, HiClist, label="HiC", ignoreSameBin=TRUE)
        Sys.sleep(3) # hack to fix problems with bplapply on MOGON
        sampCisPairs[[i]] <- addHiCobsExp(sampCisPairs[[i]], tssGR, expectedHiCList, HIC_RESOLUTION, HiClabel="HiC", label="HiCobsExp")

        sampCis100Pairs[[i]] <- addHiCfreq(sampCis100Pairs[[i]], tssGR, HiClistRaw, label="HiCRaw", ignoreSameBin=TRUE)
        Sys.sleep(3) # hack to fix problems with bplapply on MOGON
        sampCis100Pairs[[i]] <- addHiCfreq(sampCis100Pairs[[i]], tssGR, HiClist, label="HiC", ignoreSameBin=TRUE)
        Sys.sleep(3) # hack to fix problems with bplapply on MOGON
        sampCis100Pairs[[i]] <- addHiCobsExp(sampCis100Pairs[[i]], tssGR, expectedHiCList, HIC_RESOLUTION, HiClabel="HiC", label="HiCobsExp")
    }
        
    message("INFO: Finished Hi-C on sampled pairs")
    
    # add capture C
    sampClosePairs <- bplapply(1:N_RAND, function(i){
        sampClosePairs[[i]][,"captureC_raw"] <- getPairwiseMatrixScoreByName(sampClosePairs[[i]], captureHiC[["raw"]], replaceZeroByNA=TRUE)
        sampClosePairs[[i]][,"captureC_ObsExp"] <- getPairwiseMatrixScoreByName(sampClosePairs[[i]], captureHiC[["obsExp"]], replaceZeroByNA=TRUE)
        return(sampClosePairs[[i]])
    })
    Sys.sleep(3) # hack to fix problems with bplapply on MOGON
    
    sampDistalPairs <- bplapply(1:N_RAND, function(i){
        sampDistalPairs[[i]][,"captureC_raw"] <- getPairwiseMatrixScoreByName(sampDistalPairs[[i]], captureHiC[["raw"]], replaceZeroByNA=TRUE)
        sampDistalPairs[[i]][,"captureC_ObsExp"] <- getPairwiseMatrixScoreByName(sampDistalPairs[[i]], captureHiC[["obsExp"]], replaceZeroByNA=TRUE)
        return(sampDistalPairs[[i]])
    })
    Sys.sleep(3) # hack to fix problems with bplapply on MOGON
    
    sampCisPairs <- bplapply(1:N_RAND, function(i){
        sampCisPairs[[i]][,"captureC_raw"] <- getPairwiseMatrixScoreByName(sampCisPairs[[i]], captureHiC[["raw"]], replaceZeroByNA=TRUE)
        sampCisPairs[[i]][,"captureC_ObsExp"] <- getPairwiseMatrixScoreByName(sampCisPairs[[i]], captureHiC[["obsExp"]], replaceZeroByNA=TRUE)
        return(sampCisPairs[[i]])
    })
    Sys.sleep(3) # hack to fix problems with bplapply on MOGON
    
    sampCis100Pairs <- bplapply(1:N_RAND, function(i){
        sampCis100Pairs[[i]][,"captureC_raw"] <- getPairwiseMatrixScoreByName(sampCis100Pairs[[i]], captureHiC[["raw"]], replaceZeroByNA=TRUE)
        sampCis100Pairs[[i]][,"captureC_ObsExp"] <- getPairwiseMatrixScoreByName(sampCis100Pairs[[i]], captureHiC[["obsExp"]], replaceZeroByNA=TRUE)
        return(sampCis100Pairs[[i]])
    })
    Sys.sleep(3) # hack to fix problems with bplapply on MOGON
    
    message("INFO: Finished Capture-C on sampled pairs")
    
    # add common compartment and subcompartment
    randPairs <- lapply(randPairs, addCommonCompartment, tssGR, compGR, subCompGR)
    sampClosePairs <- lapply(sampClosePairs, addCommonCompartment, tssGR, compGR, subCompGR)
    sampEhClosePairs <- lapply(sampEhClosePairs, addCommonCompartment, tssGR, compGR, subCompGR)
    sampDistalPairs <- lapply(sampDistalPairs, addCommonCompartment, tssGR, compGR, subCompGR)
    sampCisPairs <- lapply(sampCisPairs, addCommonCompartment, tssGR, compGR, subCompGR)
    sampCis100Pairs <- lapply(sampCis100Pairs, addCommonCompartment, tssGR, compGR, subCompGR)
    
    # add duplication age
    randPairs <- bplapply(randPairs, addAge, pair.ac, uniPSwissToEnsgDF)
    Sys.sleep(3) # hack to fix problems with bplapply on MOGON
    sampClosePairs <- bplapply(sampClosePairs, addAge, pair.ac, uniPSwissToEnsgDF)
    Sys.sleep(3) # hack to fix problems with bplapply on MOGON
    sampEhClosePairs <- bplapply(sampEhClosePairs, addAge, pair.ac, uniPSwissToEnsgDF)
    Sys.sleep(3) # hack to fix problems with bplapply on MOGON
    sampDistalPairs <- bplapply(sampDistalPairs, addAge, pair.ac, uniPSwissToEnsgDF)
    Sys.sleep(3) # hack to fix problems with bplapply on MOGON
    sampCisPairs <- bplapply(sampCisPairs, addAge, pair.ac, uniPSwissToEnsgDF)
    Sys.sleep(3) # hack to fix problems with bplapply on MOGON
    sampCis100Pairs <- bplapply(sampCis100Pairs, addAge, pair.ac, uniPSwissToEnsgDF)
    Sys.sleep(3) # hack to fix problems with bplapply on MOGON

    # add HIPPIE
    randPairs <- lapply(randPairs, addHIPPIE, hippie)
    sampClosePairs <- lapply(sampClosePairs, addHIPPIE, hippie)
    sampEhClosePairs <-lapply(sampEhClosePairs, addHIPPIE, hippie)
    sampDistalPairs <- lapply(sampDistalPairs, addHIPPIE, hippie)
    sampCisPairs <- lapply(sampCisPairs, addHIPPIE, hippie)
    sampCis100Pairs <- lapply(sampCis100Pairs, addHIPPIE, hippie)
    
    save(randPairs, sampClosePairs, sampEhClosePairs, sampDistalPairs, sampCisPairs, sampCis100Pairs, file=paste0(WORKIMAGE_FILE, ".sampled_pairs_after_annotation.Rdata"))
    #load(paste0(WORKIMAGE_FILE, ".sampled_pairs_after_annotation.Rdata"))

    # combine all sampling replicates to one data frame
    randPairsCombined <- do.call("rbind", randPairs)
    randPairsInCisCombined <- do.call("rbind", randPairsInCis)
    sampClosePairsCombined <- do.call("rbind", sampClosePairs)
    sampEhClosePairsCombined <- do.call("rbind", sampEhClosePairs)
    sampDistalPairsCombined <- do.call("rbind", sampDistalPairs)
    sampCisPairsCombined <- do.call("rbind", sampCisPairs)
    sampCis100PairsCombined <- do.call("rbind", sampCis100Pairs)
    
    #-----------------------------------------------------------------------
    # save a work image after sampling and annotation.
    #-----------------------------------------------------------------------
    save.image(paste0(WORKIMAGE_FILE, ".sampling_and_annotation.Rdata"))
#~     stop("INFO: Stopped script HERE!")

}else{
    load(paste0(WORKIMAGE_FILE, ".sampling_and_annotation.Rdata"))
}


#=======================================================================
# Build data frames for analysis
#=======================================================================

#-----------------------------------------------------------------------
# Make data.frame for plotting with ggplot2 for all cis pairs
#-----------------------------------------------------------------------

breaksExp <- c(0, 1, 10, Inf)
breaksExpCor <- c(-1, -.5, 0, .5, 1)
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

# write data.frame as table to output file
write.table(plotDFcis100HiC, file=paste0(outPrefix, ".plotDFcis100HiC.csv"),
    sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

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

# write data.frame as table to output file
write.table(plotDFeh, file=paste0(outPrefix, ".plotDFeh.csv"),
    sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

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
# allDF for sampling by enhancers
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
# save workspace image
#=======================================================================
save.image(WORKIMAGE_FILE)
#INFO: FANTOM enhancer map: 65953 of 66942 enhancer-associated genes could be mapped uniquelly to ENSG ID
