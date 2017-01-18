########################################################################
#
# A script to analyse the co-regulation of paralog genes.
# It looks for association (compared to sampled gene pairs) 
# for a paralog pair to
# - co-localization on linear genome in same chromosome and short distance between genes
# - common enhancers associations in correlation based maps
# - co-occurances in same TAD
# - more contact to each other in Hi-C map or more contacts to common enhancers.
# - other related features
########################################################################

require(stringr)        # for some string functionality
require(GenomicRanges)  # for genomic intervals and overlap calculation
require(rtracklayer)    # for import.bed
require(BiocParallel)   # for parallel computing

#-----------------------------------------------------------------------
# Load parameters from external script
#-----------------------------------------------------------------------
# read parameter script location from command-line argument

#~ args <- commandArgs(trailingOnly = TRUE)
#~ PARAM_SCRIPT=args[1]
PARAM_SCRIPT="R/paralog_regulation.param.R"

message(paste("INFO: Load parameter from this script:", PARAM_SCRIPT))
print(PARAM_SCRIPT)

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
source("R/functions.regMap.R")
source("R/functions.GRanges.R")
source("R/functions.genePairs.R")
source("R/functions.genePairs.randomization.R")
source("R/functions.genePairs.paralog_analysis.R")
source("R/parseHiC.R")

# load ensemble data sets of genes and paralog pairs
source("R/data.ensembl.R")     # load ensambl data
source("R/data.expression.R")  # load expression data from EBI expression atlas
source("R/data.cohesinKO.R")  # load cohesin KO expression data
# combine expression data sets into one list:

source("R/data.captureHiC.R")  # load capture Hi-C data between promoters from Mifsud et al. 2015



if (!LOAD_PAIRS) {
    
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
            "stable_TADs"=getConservedTADs(RaoTADs, n=4, fraction=.9)
        )        
        
        allTADs = c(
            RaoTADs,
            stableTADs,
            DixonTADs
        )
        
        allBoundaries <- lapply(allTADs, getBoundaries)
        
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
                
        #-----------------------------------------------------------------------
        # Load Hi-C data from Rao et al. 2014
        #-----------------------------------------------------------------------
        if ( !USE_LOCAL_HIC) {
        
            # parse normalized Hi-C map from Rao et al. 
            HiClist = parseRaoHiC(CELL, HIC_RESOLUTION, HIC_DATA_DIR, seqInfo)
    
            HiClistRaw = parseRaoHiC(CELL, HIC_RESOLUTION, HIC_DATA_DIR, seqInfo, normStr=NULL)
            
            # parse expected contacts by distance for each chromosome
            expectedHiCList <- parseRaoExpected(CELL, HIC_RESOLUTION, HIC_DATA_DIR)
            
            # save data for faster loading next time
            save(HiClist, HiClistRaw, expectedHiCList, file=paste0(outDataPrefix, ".Hi-C.", CELL, ".", HIC_RESOLUTION, ".HiClist_expectedHiCList.RData"))
                    
        }else{
            load(paste0(outDataPrefix, ".Hi-C.", CELL, ".", HIC_RESOLUTION, ".HiClist_expectedHiCList.RData"))
        }
        
        #-----------------------------------------------------------------------
        # Parse Hi-C data for other species
        #-----------------------------------------------------------------------
        
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
    
        hippie$ID <- getPairIDsorted(cbind(entrezToENSG[as.character(hippieDF[,2])], entrezToENSG[as.character(hippieDF[,4])]))

        #-----------------------------------------------------------------------
        # Parse house keeping genes 
        #-----------------------------------------------------------------------
        housekeepingDF <- read.table("data/Eisenberg2013/HK_genes.txt", header=FALSE, sep="\t", stringsAsFactors=FALSE)
        housekeepingGenes <- refSeqToENSG[housekeepingDF[,2]]

        #-----------------------------------------------------------------------
        # save all loaded data as image file
        #-----------------------------------------------------------------------
        save.image(paste0(outDataPrefix, ".loaded_data.Rdata"))
    }else{
        load(paste0(outDataPrefix, ".loaded_data.Rdata"))
    }

    #=======================================================================
    # 2.) Filter and annotate paralog gene pairs
    #=======================================================================
        
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
        
    # get for each gene only one unique pair, the one with highest similarity
    # this is computed by an maximum weight matching
    paralogPairsUniqPnonOVLWithDS = paralogPairsUniqPnonOVL[!is.na(paralogPairsUniqPnonOVL[,"hsapiens_paralog_ds"]),]
    
    # choose factor to be multiplied with gene pair substitution rate which will be maximized during pair selection
    pairSelectFact = ifelse(SELECT_OLD_PAIRS, 1, -1)
    paralogPairsUniqPnonOVLUniqG = uniquePairPerGeneBySim(paralogPairsUniqPnonOVLWithDS, pairSelectFact*paralogPairsUniqPnonOVLWithDS[,"hsapiens_paralog_ds"])
    
    # get only a unique pair order (one of A-B, B-A) form the unique pairs
    paralogPairsUniq = uniquePair(paralogPairsUniqPnonOVLUniqG)

    # write unique paralog pairs per gene to output file:
    write.table(paralogPairsUniq, file=paste0(outPrefix, ".paralog_pairs.paralogPairsUniq.txt"),
        sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    
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
        message(paste("INFO: Compute overlap with TADs and boundaries from:", tadName))
        # co-occurance within the same domain
        allCisPairsGR = addWithinSubject(allCisPairsGR, allTADs[[tadName]], tadName)
        
        # separated by a TAD boundary?
		separated <- countOverlaps(allCisPairsGR, allBoundaries[[tadName]]) >= 1
        mcols(allCisPairsGR)[, paste0("Boundary_", tadName)] <- separated
    }
  
    # assign annotation in GRanges object to gene pair data.frames
    allCisPairs[,names(allTADs)] <- data.frame( mcols(allCisPairsGR)[, names(allTADs)] )
    allCisPairs[,paste0("Boundary_", names(allTADs))] <- data.frame( mcols(allCisPairsGR)[, paste0("Boundary_", names(allTADs))] )

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
        
    # add HIPPIE
    allCisPairs <- addHIPPIE(allCisPairs, hippie)

    # add housekeeping gene information
    allCisPairs <- addInGeneSet(allCisPairs, housekeepingGenes, "housekeeping")
        
    # save allCisPairs with all annotations
    write.table(allCisPairs, file=paste0(outPrefix, ".paralog_pairs.allCisPairs.annotated.txt"),
        sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
        
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
        
    #=======================================================================
    # 3.) Sample random control/background data sets
    #=======================================================================
    # There will be the following types of sampled gene pairs:
    # - randPairs               := sampled completely pairs randomly from all genes
    # - sampCisPairs            := sampled by distance of allCisPairs
    # - sampDistEhPairs         := sampled by distance and numb. of enhancer from allCisPairs
    # - sampDistEhStrandPairs   := sampled by distance and numb. of enhancer, and strand from allCisPairs
    # - sampDistEhLenPairs      := sampled by distance and numb. of enhancer, and gene length from allCisPairs
    message("INFO: Start sampling of gene pairs...")
    
    #-----------------------------------------------------------------------
    # Sample pairs with equal probability from all genes
    #-----------------------------------------------------------------------
    randPairs <- lapply(1:N_RAND, function(x){getRandomPairs(nrow(paralogPairsUniq), names(tssGR))})
    Sys.sleep(3) # hack to fix problems with bplapply on MOGON
    
    # add same strand information
    randPairs <- lapply(randPairs, addSameStrand, tssGR)
    randPairs <- lapply(randPairs, addPairDist, tssGR)
        
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

    #-----------------------------------------------------------------------
    # Sample from all cis pairs according to distance and enhancer
    #-----------------------------------------------------------------------    
    # get probability density for the number of linked enhancers for each gene
    cisEhProp <- weightPairsByEnhancers(tssGR, allCisPairs, allCisGenePairs)

    weightDistEh <- cisDistProp * cisEhProp 
    probDistEh <- weightDistEh / sum(weightDistEh)
    
    # sample according to distance, strand, and enhancers
    sampDistEhPairs <- bplapply(1:N_RAND, function(x){ 
        sampleFromAllPairsByWeight(n=nrow(allCisPairs), hitDF=allCisGenePairs, tssGR, weight=probDistEh)
        })
    Sys.sleep(3) # hack to fix problems with bplapply on MOGON

    #-----------------------------------------------------------------------
    # Sample from all cis pairs according to distance, enhancer, and strand
    #-----------------------------------------------------------------------    
    # annotate all pairs with sameStrand information
    allCisGenePairs[,"sameStrand"] <- strand(tssGR)[allCisGenePairs[,1]] == strand(tssGR)[allCisGenePairs[,2]]
        
    # get weight for strand information
    strandWeight <- weightsByFactorFreq(allCisPairs[,"sameStrand"], allCisGenePairs[,"sameStrand"])
    
    weightDistEhStrand <- cisDistProp * cisEhProp * strandWeight 
    probDistEhStrand <- weightDistEhStrand / sum(weightDistEhStrand)
    
    # sample according to distance, strand, and enhancers
    sampDistEhStrandPairs <- bplapply(1:N_RAND, function(x){ 
        sampleFromAllPairsByWeight(n=nrow(allCisPairs), hitDF=allCisGenePairs, tssGR, weight=probDistEhStrand)
        })
    Sys.sleep(3) # hack to fix problems with bplapply on MOGON

    #-----------------------------------------------------------------------
    # Sample from all cis pairs according to distance, enhancer, and gene_length
    #-----------------------------------------------------------------------
    # get probability density for gene length
    paraLen1 <- log10(tssGR[allCisPairs[,1]]$gene_size)
    paraLen2 <- log10(tssGR[allCisPairs[,2]]$gene_size)
    allLen1 <- log10( tssGR$gene_size[allCisGenePairs[,1]] ) 
    allLen2 <- log10( tssGR$gene_size[allCisGenePairs[,2]] ) 
    
    # get sampling weight for each gene in pair separately
    popGeneLength1 <- weightsByBin(paraLen1, allLen1, breaks=20)
    popGeneLength2 <- weightsByBin(paraLen2, allLen2, breaks=20)

    # combine sampling weights from both genes and normalize to sum up to 1
    weightGeneLength = popGeneLength1 * popGeneLength2 / sum(popGeneLength1 * popGeneLength2)
    
    weightDistEhLen <- cisDistProp * cisEhProp * weightGeneLength 
    probDistEhLen <- weightDistEhLen / sum(weightDistEhLen)
    
    # sample according to distance, strand, and enhancers
    sampDistEhLenPairs <- bplapply(1:N_RAND, function(x){ 
        sampleFromAllPairsByWeight(n=nrow(allCisPairs), hitDF=allCisGenePairs, tssGR, weight=probDistEhLen)
        })
    Sys.sleep(3) # hack to fix problems with bplapply on MOGON

    #-----------------------------------------------------------------------
    # save sampled gene pairs before annotation.
    #-----------------------------------------------------------------------
    save(randPairs, sampCisPairs, sampDistEhPairs, sampDistEhStrandPairs, sampDistEhLenPairs, file=paste0(WORKIMAGE_FILE, ".sampled_pairs_before_annotation.Rdata"))
    #load(paste0(WORKIMAGE_FILE, ".sampled_pairs_before_annotation.Rdata"))
    message("INFO: Finished sampling of gene pairs.")
    
    #===================================================================
    # Annotate sampled gene pairs
    #===================================================================

    # annotate random pairs with common compartment, Age, and HIPPIE
    randPairs <- lapply(randPairs, addCommonCompartment, tssGR, compGR, subCompGR)
    Sys.sleep(3) # hack to fix problems with bplapply on MOGON
    randPairs <- lapply(randPairs, addHIPPIE, hippie)
    randPairs <- lapply(randPairs, addInGeneSet, housekeepingGenes, "housekeeping")
    
    # annotate with orthologs in other species:
    for (orgStr in names(orgStr2Name)){
        
        message(paste("INFO: Annotate pairs with ortholog annotation from:", orgStr))
        
        # add ortholog information to random sampled pairs on same chromosome
        randPairs <- bplapply(randPairs, addOrthologAnnotation, orthologsSpeciesList[[orgStr]], orgStr, speciesTssGR[[orgStr]], speciesTADs[[orgStr]], speciesHiC[[orgStr]][[1]], speciesHiC[[orgStr]][[2]], inParallel=FALSE)
        Sys.sleep(3) # hack to fix problems with bplapply on MOGON
    
        # add age related boolean flags
        randPairs <- bplapply(randPairs, addAgeFlags, orthologsSpeciesList[[orgStr]], orgStr )
        Sys.sleep(3) # hack to fix problems with bplapply on MOGON

    }

    # filter for random pairs in Cis
    randPairsInCis <- lapply(randPairs, getCisPairs, tssGR)
    
    # make a list of all sampled data sets for combined annotation
    sampPairs <- list(
        "sampCisPairs"=sampCisPairs, 
        "sampDistEhPairs"=sampDistEhPairs, 
        "sampDistEhStrandPairs"=sampDistEhStrandPairs, 
        "sampDistEhLenPairs"=sampDistEhLenPairs)
    
    for (j in 1:length(sampPairs)){
    
        message(paste("INFO: Start annotating the sampled pairs:", names(sampPairs)[j]))
        
        # annotate gene symbols and strand
        sampPairs[[j]] <- lapply(sampPairs[[j]], addHGNC, tssGR)
        sampPairs[[j]] <- lapply(sampPairs[[j]], addSameStrand, tssGR)
        
        # annotate common enhancers
        sampPairs[[j]] <- bplapply(sampPairs[[j]], addCommonEnhancer, gene2ehID)
        Sys.sleep(3) # hack to fix problems with bplapply on MOGON
        sampPairs[[j]] <- bplapply(sampPairs[[j]], addRelativeEnhancerPosition, tssGR, gene2ehID, ehGR)
        Sys.sleep(3) # hack to fix problems with bplapply on MOGON
        
        # add expression in IMR90 cells
        sampPairs[[j]] <- lapply(sampPairs[[j]], addPairExp, expDFlist[["ENCODE_cell_lines"]], expCol="IMR_90", label="exp_IMR90")
        
        
        # add expression correlation
        for (expName in names(expDFlist)) {
            
            message(paste("INFO: Annotate sampled pairs with expression form:", expName))
            expDF <- expDFlist[[expName]]
        
            sampPairs[[j]] <- bplapply(sampPairs[[j]], addCor, expDF, colName=paste0(expName, "_expCor"))
            Sys.sleep(3) # hack to fix problems with bplapply on MOGON
        }
                
        # save temp
        #~ save.image(paste0(WORKIMAGE_FILE, ".sampling_and_after_expression_annotation.Rdata"))
        
        # annotate with orthologs in other species:
        for (orgStr in names(orgStr2Name)){
            
            message(paste("INFO: Annotate pairs with ortholog annotation from:", orgStr))
            
            # add ortholog information to random sampled pairs on same chromosome
            sampPairs[[j]] <- bplapply(sampPairs[[j]], addOrthologAnnotation, orthologsSpeciesList[[orgStr]], orgStr, speciesTssGR[[orgStr]], speciesTADs[[orgStr]], speciesHiC[[orgStr]][[1]], speciesHiC[[orgStr]][[2]], inParallel=FALSE)
        
            # add age related boolean flags
            sampPairs[[j]] <- bplapply(sampPairs[[j]], addAgeFlags, orthologsSpeciesList[[orgStr]], orgStr )
            Sys.sleep(3) # hack to fix problems with bplapply on MOGON

        }
        
        # add annotation of beeing located in the same TAD        
        # make GRanges objects for cis paralog pairs and random paris on same chromosome
        gplGR = bplapply(sampPairs[[j]], getPairAsGR, tssGR)
        Sys.sleep(3) # hack to fix problems with bplapply on MOGON
        
        for(tadName in names(allTADs)){
        
            message(paste("INFO: Compute overlap with TADs from:", tadName))
            # co-occurance within the same domain
            gplGR = bplapply(gplGR, addWithinSubject, allTADs[[tadName]], tadName)
            Sys.sleep(3) # hack to fix problems with bplapply on MOGON

			# separated by a TAD boundary?
            gplGR = bplapply(gplGR, function(gr){

					separated <- countOverlaps(gr, allBoundaries[[tadName]]) >= 1
			        mcols(gr)[, paste0("Boundary_", tadName)] <- separated
					return(gr)
				})
            Sys.sleep(3) # hack to fix problems with bplapply on MOGON

            sampPairs[[j]] <- bplapply(sampPairs[[j]], addSubTADmode, allTADs[[tadName]], tssGR, paste0(tadName, "_subTAD"))
            Sys.sleep(3) # hack to fix problems with bplapply on MOGON
        }
        
        # assign annotation in GRanges object to gene pair data.frames
        for (i in seq(N_RAND)){
            sampPairs[[j]][[i]][,names(allTADs)] <- data.frame( mcols(gplGR[[i]])[,names(allTADs)] )
            sampPairs[[j]][[i]][, paste0("Boundary_", tadName)] <- data.frame( mcols(gplGR[[i]])[, paste0("Boundary_", tadName)] )
        }
                
        # Adds Hi-C contact frequencies to a gene pair data set
        for (i in 1:N_RAND) {
            message(paste("Work on sample:", i))
    
            sampPairs[[j]][[i]] <- addHiCfreq(sampPairs[[j]][[i]], tssGR, HiClistRaw, label="HiCRaw", ignoreSameBin=TRUE)
            Sys.sleep(3) # hack to fix problems with bplapply on MOGON
            sampPairs[[j]][[i]] <- addHiCfreq(sampPairs[[j]][[i]], tssGR, HiClist, label="HiC", ignoreSameBin=TRUE)
            Sys.sleep(3) # hack to fix problems with bplapply on MOGON
            sampPairs[[j]][[i]] <- addHiCobsExp(sampPairs[[j]][[i]], tssGR, expectedHiCList, HIC_RESOLUTION, HiClabel="HiC", label="HiCobsExp")
    
        }
            
        message("INFO: Finished Hi-C on sampled pairs")
        
        # add capture C
        sampPairs[[j]] <- bplapply(1:N_RAND, function(i){
            sampPairs[[j]][[i]][,"captureC_raw"] <- getPairwiseMatrixScoreByName(sampPairs[[j]][[i]], captureHiC[["raw"]], replaceZeroByNA=TRUE)
            sampPairs[[j]][[i]][,"captureC_ObsExp"] <- getPairwiseMatrixScoreByName(sampPairs[[j]][[i]], captureHiC[["obsExp"]], replaceZeroByNA=TRUE)
            return(sampPairs[[j]][[i]])
        })
        Sys.sleep(3) # hack to fix problems with bplapply on MOGON
        
        message("INFO: Finished Capture-C on sampled pairs")

        # add common compartment and subcompartment
        sampPairs[[j]] <- lapply(sampPairs[[j]], addCommonCompartment, tssGR, compGR, subCompGR)
            
        # add HIPPIE
        sampPairs[[j]] <- lapply(sampPairs[[j]], addHIPPIE, hippie)

        # add house keeping gene
        sampPairs[[j]] <- lapply(sampPairs[[j]], addInGeneSet, housekeepingGenes, "housekeeping")
                    
    }

    # attach list items to search path by their names in list of sampled pairs
    sampCisPairs <- sampPairs[["sampCisPairs"]]
    sampDistEhPairs <- sampPairs[["sampDistEhPairs"]]
    sampDistEhStrandPairs <- sampPairs[["sampDistEhStrandPairs"]]
    sampDistEhLenPairs <- sampPairs[["sampDistEhLenPairs"]]

    save(randPairs, sampCisPairs, sampDistEhPairs, sampDistEhStrandPairs, sampDistEhLenPairs, file=paste0(WORKIMAGE_FILE, ".sampled_pairs_after_annotation.Rdata"))
    #load(paste0(WORKIMAGE_FILE, ".sampled_pairs_after_annotation.Rdata"))

    # combine all sampling replicates to one data frame
    randPairsCombined <- do.call("rbind", randPairs)
    randPairsInCisCombined <- do.call("rbind", randPairsInCis)
    sampCisPairsCombined <- do.call("rbind", sampCisPairs)
    sampDistEhPairsCombined <- do.call("rbind", sampDistEhPairs)
    sampDistEhStrandPairsCombined <- do.call("rbind", sampDistEhStrandPairs)
    sampDistEhLenPairsCombined <- do.call("rbind", sampDistEhLenPairs)

    #-----------------------------------------------------------------------
    # save a work image after sampling and annotation.
    #-----------------------------------------------------------------------
    save.image(paste0(WORKIMAGE_FILE, ".sampling_and_annotation.Rdata"))

}else{
    message(paste("INFO: Start loading workimage form this file:",  paste0(WORKIMAGE_FILE, ".sampling_and_annotation.Rdata")))
    load(paste0(WORKIMAGE_FILE, ".sampling_and_annotation.Rdata"))
}


#=======================================================================
# Build data.frame with all annotations for analysis
#=======================================================================

#-----------------------------------------------------------------------
# Make data.frame for plotting with ggplot2 for all cis pairs
#-----------------------------------------------------------------------

breaksExp <- c(0, 1, 10, Inf)
breaksExpCor <- c(-1, -.5, 0, .5, 1)
breaksCis = 10^(0:9) / 10^3
breaksTAD = c(0, 10, 1000, Inf)



# combine all parlogs and all randomly sampled pairs
paralogPairsUniq$g1 <- paralogPairsUniq[,1]
paralogPairsUniq$g2 <- paralogPairsUniq[,2]

commonCols <- intersect(names(paralogPairsUniq), names(randPairsCombined))

groupFactor <- factor(rep(c("paralog", "random"), c(nrow(paralogPairsUniq), nrow(randPairsCombined))), levels=c("paralog", "random"))

allPairs <- rbind.fill(paralogPairsUniq[,commonCols], randPairsCombined[,commonCols])
allPairs$group <- groupFactor
# add flag for cis pairs
allPairs$sameChrom <- !is.na(allPairs[,"dist"])

# write data.frame as table to output file
write.table(allPairs, file=paste0(outPrefix, ".allPairs.csv"),
    sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#=======================================================================
# combine all parlogs and all randomly sampled pairs
#=======================================================================

nPairs <- nrow(allCisPairs)
allCisPairs$g1 <- allCisPairs[,1]
allCisPairs$g2 <- allCisPairs[,2]

allPairDFs <- list(allCisPairs, sampCisPairsCombined, sampDistEhPairsCombined, sampDistEhStrandPairsCombined, sampDistEhLenPairsCombined)
names(allPairDFs) <- c("paralogs", "sampled by distance", "sampled by distance\nand enhancer", "sampled by distance,\nenhancer and strand", "sampled by distance,\nenhancer and gene length")

# get common columns in all data sets
cC <- Reduce(intersect, lapply(allPairDFs, names))

# combine all pairs by using only the common columns
aP <- do.call("rbind", lapply(allPairDFs, function(df) df[,cC]))

# add noZero HI-C 
aP <- addNoZero(aP)

# add further columns
aP[,"group"] <- factor(rep(c("paralog", "sampled"), c(nrow(allPairDFs[[1]]), sum(sapply(allPairDFs[2:length(allPairDFs)], nrow)))), levels=c("paralog", "sampled"))
aP[,"sampType"] <- factor(rep(names(allPairDFs), times=sapply(allPairDFs, nrow)), levels=names(allPairDFs))
aP[,"replicate"] <- c(rep(1, nPairs), replicate(length(allPairDFs)-1, rep(1:N_RAND, each=nPairs)))
aP[,"dist"] <- abs(aP[,"dist"]) / 10^3
aP[,"distGroup"] <- factor(ifelse(aP[,"dist"] <= MAX_DIST/10^3, "close", "distal"), levels=c("close", "distal"))
aP[,"dupAgeGroup"] = factor(aP[,"mmusculus_commonOrtholg"], c(TRUE, FALSE), c("Young", "Old"))
aP[,"PPI"] <- factor(aP[,"HIPPIE"] >= HIPPIE_MEDIUM, levels=c(TRUE, FALSE), labels=c("PPI", "no PPI"))
aP[,"distBin"] <- as.factor(breaksCis[.bincode(aP[,"dist"], breaksCis)])
aP[,"distTadBin"] <- factor(breaksTAD[.bincode(aP[,"dist"], breaksTAD)], levels=breaksTAD[1:3], labels=c("<10kb", "10-1000kb", ">1000kb"))
aP[,"avgExp"] <- apply(aP[,c("g1_exp_IMR90", "g2_exp_IMR90")], 1, mean)
aP[,"expBin"] <- factor(breaksExp[.bincode(apply(aP[,c("g1_exp_IMR90", "g2_exp_IMR90")], 1, mean)+1, breaksExp)], labels=c("no", "low", "high"))
aP[,"eh"] <- aP[,"commonEnhancer"] > 0

# save bright allDF data set with column for each source:
write.table(aP, file=paste0(outPrefix, ".allPairs_broad.csv"),
    sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)


# DIMENSION:
# - paralog paris cis       (1560)
# - sampling replicates     (1+10)
# - sample types            (4) dsit/dist+eh/dist+eh+strand/dist+eh+len
# - TAD source              (10) 
# - Expression              (4)
# - other species           (2)
# N=(1560 + 1560*10*4) * 10 * 4 * 2 = 5116800

# remove all TAD, expression, and species related column names
sourceSpecificCols <- c(
    names(allTADs),
    paste0(names(allTADs), "_subTAD"),
    paste0(names(expDFlist), "_expCor"),
    as.vector(sapply(names(speciesSeqInfo), grep, cC, value=TRUE))
)

aPflat <- aP[,setdiff(names(aP), sourceSpecificCols)]

# define number fo different sources 
nTAD <- length(allTADs)
nExp <- length(expDFlist)
nSpecies  <- length(orgStr2Name) 

# define columns for species specific information
speciesCols <- str_split_fixed(grep(names(orgStr2Name)[1], cC, value=TRUE), '_', 2)[,2]


# create data.frame with multiple rows per gene pair for the different data sources
allDF <- cbind(
    aPflat[rep(1:nrow(aP), nTAD*nExp*nSpecies),],
    data.frame(
        tadSource = factor(rep(rep(names(allTADs), each=nrow(aP)), nExp*nSpecies), levels=names(allTADs)),
        inTAD = factor(rep( unlist(aP[,names(allTADs)]), nExp*nSpecies), levels=c(TRUE, FALSE), labels=c("same TAD", "not same TAD")),
        subTAD = rep( unlist(aP[,paste0(names(allTADs), "_subTAD")]), nExp*nSpecies),
        expSource = factor(rep(rep(names(expDFlist), each=nrow(aP)*nTAD), nSpecies), levels=names(expDFlist)),
        expCor = rep( unlist(aP[,rep(paste0(names(expDFlist), "_expCor"), each=nTAD)]), nSpecies),
        expCorBin=factor(breaksExpCor[.bincode(unlist(aP[,rep(paste0(names(expDFlist), "_expCor"), each=nTAD)]), breaksExpCor)], labels=c("high neg", "low neg", "low pos", "high pos")),
        speciesFactor = factor(rep(orgStr2Name, each=nrow(aP)*nTAD*nExp), levels=orgStr2Name)
    ),
    do.call("rbind", lapply(names(orgStr2Name), function(spec){
        spDF <- aP[,paste0(spec, "_", speciesCols)]
        names(spDF) <- paste0("ortholog_", speciesCols)
        return(spDF[rep(1:nrow(aP), nTAD*nExp),])
    }))
)

# add specific group column for being in TAD and/or distal
allDF[,"inTADclose"] <- as.character(allDF$distTadBin)
allDF[allDF$distTadBin == "10-1000kb" ,"inTADclose"] <- as.character(allDF[allDF$distTadBin == "10-1000kb" ,"inTAD"])
allDF$inTADclose <- factor(allDF$inTADclose, levels=c("<10kb", "same TAD", "not same TAD", ">1000kb"))


# add three category sub TAD structure:
sub3TAD <- as.character(allDF[,"subTAD"])
sub3TAD[sub3TAD == "no TAD" | sub3TAD == "diff TAD"] <- "not same TAD"
sub3TAD[sub3TAD == "diff sub TAD"] <- "same TAD"
allDF[,"sub3TAD"] <- factor(sub3TAD, levels=c("not same TAD", "same TAD", "same sub TAD"))

# make ortholog_dist absolut and in kb
allDF[,"ortholog_dist"] <- abs(allDF[,"ortholog_dist"]) / 10^3
allDF[,"ortholog_TAD"] <- factor(allDF[,"ortholog_TAD"], c(TRUE, FALSE), c("same TAD", "not same TAD"))

# add NoZero column for species Hi-C data
allDF <- addNoZero(allDF, cols=c("ortholog_HiC", "ortholog_HiCnorm"))

# save as table
write.table(allDF, file=paste0(outPrefix, ".allDF.csv"),
    sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

# save both data sets as image
save(aP, allDF, file=paste0(WORKIMAGE_FILE, ".aP_allDF.Rdata"))

message("INFO: Wrote aP and allDF to output image file.")

#=======================================================================
# 6.) Mouse and dog paralog analysis
#=======================================================================
runBasicParalogAnalysis(paste0(outPrefix, ".mouse"), paralogPairsMouse, "mmusculus_paralog_ds", speciesTssGR[["mmusculus"]], speciesGenesGR[["mmusculus"]], speciesTADs[["mmusculus"]], tissueName="Mouse", HiClist=speciesHiC[["mmusculus"]][[1]], HiClistNorm=speciesHiC[["mmusculus"]][[2]])
message("Finish Mouse paralog analysis!")

runBasicParalogAnalysis(paste0(outPrefix, ".dog"), paralogPairsDog, "cfamiliaris_paralog_ds", speciesTssGR[["cfamiliaris"]], speciesGenesGR[["cfamiliaris"]], speciesTADs[["cfamiliaris"]], tissueName="Dog", HiClist=speciesHiC[["cfamiliaris"]][[1]], HiClistNorm=speciesHiC[["cfamiliaris"]][[2]])
message("Finish Dog paralog analysis!")


#=======================================================================
# save workspace image
#=======================================================================
message("INFO: Finished with sampling_and_annotation script.")
#save.image(WORKIMAGE_FILE)
#INFO: FANTOM enhancer map: 65953 of 66942 enhancer-associated genes could be mapped uniquelly to ENSG ID
