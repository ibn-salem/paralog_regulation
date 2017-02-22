########################################################################
#
#  This script implements functionality to annotate and analyse the
#  Co-regulation of gene pairs. It is used by the main script 
#  'paralog_regulation.R'
# 
########################################################################

require(stringr)        # for some string functionality
require(GenomicRanges)
require(rtracklayer)    # for import.bed
require(plyr)           # count() function
require(rPython)        # to execute python code from within R

# load other custom modules
source("R/parseHiC.R")

# load custom python script
python.load("python/maximumWeightMatching.py")


#-----------------------------------------------------------------------
# get subset of gene pairs that are located on the same chromosome
#-----------------------------------------------------------------------
getCisPairs <- function(genePairs, tssGR){

    # get chromosomes of gene pairs
    c1 = as.character(seqnames(tssGR[genePairs[,1]]))
    c2 = as.character(seqnames(tssGR[genePairs[,2]]))

    # subset of gene pairs that are located on the same chromosome
    return(genePairs[c1==c2,])
}


#-----------------------------------------------------------------------
# remove double entries of the form A-B and B-A
#-----------------------------------------------------------------------
uniquePair <- function(genePairs){
    
    # get string of sorted IDs as unique pair ID
    pairID = apply(apply(genePairs[,1:2], 1, sort), 2, paste, collapse="_")
    
    genePairs[!duplicated.random(pairID),]
}


#-----------------------------------------------------------------------
# get only one pair per unique gene
#-----------------------------------------------------------------------
uniquePairPerGene <- function(genePairs){

    seen = c()  # set off seen IDs
    uniqPairs = rbind() # initialize output pairs
    for (i in seq(nrow(genePairs))){
        g1 = as.character(genePairs[i,1])
        g2 = as.character(genePairs[i,2])
        
        # check if one of the pairs was seen before
        if (! ( (g1 %in% seen) | (g2 %in% seen) )){
            uniqPairs = rbind(uniqPairs, genePairs[i,])
        }
        # update seen ID set
        seen = c(seen, g1, g2)
    }
    
    return(uniqPairs)
}

#-----------------------------------------------------------------------
# creates an ID for each gene pair by concatenating both gene names
#-----------------------------------------------------------------------
getPairID <- function(genePairs) paste(genePairs[,1], genePairs[,2], sep="_")

#-----------------------------------------------------------------------
# get only one pair per unique gene by choosing the gene pair with highest 
# sequence similarity
#-----------------------------------------------------------------------
uniquePairPerGeneBySim <- function(genePairs, similarity){
    
    # get maximal weight matching of the graph G induced by the paris of genes
    # with similarity as weight. 
    # This command call the function form inside the python script
    matching = python.call( "getMaxWeightMatchingAsDict", genePairs[,1], genePairs[,2], similarity)

    # convert the matching to a data frame of gene pairs
    uniqPairs <- data.frame(names(matching), matching, stringsAsFactors=FALSE)
    
    # for all unique pairs get indices in the input set of pairs
    orgIDs = match(getPairID(uniqPairs), getPairID(genePairs))
    
    # add annotation columns from original data.frame
    uniqPairs = cbind(uniqPairs, genePairs[orgIDs,3:ncol(genePairs)])
    
    # set colom names to those of the input data.frame and delete column names
    names(uniqPairs) <- names(genePairs)
    row.names(uniqPairs) = NULL
    
    return(uniqPairs)
}

#-----------------------------------------------------------------------
# filter gene pairs by negative pair set
#-----------------------------------------------------------------------
containsGenePairs <- function(genePairs, negPairs, gPidx=FALSE, nPidx=FALSE, gr=NULL){
    
    # take either the index directly or get the index from the GRange object
    if (gPidx){
        gP1 = genePairs[,1]
        gP2 = genePairs[,2]
    }else{
        gP1 <- match(genePairs[,1], names(gr))
        gP2 <- match(genePairs[,2], names(gr))
    }
    if (nPidx){
        nP1 = negPairs[,1]
        nP2 = negPairs[,2]
    }else{
        nP1 <- match(negPairs[,1], names(gr))
        nP2 <- match(negPairs[,2], names(gr))
    }
    
    # sort pairs according to index
    gPmin <- apply(cbind(gP1, gP2), 1, min)
    gPmax <- apply(cbind(gP1, gP2), 1, max)

    nPmin <- apply(cbind(nP1, nP2), 1, min)
    nPmax <- apply(cbind(nP1, nP2), 1, max)
    
    # combine id of first and second to get unique ID per pair
    gPid <- paste(gPmin, gPmax, sep="|")
    nPid <- paste(nPmin, nPmax, sep="|")

    inNeg <- gPid %in% nPid
}

#-----------------------------------------------------------------------
# test if gene pairs are non-overlapping
#-----------------------------------------------------------------------
nonOverlappingGenePairs <- function(genePairs, genesGR, useIDs=FALSE){
    
    genesHitDF <- as.data.frame(findOverlaps(genesGR, genesGR))
    ovl <- containsGenePairs(genePairs, negPairs=genesHitDF, gPidx=useIDs, nPidx=TRUE, gr=genesGR)
    return( ! ovl )
    
}

#-----------------------------------------------------------------------
# Add HGNC symbols to genePairs
#-----------------------------------------------------------------------
addHGNC <- function(genePairs, tssGR){
    genePairs[,"HGNC_g1"] = tssGR[genePairs[,1]]$hgnc_symbol
    genePairs[,"HGNC_g2"] = tssGR[genePairs[,2]]$hgnc_symbol
    return(genePairs)
}

#-----------------------------------------------------------------------
# Add same starnd information
#-----------------------------------------------------------------------
addSameStrand <- function(genePairs, tssGR){
    
    # check if both genes have strand information
    s1 = as.character(strand(tssGR[genePairs[,1]]))
    s2 = as.character(strand(tssGR[genePairs[,2]]))
    hasStrandInfo = s1 != "*" & s2 != "*" 
    #check if they are equal
    sameStrand = s1==s2
    
    genePairs[,"sameStrand"] = as.vector(ifelse(hasStrandInfo, sameStrand, NA))
    
    return(genePairs)
}

#-----------------------------------------------------------------------
# Add house keeping genes annotation
#-----------------------------------------------------------------------
addInGeneSet <- function(genePairs, geneSet, colName){

    g1 <- genePairs[,1] %in% geneSet
    g2 <- genePairs[,2] %in% geneSet
    
    genePairs[,colName] <- rep(NA, nrow(genePairs))
    genePairs[,colName][g1 | g2] <- "one"
    genePairs[,colName][g1 & g2] <- "both"
    genePairs[,colName][!g1 & !g2] <- "none"
    
    genePairs[,colName] <- factor(genePairs[,colName], c("none", "one", "both"))

    return(genePairs)

}

#-----------------------------------------------------------------------
# for a pair of genes, count number of common enhancers
#-----------------------------------------------------------------------
getCommonEnhancers <- function(twoGenes, gene2ehID){
    if(length(twoGenes) < 2){ return(0) }
    length(intersect(gene2ehID[[twoGenes[1]]], gene2ehID[[twoGenes[2]]]))
}

#-----------------------------------------------------------------------
# for a set of gene pairs (in ENSG) add, the number of shared enhancers
#-----------------------------------------------------------------------
addCommonEnhancer <- function(genePair, gene2ehID){
    # add new column with number of common enhancers
    genePair[,"commonEnhancer"] = apply(genePair[,1:2], 1, getCommonEnhancers, gene2ehID)
    return(genePair)
}

#-----------------------------------------------------------------------
# add the number of shared enhancers that are upstrem, between, or downstream
# of gene pairs with same strand orientation
#-----------------------------------------------------------------------
addRelativeEnhancerPosition <- function(genePairs, tssGR, gene2ehID, ehGR, colNames=c("eh_up", "eh_cent", "eh_down")){
    
    commonEhID <- apply(genePairs[,1:2], 1, function(gP){
        intersect(gene2ehID[[gP[1]]], gene2ehID[[gP[2]]])
        })
    
    # reduce enhancer coordinates to center
    ehCenter = start(resize(ehGR, 1, fix="center"))
        
    # check if both genes have strand information
    strand1 = as.vector(strand(tssGR[genePairs[,1]]))
    strand2 = as.vector(strand(tssGR[genePairs[,2]]))
    commonStrand = ifelse(strand1==strand2, strand1, NA)
    
    allCounts <- sapply(1:nrow(genePairs), function(i){

        gP <- genePairs[i,1:2]
        st <- commonStrand[i]
        ehIDs <- commonEhID[[i]]
        
        # get start coordinate of genes
        s1 <- start(tssGR[as.character(gP[1])])
        s2 <- start(tssGR[as.character(gP[2])])
        
        # if no common strand return NA counts
        if (is.na(st)){
            ucdCounts <- rep(NA, 3)
        }else{
            # if no enhancer is linked to both return 0 for all possible combinations
            if (length(ehIDs) == 0){
                ucdCounts <- c(0,0,0)
            }else{
                # iterate over enhancers and count "upstream", "center", and "downstream"
                cEh <- ehCenter[ehIDs]
                ucdCounts <- c(
                    sum((st=="+" & cEh <= s1 & cEh <= s2) | (st=="-" & cEh > s1 & cEh > s2)),
                    sum( (cEh > s1 & cEh <= s2) | (cEh <= s1 & cEh > s2) ),
                    sum( (st=="+" & cEh > s1 & cEh > s2) | (st=="-" & cEh <= s1 & cEh <= s2) )
                    )

                # double check that all cases are covered
                stopifnot(sum(ucdCounts) == length(ehIDs))
            }
        }
        return(ucdCounts)
    })
    
    # add counts to gene pair df
    genePairs[,colNames] <- t(allCounts)
    
    return(genePairs)
    
}

#-----------------------------------------------------------------------
# adds information of common compartment to gene pair data.frame
#-----------------------------------------------------------------------
addCommonCompartment <- function(genePair, tssGR, compGR, subCompGR){

    # get compartment for each gene in the pairs
    g1Hit <- findOverlaps(tssGR[genePair[,1]], subCompGR)
    g2Hit <- findOverlaps(tssGR[genePair[,2]], subCompGR)
    g1Comp <- rep(NA, nrow(genePair))
    g2Comp <- rep(NA, nrow(genePair))
    g1Comp[queryHits(g1Hit)] <- subCompGR[subjectHits(g1Hit)]$comp
    g2Comp[queryHits(g2Hit)] <- subCompGR[subjectHits(g2Hit)]$comp
    
    # combine compartment annotation 
    compComb <-paste(g1Comp, g2Comp, sep="/")
    compComb[compComb == "B/A"] <- "A/B"
    compComb[is.na(g1Comp) | is.na(g2Comp)] <- NA
    
    # set pairs within the same compartment region to NA
    pairGR <- getPairAsGR(genePair, tssGR)
    sameCompReg <- countOverlaps(pairGR, compGR, type="within") >= 1
    
    # get boolean vector for common compartment annotation but not same compartment region
    commonComp <- compComb %in% c("A/A","B/B") & !sameCompReg
        
    #-------------------------
    # for sub-compartment
    #-------------------------

    # get sub-compartment for each gene in the pairs
    g1SubComp <- rep(NA, nrow(genePair))
    g2SubComp <- rep(NA, nrow(genePair))
    g1SubComp[queryHits(g1Hit)] <- subCompGR[subjectHits(g1Hit)]$subcomp
    g2SubComp[queryHits(g2Hit)] <- subCompGR[subjectHits(g2Hit)]$subcomp

    # combine sub-compartment annotation 
    compSubComb <-paste(g1SubComp, g2SubComp, sep="/")
    compSubComb[is.na(g1SubComp) | is.na(g2SubComp)] <- NA

    # set pairs within the same sub-compartment region to NA
    sameSubCompReg <- countOverlaps(pairGR, subCompGR, type="within") >= 1
    
    # get boolean vector for common sub-compartment annotation but not same sub-compartment region
    commonSubComp <- g1SubComp == g2SubComp & !is.na(g1SubComp) & !is.na(g2SubComp) & !sameSubCompReg
    
    # add all annotations to gene pair df:
    genePair[,"comp_combination"] <- compComb
    genePair[,"same_comp_region"] <- sameCompReg
    genePair[,"common_comp"] <- commonComp
    genePair[,"subcomp_combination"] <- compSubComb
    genePair[,"same_subcomp_region"] <- sameSubCompReg
    genePair[,"common_subcomp"] <- commonSubComp

    return(genePair)
}



#-----------------------------------------------------------------------
# distribution of shared enhancers among pairs of genes
#-----------------------------------------------------------------------
percentCounts <- function(counts){
    df = count(counts)
    names(df) = c("count", "freq")
    df$freq = 100 * df$freq / length(counts)
    return(df)
}

#-----------------------------------------------------------------------
# add linear distance between genes (by assuming them on the same chromosome)
#-----------------------------------------------------------------------
addPairDist <- function(genePairs, tssGR){
    # get chromosomes of gene pairs
    sameChrom = as.character(seqnames(tssGR[genePairs[,1]])) == as.character(seqnames(tssGR[genePairs[,2]]))
    s1 = start(tssGR[genePairs[,1]])
    s2 = start(tssGR[genePairs[,2]])
    # add a new column "dist" to the data.frame
    genePairs[, "dist"] = ifelse(sameChrom, s2-s1, NA)
    return(genePairs)
}

#-----------------------------------------------------------------------
# make GRange object of gene pairs with TSS on same chromosome
#-----------------------------------------------------------------------
getPairAsGR <- function(genePairs, tssGR){
    # get chromosomes of gene pairs
    chrom = seqnames(tssGR[genePairs[,1]])
    
    s1 = start(tssGR[genePairs[,1]])
    s2 = start(tssGR[genePairs[,2]])
    up = apply(cbind(s1, s2), 1, min)
    down = apply(cbind(s1, s2), 1, max)
    GR = GRanges(chrom, IRanges(up, down))
    # add gene IDs and other annotations
    mcols(GR) = genePairs
    return(GR)
}


#-----------------------------------------------------------------------
# add column to indicate that query lies within at least one subject object
#-----------------------------------------------------------------------
addWithinSubject <- function(query, subject, colName="inRegion"){
    mcols(query)[, colName] = countOverlaps(query, subject, type="within") >= 1
    return(query)
}

#-----------------------------------------------------------------------
# add column to indicate that two query region overlap the same subset of subject regions.
# If both query regions do not overlap any subject, FALSE is returned.
#-----------------------------------------------------------------------
addSubTADmode <- function(genePairs, tadGR, tssGR, colName="subTAD"){
    
    # compute overlap of all genes with TADs
    hitsDF <- as.data.frame(findOverlaps(tssGR, tadGR, type="within"))
    
    # get indices of genes in tssGR
    idx1 <- match(genePairs[,1], names(tssGR))
    idx2 <- match(genePairs[,2], names(tssGR))
    
    # get for each gene the set of TADs that overlap
    set1 <- lapply(idx1, function(i) hitsDF[hitsDF[,1]==i,2])
    set2 <- lapply(idx2, function(i) hitsDF[hitsDF[,1]==i,2])
    
    # test if sets are equal for each gene in pair
    commonSubset <- mapply(setequal, set1, set2)
    
    # compute overlap with any common TAD
    genePairGR <- getPairAsGR(genePairs, tssGR)
    commonTAD <- countOverlaps(genePairGR, tadGR, type="within") >= 1

    # get combinations
    subTADmode <- NA
    subTADmode[!commonTAD & commonSubset] <- "no TAD"
    subTADmode[!commonTAD & !commonSubset] <- "diff TAD"
    subTADmode[commonTAD & !commonSubset] <- "diff sub TAD"
    subTADmode[commonTAD & commonSubset] <- "same sub TAD"
    subTADmode <- factor(subTADmode, levels=c("no TAD", "diff TAD", "diff sub TAD", "same sub TAD"))

    # add column to gene pairs and return 
    genePairs[, colName] <- subTADmode
    return(genePairs)
}

#-----------------------------------------------------------------------
# returns the percent of TRUE elements in input vector v
#-----------------------------------------------------------------------
percentTrue <- function(v){
    sum(v, na.rm=TRUE)/length(v) * 100
}

#-----------------------------------------------------------------------
# Inter-chromosomal gene pairs counts matrix
#-----------------------------------------------------------------------
interChromPairMatrix <- function(genePairs, tssGR, symmetric=FALSE){
    
    # get vector of all unique chromosome names
    chroms = seqnames(seqinfo(tssGR))
    
    # initialize matrix with zero counts
    n = length(chroms)
    mat = matrix(rep(0, n*n), n, dimnames=list(chroms, chroms))
    
    # get chromsome names of gene pairs
    c1 = seqnames(tssGR[genePairs[,1]])
    c2 = seqnames(tssGR[genePairs[,2]])
    
    # count pairwise occurrences
    counts = count(data.frame(c1, c2, stringsAsFactors=FALSE))
    
    # iterate over all found pairs and increase counter
    for (i in 1:nrow(counts)){
        mat[counts[i,1], counts[i,2]] = mat[counts[i,1], counts[i,2]] + counts[i,3]
        
        # if option symmetric is FALSE, count pair as cA-cB and cB-cA
        if (!symmetric){
            mat[counts[i,2], counts[i,1]] =  mat[counts[i,2], counts[i,1]] + counts[i,3]
        }
    }
    
    # make single letter chromosome names
    rownames(mat) = gsub("chr", "", rownames(mat))
    colnames(mat) = gsub("chr", "", colnames(mat))
    
    return(mat)
}


#-----------------------------------------------------------------------
# This function returns a logical vector, the elements of which are FALSE, unless there are duplicated values in x, in which case all but one elements are TRUE (for each set of duplicates). The only difference between this function and the duplicated() function is that rather than always returning FALSE for the first instance of a duplicated value, the choice of instance is random.
# Source: https://amywhiteheadresearch.wordpress.com/2013/01/22/randomly-deleting-duplicate-rows-from-a-dataframe-2/
#-----------------------------------------------------------------------
duplicated.random = function(x, incomparables = FALSE, ...) { 
     if ( is.vector(x) ) { 
         permutation = sample(length(x)) 
         x.perm      = x[permutation] 
         result.perm = duplicated(x.perm, incomparables, ...) 
         result      = result.perm[order(permutation)] 
         return(result) 
     } 
     else if ( is.matrix(x) ) { 
         permutation = sample(nrow(x)) 
         x.perm      = x[permutation,] 
         result.perm = duplicated(x.perm, incomparables, ...) 
         result      = result.perm[order(permutation)] 
         return(result) 
     } 
     else{ 
         stop(paste("duplicated.random() only supports vectors", 
                "matrices for now.")) 
     } 
} 

#-----------------------------------------------------------------------
# Adds Hi-C contact frequencies to a gene pair data set
#-----------------------------------------------------------------------
addHiCfreq <- function(genePairs, tssGR, HiClist, label="HiC", ...){
    
    xRange <- tssGR[as.character(genePairs[,1])]
    yRange <- tssGR[as.character(genePairs[,2])]
    
    if(all(width(xRange) == 1 & width(yRange) == 1)){

        genePairs[,label] = getInteractionsPoint(
                xRange, 
                yRange, 
                HiClist,
                ...)
    }else{
        genePairs[,label] = getInteractionsMulti(
                xRange, 
                yRange, 
                HiClist,
                ...)
    }
    return(genePairs)
}


#-----------------------------------------------------------------------
# add Hi-C observed / expected ratio 
#-----------------------------------------------------------------------
addHiCobsExp <- function(genePairs, tssGR, expectedHiCList, resolution, HiClabel="HiCfreq", label="HiCobs/exp"){
    
    # assume pairs on same chromosome
    stopifnot(all(seqnames(tssGR[genePairs[,1]]) == seqnames(tssGR[genePairs[,1]])))
    chroms = seqnames(tssGR[genePairs[,1]])

    # iterate over all unique chromosomes (in parallel)
    for (chr in as.character(unique(chroms))){
    
        message(paste("INFO: Query expected contacts for chromosome:", chr))

        # get indexes of input ranges on that chrom
        onChrom = which(chroms == chr)
        
        if (chr %in% names(expectedHiCList)){

            # get vector with expected counts for this chrom
            e <- expectedHiCList[[chr]]
            
            # calculate the appropriate distance bin
            d <- abs(genePairs[onChrom, "dist"])
            distIDX <- d %/% resolution + 1
    
            # get the expected count for each gene pair
            expectedContacts <- e[distIDX]
        }else{
            expectedContacts <- rep(NA, length(onChrom))
        }
        # add obs/exp value to gene pairs
        genePairs[onChrom, label] <- genePairs[onChrom, HiClabel] / expectedContacts
    }
    
    return(genePairs)

}


#-----------------------------------------------------------------------
# Returns a unique gene ID for a pair of strings that is independent of pair order (e.g. sorted)
#-----------------------------------------------------------------------
getPairIDsorted <- function(gP){
    mapply(function(g1, g2) 
        paste(
            sort(c(as.character(g1), as.character(g2))), 
        collapse="_"), gP[,1], gP[,2])
}


#-----------------------------------------------------------------------
# returns the percentage of gene pairs in gPa that are found in gPb
#-----------------------------------------------------------------------
percentIncluded <- function(gPa, gPb){
    a <- getPairIDsorted(gPa)
    b <- getPairIDsorted(gPb)
    return( 100 * sum(a %in% b) / length(a))
}

#-----------------------------------------------------------------------
# Adds HIPPIE PPI interaction score
#-----------------------------------------------------------------------
addHIPPIE <- function(gP, hippie, colName="HIPPIE"){

    gPID <- getPairIDsorted(gP)
    
    # get index of gene pair in hippie by matching symbols
    score <- hippie[match(gPID, hippie$ID), "score"]
    
    gP[,colName] <- score
    return(gP)
}

#-----------------------------------------------------------------------
# Adds the expression values for both genes
#-----------------------------------------------------------------------
addPairExp <- function(gP, expDF, expCol, label="exp"){

    g1exp <- expDF[gP[,1], expCol]
    g2exp <- expDF[gP[,2], expCol]

    gP[,paste0("g1_", label)] <- g1exp
    gP[,paste0("g2_", label)] <- g2exp
    
    return(gP)
}

#-----------------------------------------------------------------------
# returns the Pearson correlation coefficient for expression of two input genes
#-----------------------------------------------------------------------
getCor <- function(gP, expDF){
    
    # correct intput to matrix (in case of only two element vector)
    gP = matrix(gP, ncol=2)
    
    # get correlation values of all cells/conditions
    # this will make a vector of NA's if the gene is not contained in the expression data set
    # furthermore, cbind(c(.)) guarantees that cor() will deal with column-vectors
    x = t(as.vector(expDF[gP[,1],]))
    y = t(as.vector(expDF[gP[,2],]))
    
    # return pearson correlation coefficient
    cor(x, y, method="pearson")

}


#-----------------------------------------------------------------------
# adds Pearson correlation coefficient for all gene pairs 
#-----------------------------------------------------------------------
addCor <- function(gP, expDF, colName="expCor"){
    pairsAsChars = sapply(gP[,1:2], as.character)
    gP[,colName] = apply(pairsAsChars, 1, getCor, expDF=expDF)
    return(gP)
}

#-----------------------------------------------------------------------
# Map genes to one2one orthologs in an other species
#-----------------------------------------------------------------------
getOrthologs <- function(genePairs, orthologsAll, orgStr, tssGR){
    
    # filter for one2one orthologs
    isOne2one = orthologsAll[,paste0(orgStr, "_homolog_orthology_type")] == "ortholog_one2one"
    isInTssGR = orthologsAll[,paste0(orgStr, "_homolog_ensembl_gene")] %in% names(tssGR)
    
    orthologs = orthologsAll[isOne2one & isInTssGR, c("ensembl_gene_id", paste0(orgStr, "_homolog_ensembl_gene"))]
    
    # remove duplicates (arising from one2many orthologsi in other species)
    orthologs = orthologs[!duplicated(orthologs),]
    rownames(orthologs) = orthologs[,1]
    
    # map each gene from the pair to its one2one ortholog
    o1 = orthologs[as.character(genePairs[,1]),2]
    o2 = orthologs[as.character(genePairs[,2]),2]
    
    orthPairs = data.frame(g1=o1, g2=o2, stringsAsFactors=FALSE)
    
    return(orthPairs)
}

#-----------------------------------------------------------------------
# return a boolean string wehter the two genes have a common ortholog
#-----------------------------------------------------------------------
commonOrthologs <- function(genePairs, orthologsAll, orgStr, types=c("ortholog_one2many", "ortholog_one2one")){
    
    # filter for allowed orthology types:
    orthologs <- orthologsAll[orthologsAll[,paste0(orgStr, "_homolog_orthology_type")] %in% types,]
        
    oList1 <- lapply(genePairs[,1], function(g) orthologs[orthologs[,1] == g, 2]) 
    oList2 <- lapply(genePairs[,2], function(g) orthologs[orthologs[,1] == g, 2]) 
    
    commonOrtholog <- mapply(intersect,oList1, oList2)
    
    commonOrthNumber <- sapply(commonOrtholog, length)
    
    return(commonOrthNumber > 0)
}

#-----------------------------------------------------------------------
# checks if both genes in a set of gene paris are not NA
#-----------------------------------------------------------------------
pairNotNA <- function(genePairs){
    return( !is.na(genePairs[,1]) & !is.na(genePairs[,2]) )
    
}

#-----------------------------------------------------------------------
# add information of the location of one-two-one orthologs of the gene paris
#-----------------------------------------------------------------------
addOrthologAnnotation <- function(genePairs, orthologsAll, orgStr, tssGR, TAD=NULL, HiClist=NULL, HiClistNorm=NULL, inParallel=TRUE){

    # get orthologs pairs
    orthoPairs = getOrthologs(genePairs, orthologsAll, orgStr, tssGR)
    
    # add bool flag if pair has for both gens one-to-one orthologs
    hasOne2one = pairNotNA(orthoPairs)
    genePairs[,paste0(orgStr, "_one2one")] = hasOne2one
    
    # check if orthologs are on the same chrom
    sameChrom = as.vector(seqnames(tssGR[orthoPairs[hasOne2one,1]]) == seqnames(tssGR[orthoPairs[hasOne2one,2]]))
    
    genePairs[hasOne2one, paste0(orgStr, "_sameChrom")] =  sameChrom

    # add linear distance between orthologs of pairs
    orthologsDist = addPairDist(orthoPairs[hasOne2one,], tssGR)[,"dist"]
    genePairs[hasOne2one,  paste0(orgStr, "_dist")] =  orthologsDist
    
    # add co-occurances in same TAD
    if (!is.null(TAD)){
        orthoPairsGR = getPairAsGR(orthoPairs[hasOne2one,], tssGR)
        tadColName = paste0(orgStr, "_TAD")
        orthoPairsGR = addWithinSubject(orthoPairsGR, TAD, tadColName)
        genePairs[hasOne2one, tadColName] = mcols(orthoPairsGR)[, tadColName]
    }
    
    # add Hi-C counts
    if (!is.null(HiClist) && !is.null(HiClistNorm)){    
        subOnSameChrom = which(genePairs[,paste0(orgStr, "_sameChrom")])
        
        genePairs[subOnSameChrom, paste0(orgStr, "_HiC")] = addHiCfreq(orthoPairs[subOnSameChrom,], tssGR, HiClist, label="rawHiC", inParallel=inParallel)$rawHiC
        genePairs[subOnSameChrom, paste0(orgStr, "_HiCnorm")] = addHiCfreq(orthoPairs[subOnSameChrom,], tssGR, HiClistNorm, label="HiCnorm", inParallel=inParallel)$HiCnorm
    }
    
    return(genePairs)
}

#-----------------------------------------------------------------------
# get pairwise information from matrix. Matrix is assumed to have gene names (matching the first two columns of query data.frame) as dimensions names
#-----------------------------------------------------------------------
getPairwiseMatrixScoreByName <- function(genePairs, M, replaceZeroByNA=FALSE){
    
    idx1 = match(genePairs[,1], dimnames(M)[[1]])
    idx2 = match(genePairs[,2], dimnames(M)[[2]])
    scores = M[cbind(idx1, idx2)]
    
    if (replaceZeroByNA){
        # For capture C data from Mifsud et al. 2015:
        # Due to sparse matrix data structure non available pairs will get 0 counts
        # Since no 0 count pair is in the original data, we can replace all 0 with NA
        scores[scores==0] <- NA
    }
    return(scores)
}

#-----------------------------------------------------------------------
# add flags for commonOrthologs and NotOneToOne
#-----------------------------------------------------------------------
addAgeFlags <- function(genePairs, orthologsSpecies, orgStr){    
    genePairs[, paste0(orgStr, "_commonOrtholg")] <- commonOrthologs(genePairs, orthologsSpeciesList[[orgStr]], orgStr)
    genePairs[, paste0(orgStr, "_NotOne2one")] <- ! genePairs[, paste0(orgStr, "_one2one")]
    return(genePairs)
}


#-----------------------------------------------------------------------
# help function to add additional column with zeros as NA to data frame
#-----------------------------------------------------------------------
addNoZero <- function(df, cols=c("HiCRaw", "HiC", "HiCobsExp", "captureC_raw", "captureC_ObsExp")){
    for (column in cols){
        df[,paste0(column, "NoZero")] <- df[,column]
        df[,paste0(column, "NoZero")][df[,paste0(column, "NoZero")] == 0] = NA
    }
    return(df)
}

