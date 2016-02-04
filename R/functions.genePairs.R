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
#~ require(entropy)        # for function 'mi.plugin()' to calculate mutual information
require(rPython)        # to execute python code from within R
require(minerva)        # to calculate maximal information coefficient (MIC) Reshef et al 2011


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
# returns true if only unique pairs (no A-B, B-A dups) are contained in genePairs
#-----------------------------------------------------------------------
hasDupPairs <- function(genePairs){
    
    pairID = apply(apply(genePairs[,1:2], 1, sort), 2, paste, collapse="_")
    length(pairID) != length(unique(pairID))
    
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
# get only one pair per paralog gene group by choosing the pair with highest 
# sequence similarity
#-----------------------------------------------------------------------
uniquePairPerGroupBySim <- function(genePairs, similarity){

    
    # load the igraph package
    require(igraph)
    
    # create the graph from gene paris
    g = graph.data.frame(genePairs, directed=FALSE)
    
    # take similarity as edge weight
    E(g)$weight = similarity
    
    # get connected components
    conComp = clusters(g)
    message(paste("Finished computation of", conComp$no, "connected components"))
    
    maxPairs = sapply(1:conComp$no, function(i) {
    
        
        # get nodes corresponding to this component
        nodes = V(g)[conComp$membership == i]

        # get edges:
        edges = E(g)[inc(nodes)]
        
        ## get index of edge with max edge weight
        maxEdgeIdx = as.numeric(edges)[which.max(edges$weight)][1]
                
        #message(paste("Work on cluster ", i))
        
        # get index of edges with max weight
        return(maxEdgeIdx)
    })
    
    return(genePairs[maxPairs,])
        
}
#maxSimPairs = uniquePairPerGroupBySim(paralogPairs[1:1000,], paralogPairs$hsapiens_paralog_perc_id[1:1000])

#-----------------------------------------------------------------------
# Filters for gene paris that occures only as tow-pairs (not in triplets, quads...)
#-----------------------------------------------------------------------
onlyTwoPairs <- function(genePairs, includeDuplicates=FALSE){

    allGeneNames = c(genePairs[,1], genePairs[,2])
    
    # count the occurrences of each gene
    geneCounts = count(allGeneNames)
    rownames(geneCounts) = geneCounts[,1]
    
    # for each pair test if gene are only in two paris
    if(includeDuplicates){
        isTwoPair = geneCounts[genePairs[,1], 2] == 2 & geneCounts[genePairs[,2], 2] == 2
    }else{
        isTwoPair = geneCounts[genePairs[,1], 2] == 1 & geneCounts[genePairs[,2], 2]== 1
    }

    return(genePairs[isTwoPair,])
}


#-----------------------------------------------------------------------
# get all adjacent gene pairs as gene pair data.frame
#-----------------------------------------------------------------------
getAdjacentPairs <- function(tssGR){
    
    # sort by ignoring strand
    sortedTSS = sort(tssGR, ignore.strand=TRUE)
    # make strand artificially to '+' because precede is calculated with respect to transcription direction
    strand(sortedTSS) = "+"
    
    # get the next tss for each gene along the chromsome
    # using the precede() function from GenomicRanges package
    nextGene = precede(sortedTSS)

    gP = data.frame(
        g1 = as.character(id(sortedTSS[!is.na(nextGene)])), 
        g2 = as.character(id(sortedTSS[nextGene[!is.na(nextGene)]]))
    )

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
# function to add the information if orthologs are directly adjacent to each other
#-----------------------------------------------------------------------
orthologsAdjacent <- function(gP, speciesLabel, orthologsSpeciesList, tssGRspecies){
    
    # sort by ignoring strand
    sortedTSS = sort(tssGRspecies, ignore.strand=TRUE)
    # make strand artificially to '+' because precede is calculated with respect to transcription direction
    strand(sortedTSS) = "+"

    g1 = as.character(gP[,1])
    g2 = as.character(gP[,2])
    orthoMap = orthologsSpeciesList[[speciesLabel]]
    
    g1_orthoIDX = match(g1, orthoMap[,1])
    g2_orthoIDX = match(g2, orthoMap[,1])
    
    g1_type = orthoMap[g1_orthoIDX, paste0(speciesLabel, "_homolog_orthology_type")]
    g2_type = orthoMap[g2_orthoIDX, paste0(speciesLabel, "_homolog_orthology_type")]
    
    g1_gene = orthoMap[g1_orthoIDX, paste0(speciesLabel, "_homolog_ensembl_gene")]
    g2_gene = orthoMap[g2_orthoIDX, paste0(speciesLabel, "_homolog_ensembl_gene")]

    # check that the homology type of both genes is 'one2one'
    boghOne2one = g1_type == "ortholog_one2one" & g2_type == "ortholog_one2one"
    
    # check that both orhtologs are in the set of TSS of this species
    bothInTSS = g1_gene %in% names(sortedTSS) & g2_gene %in% names(sortedTSS)

    boghOne2one = boghOne2one & !is.na(boghOne2one) & bothInTSS


    preIDX = precede(sortedTSS[g1_gene[boghOne2one]], sortedTSS)
    pre = rep(NA, length(preIDX))
    pre[!is.na(preIDX)] = names(sortedTSS[preIDX[!is.na(preIDX)]])

    folIDX = follow(sortedTSS[g1_gene[boghOne2one]], sortedTSS)
    fol = rep(NA, length(folIDX))
    fol[!is.na(folIDX)] = names(sortedTSS[folIDX[!is.na(folIDX)]])
    
    isAdjacent = (pre == g2_gene[boghOne2one]) | (fol == g2_gene[boghOne2one])

    # get values for all input gene pairs
    orthologsAdjacent = rep(NA, nrow(gP))
    orthologsAdjacent[boghOne2one] = isAdjacent
    
    return(orthologsAdjacent)
}


#-----------------------------------------------------------------------
# get mapping of gene names to enhancer IDs
#-----------------------------------------------------------------------
getGenetoEhIDmapping <- function(geneNames, enhancerIDs){
    
    stopifnot(length(geneNames) == length(enhancerIDs))
    
    # initilize empty list
    l = list()
    
    # iterate over both vectors
    for (i in seq(length(geneNames))){
        g = as.character(geneNames[i])
        e = enhancerIDs[i]
        # append e to list entry with key g
        l[[g]] = c(l[[g]], e)
    }
    return(l)
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
#~     genePair[,"commonEnhancer"] = apply(genePair[,c("HGNC_g1", "HGNC_g2")], 1, getCommonEnhancers, gene2ehID)
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
# for a set of gene pairs (in ENSG) add binary vector indicating same annotaion for each annotation column
#-----------------------------------------------------------------------
addCommonAnnotation <- function(genePair, tssGR, columns=c("comp", "subcomp")){
    
    for (colStr in columns){
        genePair[,paste0("common_", colStr)] = mcols(tssGR[genePair[,1]])[, colStr] == mcols(tssGR[genePair[,2]])[, colStr]
    }
    return(genePair)
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
# get column to indicate that query lies within at least one subject object
#-----------------------------------------------------------------------
getWithinSubject <- function(query, subject){
    countOverlaps(query, subject, type="within") >= 1
}

#-----------------------------------------------------------------------
# make a gene pairs data frame with all gene paris that overlap with the same query region
#-----------------------------------------------------------------------
getPairsFromOverlap <- function(queryGR, subjectGR, colNames=c("g1", "g2")){
    
    # calculate overlap
    hit = findOverlaps(queryGR, subjectGR)
    
    # convert to list and keep only those grups with at least 2 subjects in it
    hitList = as.list(hit)
    hitList = hitList[sapply(hitList, length) >= 2]
    
    # get all possible pairs per group
    pairsID = do.call("rbind", sapply(sapply(hitList, combn, 2), t))

    # convert indexes to names of subject
    
    pairDF = data.frame(names(subjectGR)[pairsID[,1]], names(subjectGR)[pairsID[,2]])
    names(pairDF) = colNames
    
    return(pairDF)
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
# shared enhancer number matrix
#-----------------------------------------------------------------------
sharedEnhancerMatrix <- function(geneIDs, gene2ehID, geneSymbols=geneIDs){

    n = length(geneIDs)
    mat = matrix(rep(0, n*n), n, dimnames=list(geneSymbols, geneSymbols))
    for (i in 1:n){
    
        for (j in 1:n){
            
            mat[i,j] = getCommonEnhancers(geneIDs[c(i,j)], gene2ehID)
        }
    }
    return(mat)
}

#-----------------------------------------------------------------------
# Pairwise linear distance matrix
#-----------------------------------------------------------------------
pairwiseDistMatrix <- function(gr, geneSymbols=names(gr)){

    distMat <- outer(start(gr), start(gr), function(x,y) abs(x-y))
    sameChromMat <- outer(seqnames(gr), seqnames(gr), function(x,y) as.logical(x==y))
    
    # set pairs on different chromosomes to NA
    distMat[!sameChromMat] <- NA
    dimnames(distMat) <- list(geneSymbols, geneSymbols)
    return(distMat)
    
}


#-----------------------------------------------------------------------
# Pairwise Hi-C contact counts as matrix
#-----------------------------------------------------------------------
pairwiseContacstMatrix <- function(exampleGenes, HiClist){

    n = length(exampleGenes)
    
    # combine query gr
    queryGRa = rep(exampleGenes, each=n)
    queryGRb = rep(exampleGenes, n)
    
    contacts = getInteractionsMulti(queryGRa, queryGRb, HiClist)
    
    mat = matrix(contacts, n, dimnames=list(exampleGenes$hgnc_symbol, exampleGenes$hgnc_symbol))
    
    return(mat)
}

#-----------------------------------------------------------------------
# Pairwise Hi-C contact counts as matrix
#-----------------------------------------------------------------------
pairwiseContacstMatrixSameChrom <- function(exampleTSS, HiClist){
    
    # assume only one chromosome
    stopifnot(length(unique(seqnames(exampleTSS))) == 1)
    stopifnot(all(width(exampleTSS) == 1))
    
    chr = seqnames(exampleTSS)[1]
    map = HiClist[[paste0(chr, chr)]]
 
    # get indexes of bins in Hi-C map overlapping the query regions
    idxX = subjectHits(findOverlaps(exampleTSS, x_intervals(map))) 
    idxY = subjectHits(findOverlaps(exampleTSS, y_intervals(map))) 
    
    # query interaction matrix
    mat = intdata(map)[idxX, idxY]

    return(as.matrix(mat))
}


#-----------------------------------------------------------------------
# This function returns a logical vector, the elements of which are FALSE, unless there are duplicated values in x, in which case all but one elements are TRUE (for each set of duplicates). The only difference between this function and the duplicated() function is that rather than always returning FALSE for the first instance of a duplicated value, the choice of instance is random.
# Source: https://amywhiteheadresearch.wordpress.com/2013/01/22/randomly-deleting-duplicate-rows-from-a-dataframe-2/
#-----------------------------------------------------------------------
duplicated.random = function(x, incomparables = FALSE, ...) 
{ 
     if ( is.vector(x) ) 
     { 
         permutation = sample(length(x)) 
         x.perm      = x[permutation] 
         result.perm = duplicated(x.perm, incomparables, ...) 
         result      = result.perm[order(permutation)] 
         return(result) 
     } 
     else if ( is.matrix(x) ) 
     { 
         permutation = sample(nrow(x)) 
         x.perm      = x[permutation,] 
         result.perm = duplicated(x.perm, incomparables, ...) 
         result      = result.perm[order(permutation)] 
         return(result) 
     } 
     else 
     { 
         stop(paste("duplicated.random() only supports vectors", 
                "matrices for now.")) 
     } 
} 

#-----------------------------------------------------------------------
# Adds Hi-C contact frequencies to a gene pair data set
#-----------------------------------------------------------------------
addHiCfreq <- function(genePair, tssGR, HiClist, label="HiCfreq", ...){
    
    genePair[,label] = getInteractionsMulti(
            tssGR[genePair[,1]], 
            tssGR[genePair[,2]], 
            HiClist,
            ...)
    return(genePair)
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
#getCor(cbind("ENSG00000000003", "ENSG00000000005"), expDF)
#getCor(cbind(c("ENSG00000000003", "ENSG00000000005"), c("ENSG00000000003", "ENSG00000000005")),expDF)
#getCor(cbind(c("ENSG00000000003", "ENSG00X00000005"), c("ENSG00000000003", "ENSG00000000005")),expDF)


#-----------------------------------------------------------------------
# adds Pearson correlation coefficient for all gene pairs 
#-----------------------------------------------------------------------
addCor <- function(gP, expDF, colName="expCor"){
    pairsAsChars = sapply(gP[,1:2], as.character)
    gP[,colName] = apply(pairsAsChars, 1, getCor, expDF=expDF)
#~     gP[,colName] = sapply(1:nrow(gP), function(i) {
#~         cor(
#~             x=t(expDF[gP[i,1],]), 
#~             y=t(expDF[gP[i,2],]), 
#~             method="pearson")
#~         })
    return(gP)
}

#-----------------------------------------------------------------------
# adds maximum information coefficient (MIC) or other maximal information-based nonparametric exploration (MINE) statistics (Reshef et al. Science 2011) 
#-----------------------------------------------------------------------
addMIC <- function(genePairs, expDF, statistic="MIC", colName=statistic, ...){
    
    g1exp = t(expDF[as.character(genePairs[,1]),])
    g2exp = t(expDF[as.character(genePairs[,2]),])
    
    genePairs[,colName] <- sapply(1:nrow(genePairs), function(j) 
        # check if one of the 
        if(all(is.na(g1exp[,j])) | all(is.na(g2exp[,j])) ){
            NA
        }else{
            mine(x=g1exp[,j], y=g2exp[,j], ...)[[statistic]]
        })
    return(genePairs)
}

#-----------------------------------------------------------------------
# makes a huge pdf with many dotplot of pairwise expression comparisons
#-----------------------------------------------------------------------
plotAllExp <- function(genePairs, expDF, tssGR, outFile=NA, ...){

    nRow = round(sqrt(nrow(genePairs))) 
    nCol = round(sqrt(nrow(genePairs))) + 1
    
    stopifnot(nRow * nCol >= nrow(genePairs))
    
    if (!is.na(outFile)){
        pdf(outFile, width=3*nCol, height=3*nRow)
    }
    par(mfrow=c(nRow,nCol), cex=.7, mar=c(4.1,4.1,1.1,2.1))
    
    for (i in 1:nrow(genePairs)){
        
        x = t(expDF[as.character(genePairs[i,1]), ])
        y = t(expDF[as.character(genePairs[i,2]), ])
        
        m = c("cor"=cor(x,y), unlist(mine(x,y)))
        #, collapse=" ")
        xName = paste0(tssGR[genePairs[i,1]]$hgnc_symbol, " (", genePairs[i,1], ")")
        yName = paste0(tssGR[genePairs[i,2]]$hgnc_symbol, " (", genePairs[i,2], ")")

        plot(x,y, 
            xlab=xName, ylab=yName,
            pch=21, bg="blue", ...)
        legend("topright", paste0(names(m), "=", signif(m, 2)), bty="n")
    }

    if (!is.na(outFile)){
        dev.off()
    }
    
}
#plotAllExp(cisPairs[which(mic_r2 >= .4),], expDF)
#plotAllExp(cisPairs[which(mic_r2 >= .4),], expDF, tssGR, "out.test.pdf")
#plotAllExp(cisPairs[which(mic_r2 >= .25 & r<=0),], expDF, tssGR, "out.test.pdf")


#-----------------------------------------------------------------------
# returns the mutual information for the expression of two input genes
#-----------------------------------------------------------------------
getExpMI <- function(gP, expDF){
    
    stopifnot(nrow(gP) == 1)

    # check if both genes are contained in expression data set
    if (gP[1] %in% row.names(expDF) && gP[2] %in% row.names(expDF)) {

        # get expression values of all cells/conditions
        # this will make a vector of NA's if the gene is not contained in the expression data set    
        freqs2d = expDF[c(gP[1], gP[2]) ,]
        
        # if only zeros in expression vecotrs, return NA as mutual information
        if (sum(freqs2d) == 0) {
            return(NA)
        }
        
        # return mutual information between genes in bit units
        mi.plugin(freqs2d, unit="log2")
    }else{
        return(NA)
    }
    
}

#getExpMI(cbind("ENSG00000000003", "ENSG00000000005"), expDF)
#gP = rbind(c("ENSG00000000003", "ENSG00000000005"), c("ENSG00000000003", "ENSG00000000005"))
#gP2 = rbind(c("ENSG00000000003", "ENSG00X00000005"), c("ENSG00000000003", "ENSG00000000005"))
#apply(gP, 1, getExpMI, expDF=expDF)
#apply(gP2, 1, getExpMI, expDF=expDF)


#-----------------------------------------------------------------------
# adds mutual inforamtion of gene expression for all gene pairs 
# TODO: improve runtime by precalculating entropy of x and y and compute
# MI as MI(X,Y) = H(X) + H(Y) - H(X, Y)
#-----------------------------------------------------------------------
addExpMI <- function(gp, expDF, colName="expMI"){
    gp[,colName] = apply(gp[,1:2], 1, getExpMI, expDF=expDF)
    return(gp)
}

#-----------------------------------------------------------------------
# parse matirx-scan output file and creates a table with the number of
# motif hits per gene (rows) for each motif (columns).
#-----------------------------------------------------------------------
parseMatrixScan <- function(matrixScanFile, geneIDs){

    classes <- sapply(read.table(matrixScanFile, nrows = 5, comment.char=c(';'), header=TRUE), class)
    
    motifHits = read.delim(matrixScanFile, comment.char=';', header=TRUE, colClasses = classes)
    

    # initialize zero matrix with to count motif instances per gene promoter
    motifs = unique(motifHits$ft_name)
    motifTable = matrix(0, nrow=length(geneIDs), ncol=length(motifs), dimnames=list(geneIDs, motifs))
    
    hitGenes = motifHits[,1]
    hitMotifs = motifHits$ft_name
    for (i in 1:nrow(motifHits)){
        motifTable[hitGenes[i],hitMotifs[i]] = motifTable[hitGenes[i],hitMotifs[i]]+1
    }

    return(data.frame(motifTable))
}
#matrixScanFile="results/paralog_regulation/EnsemblGRCh37_paralog_genes.promoters.bed.names.fa.jaspar.matrix-scan.uniquePos"

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
# get pairwise information from gene matrix (like promoter contacts from capture Hi-C)
#-----------------------------------------------------------------------
getPairwiseMatrixScore <- function(genePairs, M, tssGR, replaceZeroByNA=FALSE){
    
    #scores = M[genePairs[,1], genePairs[,2]]
    idx1 = match(genePairs[,1], id(tssGR))
    idx2 = match(genePairs[,2], id(tssGR))
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
# for each gene pair get the gene ID of unique common ortholog genes
#-----------------------------------------------------------------------
getUniqueCommonOrthologs <- function(genePairs, orthologsSpecies, orthoGeneCol){
    
    # convert gene ID to character
    genePairs[,1] = as.character(genePairs[,1])
    genePairs[,2] = as.character(genePairs[,2])

    # iterate over all gene pairs gP
    commonOrthologsSpecies = apply(genePairs[,1:2], 1, function(gP){
        
        # get the set of orthologouse genes to each single gene in the pair
        ortho1 = orthologsSpecies[unlist(gP[1]) == orthologsSpecies[,1], orthoGeneCol]
        ortho2 = orthologsSpecies[unlist(gP[2]) == orthologsSpecies[,1], orthoGeneCol]
        
        # check if for each gene only one ortholog was found and if they are identical
        if(length(ortho1)==1 & length(ortho2)==1 & length(intersect(ortho1, ortho2)) == 1){
            return(ortho1)
        }else{
            # if this was not the case return NA
            return(NA)
        }
        
    })
    
    return(commonOrthologsSpecies)
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
# get a boolean vector indicating whether (in case of same strand) the first 
# gene in the pair is upstream of the second.
#-----------------------------------------------------------------------
getFirstUpstream <- function(genePairs, tssGR){
    
    g1 = tssGR[genePairs[,1]]
    g2 = tssGR[genePairs[,2]]

    firstUpstream = (strand(g1) == "+" & strand(g2) == "+" & start(g1) <= start(g2)) | (strand(g1) == "-" & strand(g2) == "-" & start(g1) > start(g2))
    firstDownstream = (strand(g1) == "+" & strand(g2) == "+" & start(g1) > start(g2)) | (strand(g1) == "-" & strand(g2) == "-" & start(g1) <= start(g2))
    
    boolVec = rep(NA, nrow(cisPairs))
    boolVec[as.logical(firstUpstream)] = TRUE
    boolVec[as.logical(firstDownstream)] = FALSE
    return(boolVec)
}


#-----------------------------------------------------------------------
# get percent identety for both genes
#-----------------------------------------------------------------------
getPercentID <- function(g, o, orthologsSpecies, orthoGeneCol= "mmusculus_homolog_ensembl_gene", simCol="mmusculus_homolog_perc_id_r1"){
    
    if( is.na(o) ){ 
        return(NA)
    }
    uniqOrtholog = g == orthologsSpecies[,1] & o == orthologsSpecies[,orthoGeneCol]
    
    sim = orthologsSpecies[uniqOrtholog, simCol][1]

    return(sim)
}

#-----------------------------------------------------------------------
#
#-----------------------------------------------------------------------
runOrientationAnalysis <-function(genePairs, tssGR, orthologsSpecies, species, seqSimSuffix="_homolog_perc_id_r1"){

    commonOrthologsSpecies = getUniqueCommonOrthologs(cisPairs, orthologsSpecies, orthoGeneCol=paste0(species, "_homolog_ensembl_gene"))
    
    
    # get similarity of ortholog to first gene in pair
    g1_orthoSim = unlist(mapply(getPercentID, genePairs[,1], commonOrthologsSpecies, MoreArgs=list(orthologsSpecies=orthologsSpecies, orthoGeneCol= paste0(species, "_homolog_ensembl_gene"), simCol=paste0(species, seqSimSuffix))))
    
    g2_orthoSim = unlist(mapply(getPercentID, genePairs[,2], commonOrthologsSpecies, MoreArgs=list(orthologsSpecies=orthologsSpecies, orthoGeneCol= paste0(species, "_homolog_ensembl_gene"), simCol=paste0(species, seqSimSuffix))))
    
    
    # annotate gene paris with the locical label if the first of the pair is upstream of the other
    firstUpstream = getFirstUpstream(genePairs, tssGR)
    
    # get similarity values for upsteam and downstream homolog
    upstream = ifelse(firstUpstream, g1_orthoSim ,g2_orthoSim)
    downstream = ifelse(firstUpstream, g2_orthoSim ,g1_orthoSim)
    
    # n1 means that the upsteam paralog is more similar to the ortholog
    upstreamHigher = upstream > downstream
    
    # check if upstream or downstream is more similar to ortholog
    downstreamHigher = upstream < downstream 
    eqSim = upstream == downstream
    
    
    nPairs = c(
        "close_pairs"=nrow(genePairs), 
        "+---same_strand"=sum(genePairs$sameStrand), 
        "+---common_ortholog"=sum(!is.na(commonOrthologsSpecies)), 
        "    +---same_strand"=sum(!is.na(commonOrthologsSpecies) & cisPairs$sameStrand), 
        "        +---equal_sim"=sum(eqSim, na.rm=TRUE), 
        "        +---upstream_conserved"=sum(upstreamHigher, na.rm=TRUE), 
        "        +---downstream_conserved"=sum(downstreamHigher, na.rm=TRUE) 
        )
    
    return(list(nPairs, upstream, downstream))

}

########################################################################
# OLD and unused stuff:
########################################################################

queryGenePairs <- function(genePairs, g){
    (g == genePairs[,1]) | (g == genePairs[,2])
}

#-----------------------------------------------------------------------
# get only one pair per unique gene by choosing the gene pair with highest 
# sequence similarity
#-----------------------------------------------------------------------
uniquePairPerGeneBySimIteratively <- function(gP, similarity){
    
    genePairs = gP
    genePairs[,1] = as.character(genePairs[,1])
    genePairs[,2] = as.character(genePairs[,2])
    
    # iterate over all genes g 
    #   get all pairs p_g with gene g
    #       While size(p_g) > 1:
                # get a-b pair with highest sim
                # remove all pars with a or b from p_g
    
    # collect final pairs indexes
    resultPairs = c()
    
    for (g in unique(c(genePairs[,1], genePairs[,2]))){
        
        PgQuery = queryGenePairs(genePairs, g)
        
        #Pg = genePairs[PgQuery,]
        
        while(sum(PgQuery) > 1 ){
            # get minimum of similarity (multiplied by true/false query to get subset of pairs wich involve gene g)
            sim = similarity
            sim[!PgQuery] = NA
            idx = which.min(sim)
            
            resultPairs = c(resultPairs, idx)
            
            a = genePairs[idx,1]
            b = genePairs[idx,2]
            
            # set pairs that ivolve a and b to false for this query
            aPairs = queryGenePairs(genePairs, a)
            bPairs = queryGenePairs(genePairs, b)
            PgQuery[aPairs | bPairs] = FALSE
        }
        
    }
    
    # TODO: finish this one to see fi results differ
    # Problem: order of gene choice matters. If Pair A-B is choosen to keep pair B-C with potentially higher similarity will be removed.
    
    return(uniqPairs)
}





