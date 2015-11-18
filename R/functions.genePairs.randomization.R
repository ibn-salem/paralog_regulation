########################################################################
#
# A script provides functionality to sample random gene pairs by 
# maintaining some properties, such as linear distance and shared enhancer. 
#
# It was initially written as for the paralog_regulation.R analysis.
#
########################################################################

require(stringr)        # for some string such as paste0()
require(GenomicRanges)

# New and faster approach
# 1.) precompute all possible gene pairs within a distance of XMb
#   ID_g1   ID_g2   dist    eh1     eh2
#   3       5       -2341   2       1
#   3       1       3234    2       3
#   ...
# Do this by extending each TSS by +-X/2 Mb and get overlapping TSS

# 2.) Sample according to distance distribution and number of enhancers in source genes
# Precalculate the density d_eh(g) for all number of linked enhancers eh(g) in the 
# set of source genes.
# For each source pair take d_eh(g1,g2) = [d(eh(g1)) * d(eh(g2))] as density.
# Analogusly get the distance distribution density d_dist(pairs) from source pairs.
# Now, from all possible gene pairs sample according to combined weight
# w_p = d_eh(pair) * d_dist(pair)
# 

#-----------------------------------------------------------------------
# get all possible gene pairs within a given maximal distance as data frame
# with annotations like distance and linked enhancers
#-----------------------------------------------------------------------
getAllGenePairs <- function(tssGR, maxDist, minDist=0){

    # calculate overlap all possible gene pairs within maxDist bp
    hits = findOverlaps(tssGR, resize(tssGR, 2*maxDist+1, fix="center", ignore.strand=TRUE), ignore.strand=TRUE)
    
    # convert hits objec to data.frame
    hitDF = as.data.frame(hits)
    names(hitDF) <- c("IDg1", "IDg2")

    # remove pairs with same gene 
    hitDF = hitDF[hitDF$IDg1 != hitDF$IDg2, ]
    
    # add distance
    hitDF$dist = start(tssGR[hitDF$IDg2]) - start(tssGR[hitDF$IDg1])

    # remove pars with distance smaller than minDist
    hitDF = hitDF[abs(hitDF$dist) >= minDist,]
    
    # annotate the gene pairs with number of enhancers
    hitDF$EHg1 = tssGR[hitDF$IDg1]$linked_enhancer
    hitDF$EHg2 = tssGR[hitDF$IDg2]$linked_enhancer
    
    return(hitDF)

}

#-----------------------------------------------------------------------
# To sample random gene pairs get probability weights of all pairs 
# according to the distribution of distances and linked enhancer in input
# gene pair set
#-----------------------------------------------------------------------
getSampleWeightsByDistAndEnhancers <- function(hitDF, tssGR, sourcePairs, ...){
    
    # get probability density for the number of linked enhancers
    weightEH = weightByEnhancers(tssGR, c(sourcePairs[,1], sourcePairs[,2]))
    
    # combine the probabilites from both gene pair partners by assuming independence
    pairWeightEH = ( weightEH[hitDF[,1]] * weightEH[hitDF[,2]] )
    pairWeightEH = pairWeightEH / sum(pairWeightEH)
    
    # get observed probabilty density of distances in source pairs
    distDens = approxfun(density(abs(sourcePairs$dist), ...))
    pairWeightDist = distDens(abs(hitDF$dist))

    # if probability is NA set it to 0
    pairWeightDist[is.na(pairWeightDist)] = 0
    pairWeightDist = pairWeightDist / sum(pairWeightDist)
    
    # combine the probabilities for linked enhancer and linear distance by multiplication by assuming independence
    pairWeight = ( pairWeightEH * pairWeightDist ) 
    pairWeight  = pairWeight / sum(pairWeight)
    
    return(pairWeight)
}


#-----------------------------------------------------------------------
# Get sampling weight only based on distance (or other column by "colName")
#-----------------------------------------------------------------------
getSampleWeightsByDist <- function(hitDF, sourcePairs, colName="dist", ...){

    # get observed probabilty density of distances in real source pairs
    distDens = approxfun(density(abs(sourcePairs[,colName]), ...))
    pairWeightDist = distDens(abs(hitDF[,colName]))

    # set prob of NAs to zero
    pairWeightDist[is.na(pairWeightDist)] = 0
    return(pairWeightDist / sum(pairWeightDist, na.rm=TRUE))
}

#-----------------------------------------------------------------------
# sample pairs from the input gene pair set according to given sampling weights
#-----------------------------------------------------------------------
sampleFromAllPairsByWeight <- function(n, hitDF, tssGR, weight){
        
    # sample random pairs from input pairs according to weights
    pairIDX <- sample.int(nrow(hitDF), n, prob=weight, replace=TRUE)
    
    # construct gene pair DF from sampled pairs
    rP <- data.frame(
        g1=as.character(names(tssGR)[hitDF[pairIDX,1]]),
        g2=as.character(names(tssGR)[hitDF[pairIDX,2]]),
        dist=hitDF[pairIDX, "dist"],
    stringsAsFactors=FALSE)
    
    return(rP)
    
}
# rP = sampleFromAllPairs( n=100, hitDF=allGenePairs,tssGR, sourcePairs=cisPairs[abs(cisPairs$dist) <= MAX_DIST, ], sourceGenes=paralogs[,1])
         

#-----------------------------------------------------------------------
# sample randomly from vector, even from size one vector with integer elements
#-----------------------------------------------------------------------
sample.vec <- function(x, ...) x[sample(length(x), ...)]

#-----------------------------------------------------------------------
# get weights for all genes for sampeling according to the distribution 
# of linked enhancers in the set of sourceGenes (paralog genes)
#-----------------------------------------------------------------------
weightByEnhancers <- function(tssGR, sourceNames){

    # non-adjusted frequencies of enhancers in all genes:
    allGenesEnhancerFreqTable = table(mcols(tssGR)[, "linked_enhancer"])
    
    # get distribution function of linked enhancer elements in the paralog genes
    paralogsEnhancerFreqTable = table(mcols(tssGR[sourceNames])[,"linked_enhancer"])
    
    # convert the numbers to character
    ehCounts = as.character(mcols(tssGR)[, "linked_enhancer"])
    
    # take ratio of enhancers in paralogs and enhancers in all genes as weight for randomly sampling from all genes
    weight = paralogsEnhancerFreqTable[ehCounts] / allGenesEnhancerFreqTable[ehCounts]
    
    # remove NA's, e,g, number of enhancers not observed in paralogs but in set of all genes. Set their probability to zero 
    weight[is.na(weight)] = 0
    
    # normalize weights to 1
    weight = weight / sum(weight)
    
    return(weight)
}
# weight = weightByEnhancers(tssGR, paralogs[,1])
# randGenes = sample(genes[,1], nrow(genes), prob=weight, replace=TRUE)

#-----------------------------------------------------------------------
# sample randomly pairs from the entire gene set with equal probabilities
#-----------------------------------------------------------------------
getRandomPairs <- function(n, geneIDs){
    randomPairs = data.frame(
        t(replicate(n, sample(geneIDs, size=2, replace=FALSE))), 
        stringsAsFactors=FALSE)
    names(randomPairs) = c("g1", "g2")
    return(randomPairs)
}
