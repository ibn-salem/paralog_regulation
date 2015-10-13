########################################################################
#
# Some functions to analyse Hi-C chromatin interactions.
# Depend mostly on HiTC package from Bioconductor
#
########################################################################
require(HiTC)
require(igraph)     # to group pairs of overlapping TADs

#-----------------------------------------------------------------------
# Compute fraction of overlap (reciprocal) 
#-----------------------------------------------------------------------
reciprocalOverlaps <- function(query, subject, fraction=.5){

    # get any overlap hit object
    hits <- findOverlaps(query, subject, type="any")

    # compute width of the overlapping region
    overlaps <- width(pintersect(query[queryHits(hits)], subject[subjectHits(hits)]))
    
    # get max length of query and subject for each overlap
    maxLen <- apply(
        cbind(width(query[queryHits(hits)]), width(subject[subjectHits(hits)])),
        1,max)

    # compute fraction of overlap
    fracOverlap <- overlaps / maxLen
    
    # filter for overlap >= threshold
    return(hits[fracOverlap >= fraction])

}

#-----------------------------------------------------------------------
# Compute overlap by allowing a specific resolution of start and end matches
#-----------------------------------------------------------------------
resolutionOverlap <- function(query, subject, resolution=10^4){

    # get any overlap hit object
    hits <- findOverlaps(query, subject, type="any")
    
    qHits <- query[queryHits(hits)]
    sHits <- subject[subjectHits(hits)]
    
    # get differences for start and end coordinates from overlap pairs
    startDiff <- abs(start(qHits) - start(sHits))
    endDiff <- abs(end(qHits) - end(sHits))
    
    return(hits[startDiff <= resolution & endDiff <= resolution])
}

#-----------------------------------------------------------------------
# get most common integer from a vector
#-----------------------------------------------------------------------
mostCommon <- function(v){
    as.numeric(names(sort(table(v), decreasing=TRUE)[1]))
}

#-----------------------------------------------------------------------
# get conserved TADs
# using binary clustering to get connected components:
# http://stackoverflow.com/questions/30407769/get-connected-components-using-igraph-in-r
# 
#-----------------------------------------------------------------------
getConservedTADs <- function(TADList, n=3, fraction=.8){
    
    # make GRange object of all TADs in all cell types
    allTADsGR = unlist(GRangesList(TADList))
        
    # get all pairwise overlaps with reciprocal fraction >= threshold
    hits = reciprocalOverlaps(allTADsGR, allTADsGR, fraction=fraction)
    
    # get groups.
    g <- graph_from_edgelist(as.matrix(hits), directed=FALSE)
    cl <- clusters(g)
    
    nodes <- as.vector(V(g))
    
    # get list of TAD indices by groups that overlap each other
    groups = bplapply(seq(cl$no)[cl$csize > n], function(k){ nodes[cl$membership %in% k] })
    
    consTADs <- unlist(GRangesList(bplapply(groups, function(gr){
        
        TADs <- allTADsGR[gr]
        
        # take most common start and and positions for consensus Range
        GRanges(seqnames(TADs)[1], IRanges(mostCommon(start(TADs)), mostCommon(end(TADs))), seqinfo=seqinfo(TADs))
        
    })))
    
    return(consTADs)
}

#-----------------------------------------------------------------------
# get conserved by providing the (filtered) all against all hit object
#-----------------------------------------------------------------------
getConservedByHits <- function(hits, n=3){
    # get groups.
    g <- graph_from_edgelist(as.matrix(hits), directed=FALSE)
    cl <- clusters(g)
    
    nodes <- as.vector(V(g))
    
    # get list of TAD indices by groups that overlap each other
    groups = bplapply(seq(cl$no)[cl$csize > n], function(k){ nodes[cl$membership %in% k] })
    
    consTADs <- unlist(GRangesList(bplapply(groups, function(gr){
        
        TADs <- allTADsGR[gr]
        
        # take most common start and and positions for consensus Range
        GRanges(seqnames(TADs)[1], IRanges(mostCommon(start(TADs)), mostCommon(end(TADs))), seqinfo=seqinfo(TADs))
        
    })))
    
    return(consTADs)
}   

#-----------------------------------------------------------------------
# Write HiCT object to long range interaction format
#-----------------------------------------------------------------------
writeHiCinteractions <- function(HiCList, outFile){
    
    allTable = rbind()
    
    for(subMap in HiCList){

        allPairIDX = t(combn(1:nrow(intdata(subMap)), 2))
        n = nrow(allPairIDX)
    
        scores = intdata(subMap)[allPairIDX]
        reg1 = x_intervals(subMap)[allPairIDX[,1]]
        reg2 = y_intervals(subMap)[allPairIDX[,2]]
    
    
        tab1 = cbind(as.character(seqnames(reg1)), start(reg1), end(reg1),
         paste0(getLocStr(reg2), ",", scores), 1:n, ".")
        
        # make reverse interactions (required by track fromat)
        tab2 = cbind(as.character(seqnames(reg2)), start(reg2), end(reg2),
            paste0(getLocStr(reg1), ",", scores), (n+1):(2*n), ".")
        
        allTable = rbind(allTable, tab1, tab2)
    }
    # write table to output file
    write.table(allTable, outFile, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

}


#-----------------------------------------------------------------------
# returns genome-wide border strength as RleList object
#-----------------------------------------------------------------------
borderStrengthRle <- function(HiClist, nBins=10, seqInfo){
    
    allBins = GRanges(seqinfo=seqInfo)
    
    # iterate over all chromosomes
    for (chr in seqnames(seqInfo)){
        print(paste("DEBUG: chr: ", chr))
        # first check if for this chromosome ther is a Hi-C map available:
        if (paste0(chr, chr) %in% names(HiClist)){

            # get intra chromosomal experiment
            hiCexp = HiClist[[paste0(chr, chr)]]
            
            # get all bins
            chrBins = x_intervals(hiCexp)
            
            # add seqInfo with chromosome lengths
            seqlevels(chrBins) = seqlevels(seqInfo)
            seqinfo(chrBins) = seqInfo
            
            # compute borderstrength for each bin on the entire chrom
            score(chrBins) = intraInterLogRatio(intdata(hiCexp), nBins, ps.count=1)
    
            # add ranges to all ranges
            allBins = c(allBins, chrBins)
        }
    }
    print(paste("DEBUG: Finished run!"))
#~     return(allBins)
    
    # compute genome wide score as RleList
    return(coverage(allBins, weight="score"))
}
# bsRle = borderStrengthRle(HiClist, 10, seqInfo[c("chr21", "chr22")])


#-----------------------------------------------------------------------
# Calculate the border strength (Van Bortel et al. 2014 Genome Biology)
# as fraction of local intra- versus inter-domain interaction frequencies
# for all bins overlapping an input query region
#-----------------------------------------------------------------------
getBorderStrength <- function(gr, hiCexp, nBins=10){ 
    
    # assume only one range as input
    stopifnot(length(gr) == 1)
    
    chr = as.character(seqnames(gr))
    
        
    # get indices of bins that overlap a region in gr
    binsIndices  = which(countOverlaps(x_intervals(hiCexp), gr) >= 1)
    
    # do some looping overall bins with current bin index s    
    sapply(binsIndices, intraInterLogRatio, nBins, intdata(hiCexp))
}
#getBorderStrength(tadGR[1], HiClist[["chr22chr22"]], 10)

#-----------------------------------------------------------------------
# Calculate the mean border strength values for each range by taking the 
# mean over all overlapping bins from the Hi-C experiment hiCexp
#-----------------------------------------------------------------------
getMeanBorderStrength <- function(inputGR, HiClist, nBins=10){
    
    # initialize vector for border stength
    borderStrength = rep(NA, length(inputGR))
    
    #iterate over all chromosomes
    for ( chr in as.character(unique(seqnames(inputGR))) ){
        
        # get intra chromosomal matrix for this chrom
        hiCexp = HiClist[[paste0(chr, chr)]]

        subSet = as.character(seqnames(gr)) == chr
        
        borderStrength[subSet] =  sapply(inputGR[subSet], function(gr) mean(getBorderStrength(gr, hiCexp, nBins)))

    }
    return(borderStrength)
}
#getMeanBorderStrength(tadGR, HiClist, 10)

#-----------------------------------------------------------------------
# calculate the ratio between intra and inter domain interactions for 
# all bins along the diagonal of the interaction matrix M by taking nBins 
# as region on both sides of s into account for intra and inter region 
# contact calculations. 
#-----------------------------------------------------------------------
intraInterLogRatio <- function(M, k=10, ps.count=1){
    
    n = nrow(M)
    stopifnot(k < n)
    
    # get intra-bin interactions along the diagonal
    intraB = getIntraIntervalContacts(M, k)
    # shift these by k bins 
    intraA = c(rep(NA, k), intraB[1:(n-k)] )
    
    # get for each adjacent regios A and b (each of k bins) the number of 
    # inter-interacionts between A and B
    inter = getInterIntervalContacts(M, k)

    
    log2((intraA + intraB + ps.count ) / (inter + ps.count))
}
#bs = intraInterLogRatio(M=intdata(HiClist[["chr22chr22"]]), k=10)

#-----------------------------------------------------------------------
# Calculates the sum of interaction frequencies within a sliding window 
# of size k bins along the diagonal of the interaction matrix M.
# TODO: improve speed of this function
#-----------------------------------------------------------------------
getIntraIntervalContacts <- function(M, k=10){

    n = nrow(M)
    stopifnot(k < n)
    M[lower.tri(M)] = NA

    # initialize vector
    I = rep(NA, n)

    # start by directly get the all values from the matrix
    I[1] = sum(M[seq(1,k), seq(1,k)], na.rm=TRUE)
    for (s in 2:n){
        
        I[s] =  I[s-1] - sum(M[s-1, seq(s-1, min(s+k-2, n))], na.rm=TRUE) + sum(M[seq(s, min(s+k-1,n)), min(s+k-1, n)], na.rm=TRUE)
    
    }
    return(I)
}
#I = getIntraIntervalContacts(M=intdata(HiClist[["chr22chr22"]]), k=10)
    
#~     nextS <- function(s) {
#~             if (s==1) sum(M[seq(1,k), seq(1,k)], na.rm=TRUE)
#~             else
#~             I[s-1] - sum(M[s-1, seq(s-1, min(s+k-2, n))], na.rm=TRUE) + sum(M[seq(s, min(s+k-1,n)), min(s+k-1, n)], na.rm=TRUE)
#~     }
#~     
#~     I = rapply(list(a=1:n), nextS)
    
#~     I[2:n] = sapply(2:n, function(s) I[s-1] - sum(M[s-1, seq(s-1, min(s+k-2, n))], na.rm=TRUE) + sum(M[seq(s, min(s+k-1,n)), min(s+k-1, n)], na.rm=TRUE))

#-----------------------------------------------------------------------
# 
#-----------------------------------------------------------------------
getInterIntervalContacts <- function(M, k=10){

    n = nrow(M)
    stopifnot(k < n)
    M[lower.tri(M)] = NA
    
    # initialize vector
    I = rep(NA, n)

    #  loop over all bins (along diagonal of matrix)
    for (s in seq(k+1, n-k+1)){
        
        a = seq(s-k, s-1)
        b = seq(s, s+k-1)
        
        I[s] = sum(M[a,b], na.rm=TRUE)
        
    }
    return(I)
}
#I = getInterIntervalContacts(M=intdata(HiClist[["chr22chr22"]]), k=10)


########################################################################
# POTENTIALLY OLD STUFF:
########################################################################
#-----------------------------------------------------------------------
# calculate the ratio between intra and inter domain interactions for bin
# with index s by taking nBins as 'intra' region on both sides of s into account. 
#-----------------------------------------------------------------------
intraInterLogRatioOLD <- function(s, nBins, mat, ps.count=1){

    n = nrow(mat)
    # bins in the left region
    a = seq(max(0, s-nBins), min(s-1, n))
    # bins in the right region
    b = seq(max(0, s), min(s+nBins-1, n))
    #~     a = mapply(seq, s-nBins, s-1)
    #~     b = mapply(seq, s, s+nBins-1)
    
    inter = sum(mat[a,b]) + ps.count 
    intra = sum(mat[a,a]) + sum(mat[b,b]) + ps.count
    
    # return the log ratio
    return(log2(intra/inter))
    
}


#-----------------------------------------------------------------------
# Calculate the border strength (Van Bortel et al. 2014 Genome Biology)
# as fraction of local intra- versus inter-domain interaction frequencies
#-----------------------------------------------------------------------
borderStrengthOLD <- function(hiCexp, nBins=10){ 

    bins = x_intervals(hiCexp)
    n = length(bins)
    mat = intdata(hiCexp)
    
    # get bins for which the calculation is possible:
    #TODO: fix this to calc it for all bins
    #s = seq(nBins, n-nBins)
    
    
    # do some looping overall bins with current bin index s    
    ratio = sapply(1:n, intraInterLogRatio, nBins, mat)

    return(ratio)
}

#-----------------------------------------------------------------------
# computes the sum of local intra-window interactions along a sliding window
# where the window has size w corresponding to the number of bins
# Furthermore, it adds smaller regions at the lower end of the matrix to the vector
# this results in addional w values such that the total size of the
# return vector will be n+w instead of n.
#-----------------------------------------------------------------------
localInteractionSum <- function(mat, w){

    # get number of bins
    n <- nrow(mat)
    
    # initialize score s with 0
    s <- 0
    scores <- c()

    for (i in 1:n) {
    
        # as long as i <= n
        #if (i <= w ){
        upLim = max(1, i-w+1)

        # add sum of interactions of the next bin
#~         s = s + sum(mat[upLim:i, i])
        s = s + mean(mat[upLim:i, i])

        if (i > w){
            # subtract interactions of old bin
#~             s = s - sum(mat[i-w, (i-w):(i-1)])
            s = s - mean(mat[i-w, (i-w):(i-1)])
        }


        # add score to results
        scores <- c(scores, s)    
    }

    # add smaller regions at the lower end of the matrix to the vector
    # this results in addional w values such that the total size of the
    # return vector will be n+w instead of n.
    
    for (k in (n-w+1):(n)){

        # subtract interactions of old bin
#~         s = s - sum(mat[k, k:n])
        s = s - mean(mat[k, k:n])

        scores <- c(scores, s)    
    }
    return(scores)
}

# DEBUG: localInteractionSum(mini, 3)

#-----------------------------------------------------------------------
# Get sum of interactions between the tow flaning regions a and b of each
# position s along the genome
#-----------------------------------------------------------------------
localCrossingInteractions <- function(mat, w){

    # get number of bins
    n <- nrow(mat)
    
    # initialize score s with 0
    s <- 0
    scores <- c()

    for (i in 1:(n-1)) {
        
        # left flanking region (including current position)
        a = seq(max(1, i-w+1), i)
        # right flanking region
        b = seq(i+1, min(i+w, n))

        # add sum of interactions of the next bin
#~         s = sum(mat[a,b])
        s = mean(mat[a,b])

        scores <- c(scores, s) 
    }
    
    # for last bin add 0 
    scores <- c(scores, 0) 

    return(scores)
}

# localCrossingInteractions(m,2)

#-----------------------------------------------------------------------
# get fraction of intra- versus inter interactions based on pre-computed
# vecotrs with fixed window size
#-----------------------------------------------------------------------
logRatioIntraInter <- function(b, localInteractions, crossingInteractions, w){

    intra = localInteractions[b] + localInteractions[b+w]
    
    inter = crossingInteractions[b]
    
    return(log( (intra +1) / (inter+1) ) )

}

