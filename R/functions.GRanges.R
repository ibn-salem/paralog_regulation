#!/usr/bin/Rscript
#=======================================================================
#
#   Functions related to GRanges objects from the GenomicRanges package
#
#=======================================================================

# load some useful libraries
require(stringr)
require(GenomicRanges)
require(igraph)     # to group pairs of overlapping TADs


#-----------------------------------------------------------------------
# read chromosom names and length as Seqinfo object
#-----------------------------------------------------------------------
parseSeqInfo <- function(chromFile, genome="hg19", header=TRUE, filterForReal=FALSE){
    # read chromosom names and length as Seqinfo object
    chroms = read.delim(chromFile, header=header)
    
    if(filterForReal){
        realChromIdx =  grep("chr[[:digit:]XY]+$", chroms[,1])
        chroms = chroms[realChromIdx,]
    }
#~     # order chromosomes by names
#~     chroms = chroms[order(paste0(chroms[,1], "$")),]
    seqInfo <- Seqinfo(seqnames=as.character(chroms[,1]), seqlength=chroms[,2], isCircular=rep(FALSE, nrow(chroms)), genome=genome)
}

#-----------------------------------------------------------------------
# make GRange object from downloaded ENSEMBLE table of genes
#-----------------------------------------------------------------------
getGenesGR <- function(genes, seqInfo, annotCols=c("hgnc_symbol", "status", "gene_biotype")){
    gr = GRanges(
        paste0("chr", genes$chromosome_name),
        IRanges(genes$start_position, genes$end_position),
        strand = ifelse(genes$strand == 1, '+', '-'), 
        names = genes$ensembl_gene_id, 
        genes[,annotCols],
        seqinfo=seqInfo
        )
    names(gr) = genes$ensembl_gene_id
    return(gr)
}


#-----------------------------------------------------------------------
# create a Grange object for uniue TSS
#-----------------------------------------------------------------------
getTssGRfromENSEMBLGenes <- function(genes, seqInfo, colNames=c("hgnc_symbol", "status", "gene_biotype")){
    
    if (length(unique(genes$ensembl_gene_id)) != length(genes$ensembl_gene_id)){
        message("WARNING: ENSG IDs are not unique!")
    }
    
    # get TSS as start of the gene if on + strand and end of the gene else 
    onPosStrand <- genes$strand == 1
    tss = ifelse(onPosStrand, genes[, "start_position"], genes[, "end_position"])

    gr = GRanges(
        paste0("chr", genes$chromosome_name),
        IRanges(tss, tss),
        strand = ifelse(onPosStrand, '+', '-'), 
        seqinfo=seqInfo
        )
    names(gr) = genes$ensembl_gene_id
    mcols(gr) = subset.data.frame(genes, select=colNames)
    return(gr)

}

#-----------------------------------------------------------------------
# create a Grange object for uniue gene centers
#-----------------------------------------------------------------------
getCenterAsGR <- function(genes, seqInfo, colNames=c("hgnc_symbol", "status", "gene_biotype")){
    
    if (length(unique(genes$ensembl_gene_id)) != length(genes$ensembl_gene_id)){
        message("WARNING: ENSG IDs are not unique!")
    }
    
    onPosStrand <- genes$strand == 1

    # get center as midpoint between start and end of the gene 
    center = genes[, "start_position"] + round((genes[, "end_position"] - genes[, "start_position"]) / 2)

    gr = GRanges(
        paste0("chr", genes$chromosome_name),
        IRanges(center, center),
        strand = ifelse(onPosStrand, '+', '-'), 
        seqinfo=seqInfo
        )
    names(gr) = genes$ensembl_gene_id
    mcols(gr) = subset.data.frame(genes, select=colNames)
    return(gr)
}


#-----------------------------------------------------------------------
# create a Grange object for uniue TSS
#-----------------------------------------------------------------------
getExonGR <- function(exonDF, seqInfo, annotCols=c("ensembl_exon_id", "rank", "gene_biotype")){
    GRanges(
        paste0("chr", exonDF$chromosome_name),
        IRanges(exonDF$exon_chrom_start, exonDF$exon_chrom_end),
        strand=ifelse(exonDF$strand == 1, '+', '-'), 
        exonDF[,annotCols],
        seqinfo=seqInfo
    )
}

#-----------------------------------------------------------------------
# Filter for only those ranges that do not overlap other ranges
#-----------------------------------------------------------------------
filterForNonOverlapping <- function(gr){
    gr[countOverlaps(gr, gr) == 1]
}


#-----------------------------------------------------------------------
# Downsample a GRanges object
#-----------------------------------------------------------------------
downSampleGR <- function(gr, n){
    
    if (n >= length(gr)){
        return(gr)
    }
    
    gr[sample.int(length(gr), n)]
}

#-----------------------------------------------------------------------
# Get the locus string chr:start-end form a GRanges object
#-----------------------------------------------------------------------
getLocStr <- function(gr){
    paste0(seqnames(gr), ':', start(gr), '-', end(gr))
}

#-----------------------------------------------------------------------
# returns the input GenomicRanges object without strand annotation
#-----------------------------------------------------------------------
noStrand <- function(gr){
    strand(gr) = '*'
    return(gr)
}

#-----------------------------------------------------------------------
# filters out all ranges with '+' or '-' strand an leaves only '*' strand
#-----------------------------------------------------------------------
filterForNoStrand <- function(gr) gr[strand(gr) == '*']

#-----------------------------------------------------------------------
# Get sum of bases covered by at least on element in gr
#-----------------------------------------------------------------------
covBases <- function(gr){
    sum(as.numeric(width(reduce(noStrand(gr)))))
}



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
# Transpose a list of n GRAnges with m equal sized bins into list of all the ith bin 
#-----------------------------------------------------------------------
transposeBinList <- function(TADbins, m=20){
    allBins = unlist(TADbins)
    n = length(allBins)
    lapply(1:m, function(i) allBins[seq(i,n,m)])    
}   

#-----------------------------------------------------------------------
# get bins around TAD boundaries
#-----------------------------------------------------------------------
binsAroundTADboundaries <- function(TAD, windowSize=20*10^3, nbins=20){
    
    # get boundaries of TADs as GRanges object of widht=1
    boundaries = c(
        GRanges(seqnames(TAD), IRanges(start(TAD), start(TAD))),
        GRanges(seqnames(TAD), IRanges(end(TAD), end(TAD)))
    )
    # get region around boundaries of window size
    boundaryRegion = resize(boundaries, width=windowSize, fix="center", ignore.strand=TRUE)
    
    # make n equal sized bins for each boundary region 
    bins= tile(boundaryRegion, n=nbins)
    
    # get indices of the left boundary of each TAD (made from start coordinates)
    leftBoundary = 1:length(TAD)
    
    # flip the bin order in left TAD boundary Region to have always the part inside TAD on the left hand-side
    bins[leftBoundary] = GRangesList(lapply(bins[leftBoundary], rev))
    
    return(bins)
    
}

#-----------------------------------------------------------------------
# Plot distribution around boundaries of TADs
# (Aaccording to: https://stat.ethz.ch/pipermail/bioconductor/2013-February/050841.html)
#-----------------------------------------------------------------------
coverageInBins <- function(gr, binList){
    
    # compute coverage of gr on whole genome
    cov = coverage(gr)
    
    lapply(binList, function(binGR) sum(cov[binGR]))

}


#-----------------------------------------------------------------------
# Count number of overlaps for each bin
#-----------------------------------------------------------------------
countOverlapsInBins <- function(gr, binList, ...){
    
    # for each i-th bins count the number of overlapping elements
    lapply(binList, function(binGR) countOverlaps(binGR, gr, ...))
        
}

#-----------------------------------------------------------------------
# Coverage Profile along TAD boundaries
# retunrs a list with a nbins long vector for each function in combFunctions
#-----------------------------------------------------------------------
countsAroundTADboundaries <- function(binList, grl, combFun=mean, normalizeToRegions=FALSE){
            
    # get coverage for each bin
    #~binCov = coverageInBins(gr, binList)
    binCovList = lapply(grl, countOverlapsInBins, binList, type="any", ignore.strand=TRUE)
    message("INFO: Finished overlap calculations.")

    # combine coverage for all i-th bins with the input functions
    combBinCov = lapply(binCovList, function(binCov) sapply(binCov, combFun))
    message("INFO: Finished combining of counts over all individual TADs.")
    
    if (normalizeToRegions) {
        # normalize counts to number of ranges in each region set
        combBinCov = mapply(function(counts, gr){counts/length(gr)}, combBinCov, grl, SIMPLIFY=FALSE)
    }
    
    return(combBinCov)
}

#-----------------------------------------------------------------------
# Count number of overlaps for each bin
#-----------------------------------------------------------------------
countOverlapsInBinsStranded <- function(gr, binList, flippedBins = 1:(length(binList[[1]])/2), ...){

    # for each i-th bins count the number of overlapping elements
    posCounts = lapply(binList, function(binGR){
        countOverlaps(binGR[-flippedBins], gr[strand(gr) %in% c('+', '*')], ...) + 
        countOverlaps(binGR[flippedBins], gr[strand(gr) == '-'], ...) 
            
    })
    negCounts = lapply(binList, function(binGR){ 
        countOverlaps(binGR[-flippedBins], gr[strand(gr) == '-'], ...) + 
        countOverlaps(binGR[flippedBins], gr[strand(gr) %in% c('+', '*')], ...) 
    })
    return(list('+'=posCounts, '-'=negCounts))
}
#countStranded = countOverlapsInBinsStranded(gr, binList)
#combBinCovPos = sapply(countStranded[['+']], combFun)


#-----------------------------------------------------------------------
# Strand-specific coverage Profile along TAD boundaries
# retunrs a list for each strand with a list with a nbins long vector for each function in combFunctions
#-----------------------------------------------------------------------
strandedCountsAroundTADboundaries <- function(binList, grl, combFun=mean){
            
    # get strand-specific coverage for each bin
    binCovListStranded = lapply(grl, countOverlapsInBinsStranded, binList, type="any", ignore.strand=TRUE)
    message("INFO: Finished strand-specific overlap calculations.")

    # combine coverage for all i-th bins with the input functions separately for positive and negative strand
    combBinCovStranded = lapply(binCovListStranded, function(binCovStrandList){
        list(
            "+" = sapply(binCovStrandList[['+']], combFun),
            "-" = sapply(binCovStrandList[['-']], combFun)
            )
        })
        
    message("INFO: Finished combining of counts over all individual TADs.")
    
    return(combBinCovStranded)
}
#strandedCounts = strandedCountsAroundTADboundaries(binList[1:3], regionCountTrees)
