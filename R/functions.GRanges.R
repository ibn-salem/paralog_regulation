#!/usr/bin/Rscript
#=======================================================================
#
#   Functions related to GRanges objects from the GenomicRanges package
#
#=======================================================================

# load some useful libraries
require(stringr)
require(GenomicRanges)
require(igraph)

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
