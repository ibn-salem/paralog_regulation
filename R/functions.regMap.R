#!/usr/bin/Rscript
#=======================================================================
#
#   Functions related to regulatory interaction maps between enhancers 
#   and promoters
#
#=======================================================================

# load some useful libraries
require(stringr)
require(GenomicRanges)
require(BiocParallel)   # for parallel computing


#-----------------------------------------------------------------------
# read the regulatory interaction map from the format provied by the
# FANTOM5 consortium. It needs need a GRanges object with genes (or TSS), a GRanges object with all enhancers, a mappin of RefSeq IDs to ENSG id.
# It the regualtory map as data.frame of ids.
#-----------------------------------------------------------------------
parseFANTOM5RegMap <- function(regMapFile, tssGR, ehGR, refSeqToENSG){
    
    # filter refSeqToENSG to only those ENSG IDs in tssGR:
    refSeqToEnsgFiltered <- refSeqToENSG[refSeqToENSG %in% names(tssGR)]

    # read FANTOM5 regulatory map, which contains bed interval of associted enhancer to promoters
    FANTOM5 = read.delim(regMapFile)
    
    #chr1:67198280-67198800;NM_001037339;PDE4B;R:0.385;FDR:0
    regGene = as.data.frame(str_split_fixed(FANTOM5$name, ";", 5))
    names(regGene) = c("loc", "RefSeq", "Gene", "R", "FDR")
    regGene$R = as.numeric(sub("^R:", "", regGene$R))
    regGene$FDR = as.numeric(sub("^FDR:", "", regGene$FDR))
    refSeqList <- str_split(regGene$RefSeq, ",")
    
    # map each refSeq ID to a ENSGIDs
    ensgList <- bplapply(refSeqList, function(ids){
        ensgID <- refSeqToEnsgFiltered[ids]
        ensgID[!is.na(ensgID)]
        })
    
    # get unique ESNG for possible multiple (transcripts) refSeq IDs mapping to the same ENSG ID
    uniqueEnsgList <- lapply(ensgList, unique)
    
    # count IDs
    countEnsg <- sapply(uniqueEnsgList, length)
    
    # get indexes of only the uniquelly mapped associations
    uniqueMappable <- which(countEnsg == 1)
    
    # get unique mapped ENSG IDs
    ensgIDs <- unlist(uniqueEnsgList[uniqueMappable])
    
    message(paste0("INFO: FANTOM enhancer map: ", length(uniqueMappable), " of ",  nrow(FANTOM5), " enhancer-associated genes could be mapped uniquelly to ENSG IDs."))
    
    # make GR for enhancers in FANTOM file to get the index of ehGR
    usedEnhancerGR = GRanges(FANTOM5[uniqueMappable,1], IRanges(FANTOM5[uniqueMappable, "thickStart"], FANTOM5[uniqueMappable,"thickEnd"]), seqinfo=seqinfo(ehGR))
    indexInAllEnhancer = subjectHits(findOverlaps(usedEnhancerGR, ehGR))

    # get map as data.frame of IDs
    regMap = data.frame(

        # map ENSG to index in GRange object for genes
        tss=match(ensgIDs, names(tssGR)),

        # put index of associated enhancer in ehGR
        enhancer=indexInAllEnhancer
    )

    # add annotation fro FANTOM5 map
    regMap = cbind(regMap, regGene[uniqueMappable,])
    
    return(regMap)
}
    

#-----------------------------------------------------------------------
# Add distance to regulatory map data.frame
#-----------------------------------------------------------------------
getMapDist <- function(map, xGR, yGR, use.strand=TRUE){
    
    x_center =  mid(ranges(xGR[map[,1]]))
    y_center =  mid(ranges(yGR[map[,2]]))
    
    if (use.strand){
        strandFac = ifelse(strand(xGR[map[,1]]) == '+', 1, -1)
    }else{
        strandFac = 1
    }
    
    strandFac * (y_center - x_center)
    
}


#-----------------------------------------------------------------------
# get GRanges object for interacting pairs
# If strand.as.direction=TRUE the 'strand' will indicate the orientation
# of each pair on the chromosome, where  x >= y is '+' and x < y is '-'.
#-----------------------------------------------------------------------
getMapAsGR <- function(map, xGR, yGR, strand.as.direction=TRUE){
    
    # take the center coordinates of each range
    x_center =  mid(ranges(xGR[map[,1]]))
    y_center =  mid(ranges(yGR[map[,2]]))
    
    forward = x_center >= y_center
    

    if (strand.as.direction){
        
        # indicator variable that x < y
        strand = ifelse(forward, '+', '-')
        
    }else{
        strand = rep('*', nrow(map))
    }
    
    gr = GRanges(
        seqnames(xGR[map[,1]]),
        IRanges(
            ifelse(forward, y_center, x_center ),
            ifelse(forward, x_center, y_center )),
        strand = strand,
        seqinfo=seqinfo(xGR))
    
    # loop over additional columns and add 
    # add additional columns as meta data columns
    mcols(gr) = map[,-c(1,2)]
    names(mcols(gr)) <- names(map)[-c(1,2)]
    
    return(gr)
    
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
