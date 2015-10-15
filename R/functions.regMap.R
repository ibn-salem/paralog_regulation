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


#-----------------------------------------------------------------------
# read the regulatory interaction map from the format provied by O'Conner et al.
#-----------------------------------------------------------------------
parseOConnorRegMap <- function(regMapFile){

    # read interaction map
    regMap = read.delim(regMapFile, header=FALSE)
    
    # split the first column by '_'
    ENSEMBL_TSS = str_split_fixed(regMap[,1], "_", 2)
    
    regMap[,"ENSEMBL"] = ENSEMBL_TSS[,1]
    regMap[,"TSS"] = as.numeric(ENSEMBL_TSS[,2])
    
    # split the second column by '_'
    CELL_CRM = str_split_fixed(regMap[,2], "_", 2)
    
    regMap[,"cell"] = CELL_CRM[,1]
    regMap[,"CRM_ID"] = as.numeric(CELL_CRM[,2])
    
    # delete initial and now redundant regMap columns
    regMap = regMap[,-c(1,2)]
}

#-----------------------------------------------------------------------
# read the regulatory interaction map from the format provied by the
# FANTOM5 consortium. It needs need a data.frame with all genes from ESNEMBL (including HGNC symbol, transcription_start, strand...) as additional input.
# It returns a list which includes the regualtory map as data.frame of ids, a Grange object with with only the used TSS and a GRanges object with the used enhancers.
#-----------------------------------------------------------------------
parseFANTOM5RegMap <- function(regMapFile, allGenes, allEnhancers, seqInfo){

    # read FANTOM5 regulatory map, which contains bed interval of associted enhancer to promoters
    FANTOM5 = read.delim(regMapFile)
    
    #chr1:67198280-67198800;NM_001037339;PDE4B;R:0.385;FDR:0
    regGene = as.data.frame(str_split_fixed(FANTOM5$name, ";", 5))
    names(regGene) = c("loc", "RefSeq", "Gene", "R", "FDR")
    regGene$R = as.numeric(sub("^R:", "", regGene$R))
    regGene$FDR = as.numeric(sub("^FDR:", "", regGene$FDR))
        
    # make GR for enhancers with unique IDs
    enhancerGR = GRanges(FANTOM5[,1], IRanges(FANTOM5$thickStart, FANTOM5$thickEnd), seqinfo=seqInfo)
    indexInAllEnhancer = subjectHits(findOverlaps(enhancerGR, allEnhancers))
    
    # unique genes by HGNC symbol and remove "" as symbol
    allGenesHGNC = allGenes[!duplicated(allGenes$hgnc_symbol),]
    allGenesHGNC = allGenesHGNC[allGenesHGNC$hgnc_symbol != "",]

    # make tssGR with only the genes used in the map
    tssGR = getTssGRfromENSEMBLGenes(allGenesHGNC[allGenesHGNC$hgnc_symbol %in% regGene$Gene, ], seqInfo)
    
    # get map as data.frame of IDs
    map = data.frame(

        # map HGNC symbol to index in GRange object for genes
        tss=match(regGene$Gene, tssGR$hgnc_symbol),

        # get index in ehGR for enhancerID
        enhancer=indexInAllEnhancer
    )
    
    # add annotation fro FANTOM5 map
    map = cbind(map, regGene)

    # remove unmappable genes
    map = map[!is.na(map[,1]),]
    
    return(list(map=map, tssGR=tssGR))
    
}
#fantomData = parseFANTOM5RegMap(REGMAP_FILE_FANTOM5, allGenes, seqInfo)
#mapFANTOM=fantomData[["map"]]
#tssGRFANTOM=fantomData[["tssGR"]]
#ehGRFANTOM=fantomData[["ehGR"]]

#-----------------------------------------------------------------------
# Get the locus string chr:start-end form a GRanges object
#-----------------------------------------------------------------------
getLocStr <- function(gr, substractOneFromStart=FALSE){
    if(!substractOneFromStart){
        paste0(seqnames(gr), ':', start(gr), '-', end(gr))
    }else{
        paste0(seqnames(gr), ':', start(gr)-1, '-', end(gr))    
    }
}

#-----------------------------------------------------------------------
# Writes a long-range interaction genome browser track file in the format
# of the EpiGenome Broweser. See: http://wiki.wubrowse.org/Long-range
# TODO: add sort, bgzip, tabix -bed commands!
#-----------------------------------------------------------------------
writeInteractionTrack <- function(map, gr1, gr2, outFile, scoreColumn="R"){
    
    n = nrow(map)
    reg1 = gr1[map[,1]]
    reg2 = gr2[map[,2]]
    
    tab1 = cbind(as.character(seqnames(reg1)), start(reg1)-1, end(reg1),
     paste0(getLocStr(reg2, TRUE), ",", map[,scoreColumn]), 1:n, ".")
    
    # make reverse interactions (required by track fromat)
    tab2 = cbind(as.character(seqnames(reg2)), start(reg2)-1, end(reg2),
     paste0(getLocStr(reg1, TRUE), ",", map[,scoreColumn]), (n+1):(2*n), ".")
    
    # write table to output file
    write.table(rbind(tab1, tab2), outFile, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

}

#-----------------------------------------------------------------------
# Writes a long-range interaction genome browser text file in the text file format of the EpiGenome Broweser. See: http://wiki.wubrowse.org/Long-range
#-----------------------------------------------------------------------
writeInteractionTextFile <- function(map, gr1, gr2, outFile, scoreColumn="R"){
    
    n = nrow(map)
    reg1 = gr1[map[,1]]
    reg2 = gr2[map[,2]]
    
    tab = cbind(getLocStr(reg1, TRUE), getLocStr(reg2, TRUE), map[,scoreColumn])
    
    # write table to output file
    write.table(tab, outFile, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

}

#-----------------------------------------------------------------------
# Parse CRM from BED file (from Yip et al 2012) as GRanges object
#-----------------------------------------------------------------------
parseCRMfromYip <- function(bedFile, seqInfo){
    crmDF = read.delim(bedFile, header=FALSE)
    names(crmDF) <- c("chr", "start", "end", "type")
    #crmDF$ID = 1:nrow(crmDF)
    
    # convert CRMs into an GenomicRanges object
    CRMs = makeGRangesFromDataFrame(crmDF, seqinfo=seqInfo, keep.extra.columns=TRUE)
    names(CRMs) = 1:nrow(crmDF)
    return(CRMs)
}
#~ # get the midpoint of CRMs to associate it to a unique genomic bin
#~ CRMcenter = resize(CRMs, width=1, fix="center")
#~ CRMcenterTree = GIntervalTree(CRMcenter)


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
