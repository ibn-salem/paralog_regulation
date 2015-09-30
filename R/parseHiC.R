#!/usr/bin/Rscript
#=======================================================================
#
#   Provides functionality to parse Hi-C data 
#   from Rao et al 2014 and Capture Hi-C data from 
#   Mifsud et al. 2015
#
#=======================================================================

require(stringr)        # for paste0 function
require(gdata)
require(Matrix)         # for sparse matrix data structure
require(HiTC)           # for Hi-C experiment class and functions
require(rtracklayer)    # for import.bed
require(BiocParallel)   # for parallel computing

#-----------------------------------------------------------------------
# Creates a GRanges object with bins at a given resolution on given chrom
#-----------------------------------------------------------------------
getBinGR <- function(chr, resolution, seqInfo){
    options("scipen"=999)
    chrLen = seqlengths(seqInfo)[chr]
    starts = seq(1, chrLen, resolution)
#~     starts = seq(0, chrLen, resolution)
    ends = starts+resolution-1
    ends = ifelse(ends < chrLen, ends, chrLen)
    n = length(starts)
    gr = GRanges(rep(chr, n), 
        IRanges(starts, ends),
        names=as.character(starts-1),
        seqinfo=seqInfo
    )
    options("scipen"=0)
    return(gr)
}
# getBinGR(chr, resolution, seqInfo)


#-----------------------------------------------------------------------
# write GRange object in bed format to output file
#-----------------------------------------------------------------------
writeGRtoBED <- function(gr, outFile){
    
    options("scipen"=999)
    df <- data.frame(seqnames=seqnames(gr),
      starts=start(gr),
      ends=end(gr)-1,
      names=gr$names,
      scores=c(rep(".", length(gr))),
      strands=strand(gr))
    
    write.table(df, file=outFile, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
    options("scipen"=0)
    
}


#-----------------------------------------------------------------------
# Rewritten function to import Hi-C experiments as HiTC object without
# the need to providing bed files, but GRanges objects as bins instead
#-----------------------------------------------------------------------
importC_with_bins <- function (con, binGR, allPairwise = FALSE, 
    forceSymmetric = TRUE, tempFile) 
{
    # write bedfiles to temp files
    writeGRtoBED(binGR, tempFile)
    
    # call regular importC function:
    hiCexp = unlist(importC(con, tempFile, allPairwise=allPairwise, forceSymmetric=forceSymmetric))[[1]]
    
    # correct bins
    seqlevels(binGR) = as.character(seqnames(binGR)[1])
    x_intervals(hiCexp) = binGR
    y_intervals(hiCexp) = binGR

    return(hiCexp)

}

#-----------------------------------------------------------------------
#
#  Normalize interaction matrix using normalization vector.
#
#  Given the raw contact matrix M \in R^nxn and normalization vecotor
#  v \in R^n, the normalized counts M*_ij = M_ij / (v_i * v_j) or written
#  as matrix product M* = diag(v^-1) * M * diag(v^-1) 
#
#-----------------------------------------------------------------------
normalizeMatrix <- function(M, v, trancateVector=TRUE){
    
    # fix entries of normalization vector
    #TODO check why this happens
    if (trancateVector) {
        v = v[1:nrow(M)]
    }

    # remove NaN values and replace it by 1 (This should not influence the normalization)
    # TODO recheck 
    v[is.nan(v)] = 1

    # create a diagonal matrix of inverted normalization values
    Vinv = Diagonal(x=1/v)

    # normalize interaction matrix such as M*_ij = M_ij / v_i*v_j
    Mnorm = Vinv %*% M %*% Vinv

}


#-----------------------------------------------------------------------
# Normalize an interaction matrix by expected counts for each bin distances
#-----------------------------------------------------------------------
normalizeByDistanceExpected <- function(hiCexp, expectedFile){

    M = intdata(hiCexp)
    
    # read exp
    e = read.delim(expectedFile, header=FALSE, colClasses="numeric")[,1]
    
    # create matrix with expected contacts
    n = length(x_intervals(hiCexp))
    
    # for each bin distance i-j add expected counts e[i-j]
    E <- sapply(1:n, function(i, j) e[abs(i-j)+1], 1:n)
    
    # get observe/expected matrix
    intdata(hiCexp) = M / E
    
    return(hiCexp)

}
#~ normedHiCexp = normalizeByDistanceExpected(hiCexp, "data/Rao2014/IMR90/50kb_resolution_intrachromosomal/chr1/MAPQGE30/chr1_50kb.KRexpected")


#-----------------------------------------------------------------------
# parse and normalize Hi-C data for one experiment and chromosome 
#-----------------------------------------------------------------------
parseSingleHiCexp <- function(interactionFile, normalizationFile, chr, resolution, seqInfo, tempFile){

    # create bins for given resolution and chromosome size as GRanges object:
    binGR = getBinGR(chr, resolution, seqInfo)
    
    # write bin coordinates as temporary BED file 
    tempBinFile=paste0(tempFile, ".temp.res_bins.bed")
    
    # parse Hi-C data as HiTC experiment object
    hiCexp = importC_with_bins(interactionFile, binGR, forceSymmetric=TRUE, tempFile=tempBinFile)
    
    # parse normalization vector
    v = read.delim(normalizationFile, header=FALSE, colClasses="numeric")[,1]
    
    # normalize interaction counts    
    intdata(hiCexp) = normalizeMatrix(intdata(hiCexp), v)
    
    return(hiCexp)
}

#~ interactionFile="data/Rao2014/IMR90/50kb_resolution_intrachromosomal/chr22/MAPQGE30/chr22_50kb.RAWobserved"
#~ normalizationFile="data/Rao2014/IMR90/50kb_resolution_intrachromosomal/chr22/MAPQGE30/chr22_50kb.KRnorm"
#~ chr="chr22"
#~ resolution=50000
#~ tempFile=paste0(interactionFile, ".temp.bin.bed")


#-----------------------------------------------------------------------
# parse Hi-C map for all chromosomes (only intra-chrom) at given resolution 
#-----------------------------------------------------------------------
parseRaoHiC <- function(cell, resolution, dirPrefix, seqInfo, normalizeByExpected=FALSE, mapQualStr="MAPQGE30", normSuffix=".KRnorm", expectedSuffix=".KRexpected"){

    # get all the proper file paths
    # An example path is: 
    # IMR90/100kb_resolution_intrachromosomal/chr1/MAPQGE30/chr1_100kb.RAWobserved
    
    # map resolution in bp to kb/mb substring:
    res2str = c(
        "1000"="1kb",
        "5000"="5kb",
        "10000"="10kb",
        "25000"="25kb",
        "50000"="50kb",
        "100000"="100kb",
        "250000"="250kb",
        "500000"="500kb",
        "1000000"="1mb"
        )

    options("scipen"=999)
    resStr = res2str[as.character(resolution)]
    options("scipen"=0)
        
    # build path to directory with chromosome subdirectories
    chromDir = file.path(dirPrefix, cell, paste0(resStr, "_resolution_intrachromosomal"))

    # get all available chromosome names:
    chromosomes = list.dirs(path =chromDir , full.names = FALSE, recursive = FALSE)
    
    # iterate over each available chromosome to parse intra-chromosomal map
    hiClist = bplapply(chromosomes, function(chr){
        
        message(paste("Begin to parse data for chromosome", chr, "..."))
        rawMatrixFile = file.path(chromDir, chr, mapQualStr, paste0(chr, "_", resStr, ".RAWobserved"))
        
        normFile = file.path(chromDir, chr, mapQualStr, paste0(chr, "_", resStr, normSuffix))
        
        # parse and normalize intra-chromsomal contacts 
        exp = parseSingleHiCexp(rawMatrixFile, normFile, chr, resolution, seqInfo=seqInfo, tempFile=paste0(rawMatrixFile, ".temp.bin.bed"))
        
        # normalize by expected counts given each bin distance
        if (normalizeByExpected){
            expectedfile = file.path(chromDir, chr, mapQualStr, paste0(chr, "_", resStr, expectedSuffix))
            exp = normalizeByDistanceExpected(exp, expectedfile)
        }
        return(exp)
        
    })

    return(HTClist(hiClist))
    
}


#-----------------------------------------------------------------------
# parse and normalize Hi-C data for one experiment and chromosome 
#-----------------------------------------------------------------------
parseRudanHiC <- function(inFile, seqInfo, resolution=50*10^3){
        
    # parse input file as data frame
    classes <- c(rep(c("character", "integer", "integer"),2), "double", "double" )
    intDF <- read.table(inFile, header = TRUE, colClasses = classes)
    
    message("Finished reading of input file.")

    # initialize list variable with false (because empty list is not allowed)
    hiClist = list()
    hiClistNormed = list()
        
    # iterate over chromosomes
    for (chr in seqnames(seqInfo)){

        message(paste("Begin to parse data for chromosome", chr, "..."))
        
        # get subset of pairwise interaction on this chromosome
        subDF = intDF[intDF$chrom1 == sub("chr", "", chr),]
        subSeqInfo = seqInfo[chr]
        
        # make GenomicRanges object for bins at given resolution
        binGR = getBinGR(chr, resolution, subSeqInfo)            
        
        # initialize matrix
        M = Matrix(0, length(binGR), length(binGR))
        rownames(M) = binGR$names
        colnames(M) = binGR$names
        
        # copy empty matrix for normalized values
        Mnormed = M

        rowIdx = match(as.character(subDF[, "start1"]), binGR$names)
        colIdx = match(as.character(subDF[, "start2"]), binGR$names)

        contacts = subDF[, "observed_count"]
        expected = ifelse(subDF[, "expected_count"] != 0, subDF[, "expected_count"], NA)
        expected = subDF[, "expected_count"]
        
        normedContacts = log2( contacts / expected )
        # correct -Inf to NA
        normedContacts[is.infinite(normedContacts)] = NA
        
        # two column with indices to select cells in matrix:
        # (see: http://adv-r.had.co.nz/Subsetting.html)
        select = matrix(cbind(rowIdx, colIdx), ncol=2)
        
        # replace values with scores
        M[select] = contacts
        
        # add only those values that are not NA or -Inf (other will be still 0 which is wrong!) TODO fix this
        validSubset = !is.na(normedContacts) & !is.infinite(normedContacts)
        
        Mnormed[select[validSubset,]] = normedContacts[validSubset]
        

        # create HTCexp object    
        exp = HTCexp(M, binGR, binGR)
        expNormed = HTCexp(Mnormed, binGR, binGR)

        # append map to list (or create list if not exists)
        if (length(hiClist) > 0){
            hiClist = c(hiClist, exp )
            hiClistNormed = c(hiClistNormed, expNormed )
        }else{
            hiClist = HTClist(exp)
            hiClistNormed = HTClist(expNormed)
        }
    }
    return(list(hiClist, hiClistNormed))
}
#~ inFile = "data/Rudan2015/GSE65126_HiC_dog_liver_merged_50000.txt"
#~ seqInfo=seqInfoDog
#~ resolution=50*10^3


#-----------------------------------------------------------------------
# parse promoter-promoter interaction from Capture Hi-C (Mifsud2015a)
#-----------------------------------------------------------------------
parseCaptureHiC <- function(inFile, tssGR){

    # parse input file as data frame
    classes <- sapply(read.delim(inFile, nrows = 5, header=TRUE, stringsAsFactors=FALSE ), class)
    inData = read.delim(inFile, header=TRUE, colClasses = classes, stringsAsFactors=FALSE)
    
    # iterate over rows in indata and parse pairwise contacts
    # This might run 35 min
#~     pwList <- lapply(1:100, function(i) {
    pwList <- lapply(1:nrow(inData), function(i) {
        
        if (i %% 1000 == 0){
            message(paste("Parsing line", i, "..."))
        }
        
        genes1 = strsplit(inData[i, 5], split="[|]")[[1]]
        genes2 = strsplit(inData[i,11], split="[|]")[[1]]

        # filter for only those genes that are in tssGR 
        genes1 = genes1[genes1 %in% names(tssGR)]
        genes2 = genes2[genes2 %in% names(tssGR)]
        
        # get all possible pairs
#~         allPairs = rbind(expand.grid(genes1, genes2), expand.grid(gene2, genes1))
        allPairs = expand.grid(genes1, genes2)

        nPairs = nrow(allPairs)
        
        cbind(allPairs, rawCounts=rep(inData[i,13], nPairs), obsExp=rep(inData[i,14], nPairs))
                
    })
    
    message("Construct temorary data.frame ...")
    pairwiseDF = do.call("rbind", pwList)
    

    message("construct the interaction matrices ...")

    # initialize matrix
    n = length(tssGR)
    
    # two column with indices to select cells in matrix:
    # (see: http://adv-r.had.co.nz/Subsetting.html)
    select = apply(pairwiseDF[, 1:2], 2, match, names(tssGR))

#~     Mraw = sparseMatrix(i=select[,1], j=select[,2], x=pairwiseDF[,3], dims=c(n,n), dimnames=list(names(tssGR), names(tssGR)), symmetric = TRUE)
    Mraw = sparseMatrix(i=select[,1], j=select[,2], x=pairwiseDF[,3], dims=c(n,n), dimnames=list(names(tssGR), names(tssGR)))
    Mraw = forceSymmetric(symmpart(Mraw) * 2)

    MobsExp = sparseMatrix(i=select[,1], j=select[,2], x=pairwiseDF[,4], dims=c(n,n), dimnames=list(names(tssGR), names(tssGR)))
    MobsExp = forceSymmetric(symmpart(MobsExp) * 2)
    
    return(list("raw"=Mraw, "obsExp"=MobsExp))
}
#captureHiC <- parseCaptureHiC(inFile="data/Mifsud2015/TS5_GM12878_promoter-promoter_significant_interactions.txt", tssGR)

#-----------------------------------------------------------------------
# Query an Hi-C interaction matrix with two (or more) genomic intervals.
# E.g to get the interaction frequencies between two TSS coordinates.
# This function takes around 1 second for 10 queries.
#-----------------------------------------------------------------------
getInteractions <- function(xRange, yRange, HiClist){
    
    # assume only one interval in xRanges and yRange
    stopifnot(length(xRange) == 1 & length(yRange) == 1)
    
    # get chromosome names and check if map exists
    chrX = seqnames(xRange)
    chrY = seqnames(xRange)
    mapName = paste0(chrX, chrY)
    
    if ( mapName %in% names(HiClist) ){
        
        # get the X and Y interval indices in the contact matrix
        idxX = subjectHits(findOverlaps(xRange, x_intervals(HiClist[[mapName]])))
        idxY = subjectHits(findOverlaps(yRange, y_intervals(HiClist[[mapName]])))
        
        # TODO: handle input query ranges without overlap in map
        
        freq = intdata(HiClist[[mapName]])[idxX, idxY]
        
        return(freq)
    }
}

#-----------------------------------------------------------------------
# Query an Hi-C interaction matrix with two (or more) genomic intervals.
# E.g to get the interaction frequencies between two TSS coordinates.
#-----------------------------------------------------------------------
getInteractionsMulti <- function(xRange, yRange, HiClist, combineFun=sum){
    
    # assume equal length of the input ranges
    stopifnot(length(xRange) == length(yRange))
    stopifnot(as.character(seqnames(xRange)) == as.character(seqnames(yRange)))

    n = length(xRange)
    chroms = seqnames(xRange)
    
    # initilize zero vector
    freq = rep(NA, n)
    
    # iterate over all unique chromosomes (in parallel)
#~     for (chr in unique(chroms)){
    freqValues = bplapply(as.character(unique(chroms)), function(chr){
    
        message(paste("INFO: Query interactions for chromosome:", chr))

        # get indexes of input ranges on that chrom
        onChrom = which(chroms == chr)

        # get name of map in Hi-C data
        mapName = paste0(chr, chr)
        
        # check if interaction map is present, if not put NA
        if ( ! mapName %in% names(HiClist) ){
            
            message(paste("WARNING: Could not find map with name:", mapName))
#~             freq[onChrom] = NA
            return(NA)
        }else{
        
            map = HiClist[[mapName]]        
        
            # get indexes of bins in Hi-C map overlapping the query regions
            idxX = as.list( findOverlaps(xRange[onChrom], x_intervals(map)) ) 
            idxY = as.list( findOverlaps(yRange[onChrom], y_intervals(map)) ) 
    
            # query the interaction matrix with the indexes of bins
#~             freq[onChrom] = mapply(
            return(
                mapply(function(i,j){ combineFun( intdata(map)[i,j] )}, 
                    idxX, idxY)
                )
        }
        
        
    })

    # add chrom as names
    names(freqValues) <- unique(chroms)
    
    # now add values to data frame
    for (chr in unique(chroms)){
        # get indexes of input ranges on that chrom
        onChrom = which(chroms == chr)
        freq[onChrom] = freqValues[[chr]]
    }
    return(freq)
}
# xRange = tssGR[cisPairs[abs(cisPairs$dist) > 10^5,1][1:10]]
# yRange = tssGR[cisPairs[abs(cisPairs$dist) > 10^5,2][1:10]]
# freq = getInteractionsMulti(xRange, yRange, HiClist)

#-----------------------------------------------------------------------
# parse Rao Domains as GRange object.
# If disjoined is true, it returns a disjoined subset of non-overlapping domains
#-----------------------------------------------------------------------
parseDomainsRao <- function(inFile, disjoin=TRUE, ...){

    # parse IMR90 domains from Rao et al 2014:
    raoDF = read.delim(inFile)
    chr = paste0("chr", raoDF$chr1)

    up = raoDF[,"x1"]
#~     down = raoDF[,"x2"]
    # substract 1 bp from down coordinate to have inclusive interval ranges
    down = raoDF[,"x2"] - 1
    
    raoDomains = GRanges(chr, IRanges(up, down), ...)
    
    if(disjoin){
        raoDomains =  disjoin(raoDomains)
        raoDomains[width(raoDomains)>1,]
    }

    return(raoDomains)
}

#-----------------------------------------------------------------------
# parse loops from Rao et al 2015 as two GRange object.
#-----------------------------------------------------------------------
parseLoopsRao <- function(inFile, ...){

    # parse IMR90 domains from Rao et al 2014:
    raoDF = read.delim(inFile)
    chr = paste0("chr", raoDF$chr1)

    # substract 1 bp from down coordinate to have inclusive interval ranges    
    upAnchor = GRanges(chr, IRanges(raoDF[,"x1"], raoDF[,"x2"]-1), ...)
    downAnchor = GRanges(chr, IRanges(raoDF[,"y1"], raoDF[,"y2"]-1), ...)
    
    # define and add annotations:
    annotUP = c("o", "e_bl",  "e_donut", "e_h", "e_v", "fdr_bl", "fdr_donut", "fdr_h", "fdr_v", "num_collapsed", "centroid1", "centroid2", "radius", "motif_x1", "motif_x2", "sequence_1", "orientation_1", "uniqueness_1")
    annotDOWN = c("motif_y1", "motif_y2", "sequence2", "orientation_2", "uniqueness_2")

    mcols(upAnchor) = raoDF[,annotUP]
    mcols(downAnchor) = raoDF[,annotDOWN]

    return(list(upAnchor, downAnchor))
}

#-----------------------------------------------------------------------
# parse TADs from excel sheed provided as sup. table by Rudan et al. 2015
#-----------------------------------------------------------------------
parseRudanTADs <- function(inFile, disjoin=TRUE, sheet=1, ...){

    # use gdata package to read .xlsx file (See: http://www.r-bloggers.com/read-excel-files-from-r/)

    df = read.xls(inFile, sheet=sheet)
    
    # build GRanges object
    gr = GRanges(df[,1], IRanges(df[,2], df[,3]), ...)
    
    if(disjoin){
        gr =  disjoin(gr)
        gr[width(gr)>1,]
    }
    
    return(gr)

}


#-----------------------------------------------------------------------
# get the TAD boundareis as GRanges object from input TADs
# They are here defined as regions between TADs of some maximal size.
#-----------------------------------------------------------------------
getBoundariesFromDisjoinedTAD <- function(tadGR, maxSize=4*10^5){
    
    # assume non-overlaping input TADs
    stopifnot(isDisjoint(tadGR))

    #boundaries <- gaps(tadGR)
    boundaries <- gaps(resize(tadGR, width=width(tadGR)-1, ignore.strand=TRUE))
    boundaries <- boundaries[width(boundaries) <= maxSize]
    return(boundaries)
}

#-----------------------------------------------------------------------
# get the TAD boundareis as GRanges object from input TADs
# They are here defined as regions between TADs of some maximal size.
#-----------------------------------------------------------------------
getNonOverlappingBoundaries <- function(tadGR, min.gapwidth=1){
    
    # get regions around start and end pints of TADs
    borders = c(
        resize(tadGR, width=1, fix="start"),
        resize(tadGR, width=1, fix="end")
    )
    
    # remove overlapping (or nearby) boundaries
    mergedBoundaries = reduce(borders, min.gapwidth=min.gapwidth)
    
    resize(mergedBoundaries, width=1, fix="center")
    
}

#-----------------------------------------------------------------------
# get conserved TAD boundaries
#-----------------------------------------------------------------------
conservedBoundaries <- function(boundaryList, min.gapwidth=1){
    
    # make GRange object of all boundaries in all cell types
    allBoundaries = unlist(GRangesList(boundaryList))
    
    # reduce the set to non-overlapping ranges
    redBoundaries = reduce(allBoundaries, min.gapwidth=min.gapwidth)
    
    # count for each boundary the number of data sets wich supports (overlap) it
    counts = countOverlaps(redBoundaries, GRangesList(boundaryList))
    
    # add the counts as "conservation" column to the set of all reduced boundaries
    redBoundaries$conservation = counts
    
    return(redBoundaries)
}

#-----------------------------------------------------------------------
# get conserved TADs
#-----------------------------------------------------------------------
getConservedTADs <- function(TADList, maxgap=10^4, n=3){
    
    # make GRange object of all TADs in all cell types
    allTADs = unlist(GRangesList(TADList))
        
    # count for each TAD the number TADs from all data sets that overlap it
    counts = countOverlaps(allTADs, allTADs, type="equal", maxgap=maxgap)
    
    # filter for those TADs that have >= n overlaps with other TADs
    consTADs = allTADs[counts >= n]
    
    # reduce TADs to not have overlapping duplicates
    consTADs = reduce(consTADs, min.gapwidth=maxgap+1)

    return(consTADs)
}

