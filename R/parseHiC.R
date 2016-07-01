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
# rewritten function from HiTC package (See HiTC:::ImportC)
#-----------------------------------------------------------------------
myImportC <- function (con, binGR) {
    
    # assign intervals 
    xgi <- binGR
    ygi <- binGR
    
    # parse interaction counts
    message(paste("INFO: Sart reading file:", con))
    cdata <- read.table(con, comment.char = "#", colClasses = c("character", 
        "character", "numeric"), check.names = FALSE)
    stopifnot(ncol(cdata) == 3)
    message("INFO: Finished reading file.")
    
    # map coordinates to IDs in binGR
    id1 <- cdata[, 1]
    id2 <- cdata[, 2]
    pos1 <- match(id1, id(ygi))
    pos2 <- match(id2, id(xgi))

    # build a sparseMatrix (with zeros for missing data9
    bigMat <- Matrix::sparseMatrix(i = pos1, j = pos2, x = cdata[, 
        3], dims = c(length(ygi), length(xgi)), dimnames = list(id(ygi), 
        id(xgi)))
    
    
    stopifnot(length(unique(seqnames(binGR))) == 1)
    chr <- as.character(unique(seqnames(ygi)))
    
    seqlevels(ygi) <- chr
    seqlevels(xgi) <- chr

    message(paste("INFO: Convert matrix into HTCexp object(s) for chr:", chr))
    hiCexp <- HTCexp(bigMat, xgi, ygi, lazyload = FALSE)   
    
    # delete input data
    rm(cdata)
    rm(bigMat)

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
# parse Hi-C map for all chromosomes (only intra-chrom) at given resolution 
#-----------------------------------------------------------------------
parseRaoHiC <- function(cell, resolution, dirPrefix, seqInfo, normalizeByExpected=FALSE, mapQualStr="MAPQGE30", normStr="KR"){

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
    hiClist = lapply(chromosomes, function(chr){
        
        message(paste("INFO: Begin to parse data for chromosome", chr, "..."))
        rawMatrixFile = file.path(chromDir, chr, mapQualStr, paste0(chr, "_", resStr, ".RAWobserved"))

        # create bins for given resolution and chromosome size as GRanges object:
        binGR = getBinGR(chr, resolution, seqInfo)

        # parse intra-chromsomal contacts 
        hiCexp <- myImportC(rawMatrixFile, binGR)
        
        if (!is.null(normStr)){

            normFile = file.path(chromDir, chr, mapQualStr, paste0(chr, "_", resStr, ".", normStr, "norm"))
            expectedfile = file.path(chromDir, chr, mapQualStr, paste0(chr, "_", resStr, ".", normStr, "expected"))
                
            # if KR normalization vector file is empty, the normalization did not converge
            # Rao et al suggest to take the VC or SQRTVC normalization in this case
            if (file.size(normFile) == 0 & normStr=="KR"){
                normFile = file.path(chromDir, chr, mapQualStr, paste0(chr, "_", resStr, ".VCnorm"))
                expectedfile = file.path(chromDir, chr, mapQualStr, paste0(chr, "_", resStr, ".VCexpected"))
            }

            # parse normalization vector
            v = read.delim(normFile, header=FALSE, colClasses="numeric")[,1]

            # normalize interaction counts    
            message(paste("INFO: Begin normalization..."))
            intdata(hiCexp) = normalizeMatrix(intdata(hiCexp), v)
            
            # normalize by expected counts given each bin distance
            if (normalizeByExpected){
                message(paste("INFO: Begin normalization by expected counts..."))
                hiCexp = normalizeByDistanceExpected(hiCexp, expectedfile)
            }
            message(paste("INFO: Finished normalization."))
        }

        return(hiCexp)
        
    })
    
    Sys.sleep(1) # hack to fix problems with bplapply on MOGON


    return(HTClist(hiClist))
    
}
#~ rawMatrixFile
#~ cell=CELL
#~ resolution=HIC_RESOLUTION
#~ dirPrefix=HIC_DATA_DIR
#~ normalizeByExpected=FALSE
#~ mapQualStr="MAPQGE30"
#~ normStr="KR"
#~ expectedSuffix=".KRexpected"

#-----------------------------------------------------------------------
# parse expected counts by distance
#-----------------------------------------------------------------------
parseRaoExpected <- function(cell, resolution, dirPrefix, mapQualStr="MAPQGE30", normStr="KR"){


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

    expectedList = lapply(chromosomes, function(chr){

        expectedfile = file.path(chromDir, chr, mapQualStr, paste0(chr, "_", resStr, ".", normStr, "expected"))

        # if KR normalization vector file is empty, the normalization did not converge
        # Rao et al suggest to take the VC or SQRTVC normalization in this case
        if (file.size(expectedfile) == 0 & normStr=="KR"){
            expectedfile = file.path(chromDir, chr, mapQualStr, paste0(chr, "_", resStr, ".VCexpected"))
        }

        # parse expected counts by distance
        e <- read.delim(expectedfile, header=FALSE, colClasses="numeric")[,1]
    
    })
    names(expectedList) <- chromosomes

    return(expectedList)
}


#-----------------------------------------------------------------------
# parse and normalize Hi-C data for one experiment and chromosome 
#-----------------------------------------------------------------------
parseRudanHiC <- function(inFile, seqInfo, resolution=50*10^3){
        
    # parse input file as data frame
    classes <- c(rep(c("character", "integer", "integer"),2), "double", "double" )
    intDF <- read.table(inFile, header = TRUE, colClasses = classes)
    
    message("INFO: Finished reading of input file.")

    # initialize list variable with false (because empty list is not allowed)
    hiClist = list()
    hiClistNormed = list()
        
    # iterate over chromosomes
    for (chr in seqnames(seqInfo)){

        message(paste("INFO: Begin to parse data for chromosome", chr, "..."))
        
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
# This function sums up all pairwise interactions counts in case of duplicated ENSG pair entries in the input file
#-----------------------------------------------------------------------
parseCaptureHiC <- function(inFile, tssGR){

    # parse input file as data frame
    classes <- sapply(read.delim(inFile, nrows = 5, header=TRUE, stringsAsFactors=FALSE ), class)
    inData <- read.delim(inFile, header=TRUE, colClasses = classes, stringsAsFactors=FALSE)
    
    # iterate over rows in indata and parse pairwise contacts
    # This runs in less than 9 min
    pwList <- lapply(1:nrow(inData), function(i) {
        
        if (i %% 1000 == 0){
            message(paste("INFO: Parsing line", i, "..."))
        }
        
        genes1 <- strsplit(inData[i, 5], split="[|]")[[1]]
        genes2 <- strsplit(inData[i, 11], split="[|]")[[1]]

        # filter for only those genes that are in tssGR 
        genes1 <- genes1[genes1 %in% names(tssGR)]
        genes2 <- genes2[genes2 %in% names(tssGR)]
        
        # get all possible pairs between the two gene sets
        allPairs <- expand.grid(genes1, genes2)

        nPairs <- nrow(allPairs)
        
        # add counts raw.count and log.observed.expected
        cbind(allPairs, rawCounts=rep(inData[i,13], nPairs), obsExp=rep(inData[i,14], nPairs))
                
    })
    
    message("INFO: Construct temorary data.frame ...")
    pairwiseDF = do.call("rbind", pwList)
    

    message("INFO: construct the interaction matrices ...")

    # initialize matrix
    n = length(tssGR)
    
    # two column with indices to select cells in matrix:
    # (see: http://adv-r.had.co.nz/Subsetting.html)
    select = apply(pairwiseDF[, 1:2], 2, match, names(tssGR))

    # create sparse matrix 
    Mraw = sparseMatrix(i=select[,1], j=select[,2], x=pairwiseDF[,3], dims=c(n,n), dimnames=list(names(tssGR), names(tssGR)))
    
    # make matrix symmetric by using the 2 times the average of upper and lower triangular matrix
    Mraw = forceSymmetric(symmpart(Mraw) * 2)

    # create sparse matrix for obs/exp
    MobsExp = sparseMatrix(i=select[,1], j=select[,2], x=pairwiseDF[,4], dims=c(n,n), dimnames=list(names(tssGR), names(tssGR)))
    MobsExp = forceSymmetric(symmpart(MobsExp) * 2)
    
    return(list("raw"=Mraw, "obsExp"=MobsExp))
}
#captureHiC <- parseCaptureHiC(inFile="data/Mifsud2015/TS5_GM12878_promoter-promoter_significant_interactions.txt", tssGR)

#-----------------------------------------------------------------------
# Query an Hi-C interaction matrix with two (or more) genomic intervals of width 1.
# E.g to get the interaction frequencies between two TSS coordinates.
#-----------------------------------------------------------------------
getInteractionsPoint <- function(xRange, yRange, HiClist, inParallel=TRUE, ignoreSameBin=FALSE){
    
    # assume equal length of the input ranges and same chromosome
    stopifnot(length(xRange) == length(yRange))
    stopifnot(as.character(seqnames(xRange)) == as.character(seqnames(yRange)))
    
    # assume intervals of length 1
    stopifnot(width(xRange) == 1 & width(yRange) == 1)

    n = length(xRange)
    chroms = seqnames(xRange)
    
    # initilize zero vector
    freq = rep(NA, n)
    
    # iterate over all unique chromosomes (in parallel)
    freqValues = get(ifelse(inParallel, "bplapply", "lapply"))(as.character(unique(chroms)), function(chr){
    
        message(paste("INFO: Query interactions for chromosome:", chr))

        # get indexes of input ranges on that chrom
        onChrom = which(chroms == chr)

        # get name of map in Hi-C data
        mapName = paste0(chr, chr)
        
        # check if interaction map is present, if not put NA
        if ( ! mapName %in% names(HiClist) ){
            
            message(paste("WARNING: Could not find map with name:", mapName))
            return(NA)
        }else{
        
            map = HiClist[[mapName]]        
        
            # get indexes of bins in Hi-C map overlapping the query regions
            idxX = subjectHits( findOverlaps(xRange[onChrom], x_intervals(map)) ) 
            idxY = subjectHits( findOverlaps(yRange[onChrom], y_intervals(map)) )
    
            # query the interaction matrix with the indexes of bins
            contacts <- intdata(map)[cbind(idxX,idxY)]
            
            # ignore contact counts if both regions map to the same bin
            if (ignoreSameBin){
                
                sameBin <- idxX == idxY
                
                contacts[sameBin] <- NA
            }
                        
            return(contacts)
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


#-----------------------------------------------------------------------
# Query an Hi-C interaction matrix with two (or more) genomic intervals.
#-----------------------------------------------------------------------
getInteractionsMulti <- function(xRange, yRange, HiClist, combineFun=mean, inParallel=TRUE, ignoreSameBin=FALSE){
    
    # assume equal length of the input ranges
    stopifnot(length(xRange) == length(yRange))
    stopifnot(as.character(seqnames(xRange)) == as.character(seqnames(yRange)))

    n = length(xRange)
    chroms = seqnames(xRange)
    
    # initilize zero vector
    freq = rep(NA, n)
    
    # iterate over all unique chromosomes (in parallel)
    freqValues = get(ifelse(inParallel, "bplapply", "lapply"))(as.character(unique(chroms)), function(chr){
    
        message(paste("INFO: Query interactions for chromosome:", chr))

        # get indexes of input ranges on that chrom
        onChrom = which(chroms == chr)

        # get name of map in Hi-C data
        mapName = paste0(chr, chr)
        
        # check if interaction map is present, if not put NA
        if ( ! mapName %in% names(HiClist) ){
            
            message(paste("WARNING: Could not find map with name:", mapName))
            return(NA)
        }else{
        
            map = HiClist[[mapName]]        
        
            # get indexes of bins in Hi-C map overlapping the query regions
            idxX = as.list( findOverlaps(xRange[onChrom], x_intervals(map)) ) 
            idxY = as.list( findOverlaps(yRange[onChrom], y_intervals(map)) ) 
    
            # query the interaction matrix with the indexes of bins
            contacts <- mapply(function(i,j){ 
                        combineFun( intdata(map)[i,j] )
                    }, 
                    idxX, idxY)
            
            # ignore contact counts if both regions map to the same bin
            if (ignoreSameBin){
                
                sameBin <- mapply(function(i,j){
                    length(i) == length(j) & all(i == j)
                }, idxX, idxY)
                
                contacts[sameBin] <- NA
            }
                        
            return(contacts)
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
#~ xRange = tssGR[gP[,1]]
#~ yRange = tssGR[gP[,2]]
# freq = getInteractionsMulti(xRange, yRange, HiClist)

#-----------------------------------------------------------------------
# parse Rao Domains as GRange object.
# If disjoined is true, it returns a disjoined subset of non-overlapping domains
# Assumes coordinate in zero-based half open format.
#-----------------------------------------------------------------------
parseDomainsRao <- function(inFile, disjoin=TRUE, ...){

    # parse IMR90 domains from Rao et al 2014:
    raoDF = read.delim(inFile)
    chr = paste0("chr", raoDF$chr1)

    # add 1 bp to start coordinate to convert to one-based inclusive interval format used in GenomicRanges
    up = raoDF[,"x1"] + 1
    down = raoDF[,"x2"] 
    
    raoDomains = GRanges(chr, IRanges(up, down), ...)
    
    if(disjoin){
        raoDomains =  disjoin(raoDomains)
        raoDomains[width(raoDomains)>1,]
    }

    return(raoDomains)
}

#-----------------------------------------------------------------------
# parse TADs from excel sheed provided as sup. table by Rudan et al. 2015
#-----------------------------------------------------------------------
parseRudanTADs <- function(inFile, disjoin=TRUE, sheet=1, ...){

    # use gdata package to read .xlsx file (See: http://www.r-bloggers.com/read-excel-files-from-r/)

    df = read.xls(inFile, sheet=sheet)
    
    # build GRanges object
    # add 1 bp to start coordinate to convert to one-based inclusive 
    gr = GRanges(df[,1], IRanges(df[,2]+1, df[,3]), ...)
    
    if(disjoin){
        gr =  disjoin(gr)
        gr[width(gr)>1,]
    }
    
    return(gr)

}
