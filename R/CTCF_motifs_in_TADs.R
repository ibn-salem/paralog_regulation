########################################################################
# TODO: check +/- 1 error in output BED!
########################################################################

require(stringr)        # for some string functionality
require(RColorBrewer)   # for nice colors
require(plyr)           # count() function
require(gplots)         # heatmap.2 function
require(ggplot2)        # for nice plots
require(scales)         # for alpha() function for transparent colors

require(Biostrings)    # to work with genomic sequences
require(BSgenome.Hsapiens.UCSC.hg19)    # hg19 genome sequence
require(TFBSTools)     # for TF sequence motifs
require(JASPAR2014)     # for TF sequence motifs


# load some custom functions
source("R/functions.plot.R")
source("R/functions.Hi-C.R")
source("R/parseHiC.R")
source("R/chip-exo.functions.R")
source("R/functions.GRanges.R")

# load ENCODE datas
source("R/data.ensembl.R")
#source("R/data.encode.R")

########################################################################


#-----------------------------------------------------------------------
# Define some parameters
#-----------------------------------------------------------------------
#~ chromFile = "data/hg19/hg19.genome"
#~ realChromFile = "data/hg19/hg19.genome.realChroms"


RaoDomainFiles = c(
    Rao_HeLa="data/Rao2014/GSE63525_HeLa_Arrowhead_domainlist.txt",
    Rao_HUVEC="data/Rao2014/GSE63525_HUVEC_Arrowhead_domainlist.txt",
    Rao_K562="data/Rao2014/GSE63525_K562_Arrowhead_domainlist.txt",
    Rao_KBM7="data/Rao2014/GSE63525_KBM7_Arrowhead_domainlist.txt",
    Rao_NHEK="data/Rao2014/GSE63525_NHEK_Arrowhead_domainlist.txt",
    Rao_IMR90="data/Rao2014/GSE63525_IMR90_Arrowhead_domainlist.txt",
    Rao_GM12878="data/Rao2014/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt"
)

DixonDomainFiles = c(
    Dixon_hESC="data/Dixon2012/hESC.hg18.bed.hg19.bed",
    Dixon_IMR90="data/Dixon2012/IMR90.hg18.bed.hg19.bed"
)

CTCF_PEAK_DIR="data/ENCODE"
ENCODE_CTCF_METADATA_FILE="data/metadata.tsv.idr.CTCF.withHeader"


#ctcfPeakMotifFile = paste0("results/loop_prediction", "/MA0139.1/CTCF_HeLa_MA0139.1_pval0.0001_s100.bed")


VERSION="v02"
BOUNDARY_SIZE=2*10^4
outDataPrefix = "results/CTCF_motif_in_TAD_boundary/"
outPrefix = paste0("results/CTCF_motif_in_TAD_boundary/", VERSION)

# create directory, if not exist
dir.create(dirname(outPrefix), recursive=TRUE, showWarnings = FALSE)

WORKIMAGE_FILE = paste0(outPrefix, "workspace.Rdata")
#load(WORKIMAGE_FILE)

#-----------------------------------------------------------------------
# Parse TADs from Rao et al 2014 and Dixon et al 2012
#-----------------------------------------------------------------------
genome <- BSgenome.Hsapiens.UCSC.hg19

# parse TAD data sets as list of GRanges
RaoTADs = lapply(RaoDomainFiles, parseDomainsRao, disjoin=FALSE, seqinfo=seqinfo(genome))
DixonTADs <- lapply(DixonDomainFiles, import.bed, seqinfo=seqinfo(genome))

stableTADs <- list(
#~     "stable_TADs_n3_f80"=getConservedTADs(RaoTADs, n=3, fraction=.8),
#~     "stable_TADs_n3_f90"=getConservedTADs(RaoTADs, n=3, fraction=.9),
#~     "stable_TADs_n4_f80"=getConservedTADs(RaoTADs, n=4, fraction=.8),
    "stable_TADs"=getConservedTADs(RaoTADs, n=4, fraction=.9)
)

#~ stableTADsRes <- list(
#~     "stable_TADs_n3_10kb"=getConservedByHits(resolutionOverlap(allTADsRaoGR, allTADsRaoGR, resolution=10^4), n=3),
#~     "stable_TADs_n3_50kb"=getConservedByHits(resolutionOverlap(allTADsRaoGR, allTADsRaoGR, resolution=5*10^4), n=3),
#~     "stable_TADs_n4_10kb"=getConservedByHits(resolutionOverlap(allTADsRaoGR, allTADsRaoGR, resolution=10^4), n=4),
#~     "stable_TADs_n4_50kb"=getConservedByHits(resolutionOverlap(allTADsRaoGR, allTADsRaoGR, resolution=10^4), n=4)
#~ )

allTADs = c(
    RaoTADs,
    stableTADs,
    DixonTADs
)
message("INFO: Finshed parsing of TADs.")

# calculate bounaries between TADs
boundatryCoords = lapply(allTADs, getNonOverlappingBoundaries, min.gapwidth=1)

allBoundaries = lapply(boundatryCoords, resize, fix="center", BOUNDARY_SIZE)

disjoinedTADs <- lapply(allTADs, disjoin)
disjoinedBoundaries <- lapply(disjoinedTADs, getBoundariesFromDisjoinedTAD, maxSize=4*10^5)

#-----------------------------------------------------------------------
# Parse CTCF peaks form ENCODE
#-----------------------------------------------------------------------
ctcfMetaData <- read.table(ENCODE_CTCF_METADATA_FILE, sep="\t", header=TRUE)

# define extra columns in BED file from ENCODE (See: http://genome.ucsc.edu/FAQ/FAQformat.html#format12)
colClasses <- c("numeric", "integer", "numeric", "integer")
names(colClasses) <- c("signalValue", "pValue", "qValue", "peak")

# iterate over all BED files and pares them as GRange object
CTCFPeaksGRL <- lapply(ctcfMetaData$File.accession, function(ac){
    message(paste("INFO: parse CTCF peak:", ac))
    rtracklayer::import.bed(con=paste0(CTCF_PEAK_DIR, "/", ac, ".bed"), extraCols=colClasses, seqinfo=seqinfo(genome))
})

names(CTCFPeaksGRL) <- ctcfMetaData$File.accession
CTCFPeaksGRL <- GRangesList(CTCFPeaksGRL)

########################################################################

# get CTCF motif from JASPAR 
optCTCF <- list(species=9606, name="CTCF")
pfm <- getMatrixSet(JASPAR2014, optCTCF)[[1]]
pwm <- toPWM(getMatrixSet(JASPAR2014, optCTCF)[[1]])

pdf(paste0(outPrefix, ".CTCF_Motif_JASPAR.", ID(pfm), ".logo.pdf"), w=7, h=3.5)
    seqLogo(toICM(pfm))
dev.off()

# iterate over all chromosomes 
siteList <- bplapply(seqnames(genome), function(chr){

    chrSeq <- genome[[chr]]
    message(paste("INFO: Search motif on chromosome:", chr))
    
    searchSeq(pwm, chrSeq, seqname="CTCF", min.score="80%", strand="*")
})

# convert to GRanges
siteGRlist <- lapply(1:length(siteList), function(i){
    sites <- siteList[[i]]
    chr <- seqnames(genome)[i]
    
    gr <- GRanges(chr, 
        as(views(sites), "IRanges"), 
        strand=strand(sites),
        score=score(sites),
        relScore=relScore(sites), 
        seqinfo=seqinfo(genome)
        )
})

# convert to a single GR
siteAllGR <- unlist(GRangesList(siteGRlist))

#-----------------------------------------------------------------------
# get overlap with CTCF peaks
#-----------------------------------------------------------------------
siteAllGR$peaks <- countOverlaps(siteAllGR, CTCFPeaksGRL)
siteAllGR$peaks_percent <- 100 * siteAllGR$peaks / length(CTCFPeaksGRL)

for (PEAK_TH in c(0, 50, 75)){
    
    siteInPeakSubset <- siteAllGR[siteAllGR$peaks_percent >= PEAK_TH]

    for (MOTIF_TH in c(80, 90, 95)){
        
        siteSubset <- siteInPeakSubset[100*siteInPeakSubset$relScore >= MOTIF_TH]
        
        allBoundariesCTCF <- lapply(allBoundaries, function(bGR){
            subsetByOverlaps(siteSubset, bGR)
        })
    
        # write all TADs to BED files:
        for (tadName in names(allTADs)){
        
            # write CTCF motifs in all boundarys
            export(sort(allBoundariesCTCF[[tadName]]), paste0(outPrefix, ".TAD_Boundary_", BOUNDARY_SIZE, ".CTCFpeaks_", PEAK_TH, "MotifScore_", MOTIF_TH, ".", tadName, ".CTCF_motifs.bed"))
        
        }
        
        # combine all TAD boundaries
        allBoundariesCTCFallTADs <-reduce(sort(unlist(GRangesList(allBoundariesCTCF))))

        export(allBoundariesCTCFallTADs, paste0(outPrefix, ".TAD_Boundary_", BOUNDARY_SIZE, ".CTCFpeaks_", PEAK_TH, "MotifScore_", MOTIF_TH, ".ALL_TADs.CTCF_motifs.bed"))
        
    }
}


# write all TADs to BED files:
for (tadName in names(allTADs)){

    # wirte TAD
    export(sort(allTADs[[tadName]]), paste0(outPrefix, ".TAD_data.", tadName, ".bed"))

    # write all boundarys
    export(sort(allBoundaries[[tadName]]), paste0(outPrefix, ".TAD_Boundary_", BOUNDARY_SIZE, ".", tadName, ".bed"))

    # write all disjoined TAD boundarys
    export(sort(disjoinedBoundaries[[tadName]]), paste0(outPrefix, ".TAD_disjoined_Boundary.", tadName, ".bed"))
}

#=======================================================================
# save workspace image
#=======================================================================
save.image(WORKIMAGE_FILE)

#~ stop("INFO: Finished. Stopped script here!")

