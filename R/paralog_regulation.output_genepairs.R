

PARAM_SCRIPT="R/paralog_regulation.param.R"

message(paste("INFO: Load parameter from this script:", PARAM_SCRIPT))
source(PARAM_SCRIPT)

# load some other functions from this project
source("R/functions.genePairs.R")

#=======================================================================
# get mapping of ENSG to Entrez gene IDs
#=======================================================================
require(biomaRt)        # to retrieve human paralogs from Ensembl
ensemblGRCh37 <- useMart(host="grch37.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", verbose=FALSE)

geneAttributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene")
geneIDs = getBM(attributes=geneAttributes, mart=ensemblGRCh37)


#=======================================================================
# load file with all pairs
#=======================================================================

#~ # save as table
#~ write.table(allDF, file=paste0(outPrefix, ".allDF.csv"),
#~     sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)


#~ allDF <- read.table(paste0(outPrefix, ".allDF.csv"), sep="\t", 
#~ 	header=TRUE)

# load aP and allDF
load(file=paste0(WORKIMAGE_FILE, ".aP_allDF.Rdata"))


# filter for paralogs and sampled genes with features:
# 	- abs(dist)
#	- same TAD 
#	- group


# define some often used subsetting vectors for the allDF
sSamp <- allDF$sampType %in% c("paralogs", "sampled by distance")
sTAD <- allDF$tadSource == "Dixon_hESC"
sExp <- allDF$expSource == "GTEx"
sSpec <- allDF$species == "mouse"
sDist <- allDF$distGroup == "close"

subDF <- subset(allDF, sSamp & sTAD & sExp & sSpec & sDist)

subDF <- subDF[,c("g1", "g2", "group", "replicate", "inTAD", "dist")]

# map ENSG to entrezgene

subDF$g1_entrez <- geneIDs$entrezgene[match(subDF$g1, geneIDs$ensembl_gene_id)]
subDF$g2_entrez <- geneIDs$entrezgene[match(subDF$g2, geneIDs$ensembl_gene_id)]

subDF$g1_hgnc <- geneIDs$hgnc_symbol[match(subDF$g1, geneIDs$ensembl_gene_id)]
subDF$g2_hgnc <- geneIDs$hgnc_symbol[match(subDF$g2, geneIDs$ensembl_gene_id)]


write.table(subDF, file=paste0(outPrefix, ".subDF_pairs_Dixon_hESC_TADs.entrez_and_hgnc.csv"),
    sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)


#=======================================================================
# check distance distribution of paralog and sampled genes
#=======================================================================
subDF <- read.table(paste0(outPrefix, ".subDF_pairs_Dixon_hESC_TADs.entrez_and_hgnc.csv"), header=TRUE, sep="\t", quote="")
ddply(subDF, .(group), summarize, mean=mean(dist), median=median(dist), n=length(dist))

ppiDF <- read.table("data/temp/TAD_h2_symbol_v2.tsv", header=TRUE, sep="\t", quote="")
ddply(ppiDF, .(group), summarize, mean=mean(dist), median=median(dist), n=length(dist))



#=======================================================================
# Analyse genomic distance distribution of PPI and non PPI gene pairs
#=======================================================================

#-------------------------------------------------------------------
# get tssGR for ENSG
#-------------------------------------------------------------------
require(biomaRt)        # to retrieve human paralogs from Ensembl
require(BSgenome.Hsapiens.UCSC.hg19)
seqInfo <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)

ensemblGRCh37 <- useMart(host="grch37.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", verbose=FALSE)

geneAttributes = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position", "strand", "status", "gene_biotype")
geneFilters="chromosome_name"
# read "normal" human chromosome names (without fixes and patches)
geneValues=c(1:22, "X", "Y")
allGenes = getBM(attributes=geneAttributes, mart=ensemblGRCh37, filters=geneFilters, values=geneValues)

# unique gene entry by ENSG ID symbol:
#genes = knownCodingGenes[!duplicated(knownCodingGenes$ensembl_gene_id),]
genes = allGenes[!duplicated(allGenes$ensembl_gene_id),]
#~ rownames(genes) <- genes$ensembl_gene_id

# make GRanges object for all known prot coding genes

tssGR = GRanges(
        paste0("chr", genes$chromosome_name),
        IRanges(genes$start_position, genes$start_position),
        strand = ifelse(genes$strand == 1, '+', '-'), 
        names = genes$ensembl_gene_id, 
        genes[,c("hgnc_symbol", "status", "gene_biotype")],
        seqinfo=seqInfo
        )
names(tssGR) = genes$ensembl_gene_id
tssGR <- sort(tssGR)



#-----------------------------------------------------------------------
# add linear distance between genes (by assuming them on the same chromosome)
#-----------------------------------------------------------------------
addPairDistKb <- function(genePairs, tssGR){
    # get chromosomes of gene pairs
    sameChrom = as.character(seqnames(tssGR[genePairs[,1]])) == as.character(seqnames(tssGR[genePairs[,2]]))
    s1 = start(tssGR[genePairs[,1]])
    s2 = start(tssGR[genePairs[,2]])
    # add a new column "dist" to the data.frame
    genePairs[, "dist"] = ifelse(sameChrom, abs(s2-s1)/1000, NA)
    return(genePairs)
}

#-----------------------------------------------------------------------
# Parse HIPPIE 
#-----------------------------------------------------------------------
SCORE_TH <- 0.72


# load mapping of entrez IDs to ENSEMBL IDs
VERSION_DATA_ENSEMBL="v01"
outPrefixDataEnsembl=paste0("results/data/ensembl.data.", VERSION_DATA_ENSEMBL)
load(paste0(outPrefixDataEnsembl, ".ENSEMBL_GRCh37.entrezToENSG.RData"))

# parse hippie and map IDs to ENSG
hippieURL <- "http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/hippie_current.txt"

#~ hippieDF <- read.table("data/HIPPIE/hippie_current.txt", header=FALSE, sep="\t")
hippieDF <- read.table(hippieURL, header=FALSE, sep="\t")
hippie <- data.frame(
	g1_ENGS = entrezToENSG[as.character(hippieDF[,2])],
	g2_ENSG = entrezToENSG[as.character(hippieDF[,4])],
	symbol1 = str_split_fixed(as.character(hippieDF[,1]), "_", 2)[,1],
	symbol2 = str_split_fixed(as.character(hippieDF[,3]), "_", 2)[,1],
	score = hippieDF[,5],
	stringsAsFactors=FALSE)

message("INFO: After parsing: ", nrow(hippie))

# filter out interactions that could not be mapped to ENSG
hippie <- hippie[!is.na(hippie[,1]) & !is.na(hippie[,2]),]
message("INFO: After ENSG mapping: ", nrow(hippie))

# filter out interaction bellow score threshold
hippie <- hippie[hippie$score >= SCORE_TH,]


# generate random interaction network
randNet <- hippie[rep(1:nrow(hippie), 10) ,c("g1_ENGS", "g2_ENSG", "score")]
randNet[,2] <- sample(randNet[,2])
#~ randNet <- addPairDistKb(randNet, tssGR)



pairsDF <- rbind(
	hippie[,c("g1_ENGS", "g2_ENSG", "score")],
	randNet
	)
	
pairsDF$group <- rep(c("PPI", "shuffled"), c(nrow(hippie), nrow(randNet)))

message("INFO: After filtering score >= ", SCORE_TH, " : ", sum(pairsDF$group == "PPI"), " and shuffled: ",sum(pairsDF$group == "shuffled"))




pairsDF <- addPairDistKb(pairsDF, tssGR)

pairsDF <- pairsDF[!is.na(pairsDF$dist),]

message("INFO: After filtering out different chromosomes : ", sum(pairsDF$group == "PPI"), " and shuffled: ",sum(pairsDF$group == "shuffled"))


pairsDF <- pairsDF[pairsDF$dist >= 10 & pairsDF$dist <= 1000,]
message("INFO: After filtering distance within 10 and 1000kb: ", sum(pairsDF$group == "PPI"), " and shuffled: ",sum(pairsDF$group == "shuffled"))


#-----------------------------------------------------------------------
# annotate to be in same TAD
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



tadName <- "Rao_IMR90"
tadGR <- import(paste0(outPrefix, ".TAD_data.", tadName, ".bed"), seqinfo=seqInfo)
#boundaryGR <- getBoundaries(tadGR)


pairsGR <- getPairAsGR(pairsDF, tssGR)

pairsDF$inTAD <- countOverlaps(pairsGR, tadGR, type="within") >= 1
pairsDF$inTAD <- factor(pairsDF$inTAD, c(TRUE, FALSE), c("Same TAD", "Not same TAD"))

#-----------------------------------------------------------------------
# plot distance
#-----------------------------------------------------------------------



require(ggplot2)
pVal <- wilcox.test(dist ~ group, data=pairsDF)$p.value


p <- ggplot(pairsDF, aes(dist, ..density.., fill=group)) + 
#~ 	geom_density(adjust = 1/4, alpha=.5) + 
	geom_histogram(binwidth=50, alpha=.5, position="identity") + 
	labs(title=paste("p =", signif(pVal, 3)), x="Genomic distance [kb]")  + 
	theme_bw()	 
ggsave(p, file=paste0(outPrefix, ".hippie_genomic_distance.v02.hist.pdf"), w=7, h=3.5)

p <- p + facet_grid(inTAD~., margins=TRUE, scales="free_y")
ggsave(p, file=paste0(outPrefix, ".hippie_genomic_distance.v02.hist.byTAD.pdf"), w=7, h=7)



#~ p <- ggplot(pairsDF, aes(y=dist, x=group, color=group)) + 
#~ 	geom_boxplot() + 
#~ 	labs(title=paste("p =", signif(pVal, 3)), x="Genomic distance [kb]")  + 
#~ 	theme_bw()


