########################################################################
#
# A script to analyse cohesin KO w.r.t the co-regulation of paralog genes.
#
########################################################################

require(stringr)        # for some string functionality
require(GenomicRanges)  # for genomic intervals and overlap calculation
require(rtracklayer)    # for import.bed
require(BiocParallel)   # for parallel computing
require(reshape2)		# for melt() and cast()
require(tidyr)			# to reformat data.frame's

#-----------------------------------------------------------------------
# Load parameters from external script
#-----------------------------------------------------------------------
# read parameter script location from command-line argument

#~ args <- commandArgs(trailingOnly = TRUE)
#~ PARAM_SCRIPT=args[1]
PARAM_SCRIPT="R/paralog_regulation.param.R"

message(paste("INFO: Load parameter from this script:", PARAM_SCRIPT))
print(PARAM_SCRIPT)

source(PARAM_SCRIPT)

# TADs from Rao et al:
raoMousTADfile = "data/Rao2014/GSE63525_CH12-LX_Arrowhead_domainlist.txt.bed.mm9.bed.mm10.bed"


DixonMouseDomainFiles = c(
    Dixon_mESC="data/Dixon2012/mouse.mESC.mm9.bed.mm10.bed",
    Dixon_cortex="data/Dixon2012/mouse.cortex.mm9.bed.mm10.bed"
)



#-----------------------------------------------------------------------
# Options for parallel computation
#-----------------------------------------------------------------------
# use all available cores but generate random number streams on each worker
multicorParam <- MulticoreParam(RNGseed=RANDOM_SEED)   
# set options
register(multicorParam)  
# bpparam() # to print current options


#-----------------------------------------------------------------------
# load some custom functions
#-----------------------------------------------------------------------
source("R/functions.GRanges.R")
source("R/functions.genePairs.R")
source("R/functions.genePairs.randomization.R")
source("R/functions.genePairs.paralog_analysis.R")
source("R/parseHiC.R")

# load ensemble data sets of genes and paralog pairs
source("R/data.ensembl.R")     # load ensambl data
source("R/data.cohesinKO.R")  # load cohesin KO expression data
# combine expression data sets into one list:

#=======================================================================
# load some other data
#=======================================================================
#-----------------------------------------------------------------------
# Parse Hi-C data for other species
#-----------------------------------------------------------------------

# parse mouse and dog TADs
speciesSeqInfo = list("mmusculus"=seqInfoMouse, "cfamiliaris"=seqInfoDog)


speciesTADs = lapply(1:2, function(i) parseRudanTADs(RudanFile, sheet=i, disjoin=FALSE, seqinfo=speciesSeqInfo[[i]]))
names(speciesTADs) = names(speciesSeqInfo)

speciesHiC = list(
	"mmusculus" = parseRudanHiC("data/Rudan2015/GSE65126_HiC_mouse_liver_merged_50000.txt", seqInfoMouse, resolution=50*10^3)
)

#=======================================================================
# annotate gene pairs
#=======================================================================
#~ runBasicParalogAnalysis <- function(outPrefix, paralogPairs, pairScoreCol, tssGR, genesGR, TAD, tissueName, HiClist, HiClistNorm){

outPrefix=paste0(outPrefix, ".MouseCohesinKO")
paralogPairs=paralogPairsMouse
pairScoreCol="mmusculus_paralog_ds"
tssGR=speciesTssGR[["mmusculus"]]
genesGR=speciesGenesGR[["mmusculus"]]
TAD=speciesTADs[["mmusculus"]]
tissueName="Mouse"
HiClist=speciesHiC[["mmusculus"]][[1]]
HiClistNorm=speciesHiC[["mmusculus"]][[2]]
seqInfo <- speciesSeqInfo[["mmusculus"]]

# parse TADs from Rao et al in mouse CH12-LX cells
RaoTADs = import(raoMousTADfile, seqinfo=seqInfo)

DixonTADs <- lapply(DixonMouseDomainFiles, import.bed, seqinfo=seqInfo)


allTADs = list(
	"VietriRudan_liver"=speciesTADs[["mmusculus"]],
	"Rao_CH12-LX"=RaoTADs
)
allTADs <- c(allTADs, DixonTADs)


# get boundaries
allBoundaries <- lapply(allTADs, getBoundaries)

LOAD_PAIRS=FALSE
if (!LOAD_PAIRS) {


	#-------------------------------------------------------------------
	# Filtering:
	#-------------------------------------------------------------------
	message("Start to filter paralog pairs...")

	# remove double entries of the form A-B and B-A
	paralogPairsUniqP = uniquePair(paralogPairs)

	# filter out overlapping gene pairs
	nonOVL <- nonOverlappingGenePairs(paralogPairsUniqP, genesGR)
	paralogPairsUniqPnonOVL <- paralogPairsUniqP[nonOVL,]

	# write all paralog pairs to output file:
	write.table(paralogPairsUniqPnonOVL, file=paste0(outPrefix, ".paralog_pairs.paralogPairsUniqPnonOVL.txt"),
		sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

	# get for each gene only one unique pair, the one with highest similarity
	# this is computed by an maximum weight matching
	paralogPairsUniqPnonOVLWithDS = paralogPairsUniqPnonOVL[!is.na(paralogPairsUniqPnonOVL[,pairScoreCol]),]
	
	pairSelectFact = ifelse(SELECT_OLD_PAIRS, 1, -1)
	paralogPairsUniqPnonOVLUniqG = uniquePairPerGeneBySim(paralogPairsUniqPnonOVLWithDS, pairSelectFact*paralogPairsUniqPnonOVLWithDS[,pairScoreCol])
	
	# get only a unique pair order (one of A-B, B-A) form the unique pairs
	paralogPairsUniq = uniquePair(paralogPairsUniqPnonOVLUniqG)
	
	# subset of paralog pairs that are located on the same chromosome
	allCisPairs = getCisPairs(paralogPairsUniq, tssGR)
	
	#-----------------------------------------------------------------------
	# annotate gene pairs
	#-----------------------------------------------------------------------
	# add same Strand info to all paralogs
	paralogPairsUniq = addSameStrand(paralogPairsUniq, tssGR)
	
	# add linear distance between TSS
	paralogPairsUniq = addPairDist(paralogPairsUniq, tssGR)
	allCisPairs = addPairDist(allCisPairs, tssGR)
	
	# add HGNC symbols
	allCisPairs = addHGNC(allCisPairs, tssGR)
	
	# add same Strand info:
	allCisPairs = addSameStrand(allCisPairs, tssGR)
	
	# make GRanges objects for cis paralog pairs and random paris on same chromosome
	allCisPairsGR = getPairAsGR(allCisPairs, tssGR)
	
	# iterate over different TAD data sets
	for (tadName in names(allTADs)){
		
		# to check whether there are in the same TAD
		allCisPairsGR = addWithinSubject(allCisPairsGR,  allTADs[[tadName]], paste0("TAD_", tadName))
		
		# separated by a TAD boundary?
		separated <- countOverlaps(allCisPairsGR, allBoundaries[[tadName]]) >= 1
		mcols(allCisPairsGR)[, paste0("Boundary_", tadName)] <- separated
		 
	}


	# assign annotation in GRanges object to gene pair data.frames
	colNames <- c(paste0("TAD_", names(allTADs)), paste0("Boundary_", names(allTADs)))
	allCisPairs[,colNames] <- data.frame( mcols(allCisPairsGR)[,colNames] )

	
	# Adds Hi-C contact frequencies to a gene pair data set
	allCisPairs = addHiCfreq(allCisPairs, tssGR, HiClist)
	allCisPairs = addHiCfreq(allCisPairs, tssGR, HiClistNorm, label="HiCnorm")

    # add pairwise correlations of gene expression over all tissues
    for (expName in names(expCoDFlist)) {
        
        message(paste("INFO: annotate pairs with expression correlation form:", expName))
        expDF = expCoDFlist[[expName]]
        
        allCisPairs = addCor(allCisPairs, expDF, colName=paste0(expName, "_expCor"))
    }


	# save allCisPairs with all annotations
	write.table(allCisPairs, file=paste0(outPrefix, ".paralog_pairs.allCisPairs.annotated.txt"),
		sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
	
	#-----------------------------------------------------------------------
	# Separate pairs by distance
	#-----------------------------------------------------------------------

	# get close cis pairs
	closePairs = allCisPairs[abs(allCisPairs$dist) <= MAX_DIST,] 
	
	# get distal pairs
	distalPairs = allCisPairs[abs(allCisPairs$dist) > MAX_DIST,] 

	#-----------------------------------------------------------------------
	# paralog pair filtering numbers
	#-----------------------------------------------------------------------
	nPairs = c(
		"paralogPairs"=nrow(paralogPairs), 
		"paralogPairsUniqP"=nrow(paralogPairsUniqP), 
		"paralogPairsUniqPnonOVL"=nrow(paralogPairsUniqPnonOVL), 
		"paralogPairsUniqPnonOVLWithDS"=nrow(paralogPairsUniqPnonOVLWithDS), 
		"paralogPairsUniqPnonOVLUniqG"=nrow(paralogPairsUniqPnonOVLUniqG), 
		"paralogPairsUniq"=nrow(paralogPairsUniq),
		"allCisPairs"=nrow(allCisPairs),
		"closePairs"=nrow(closePairs),
		"distalPairs"=nrow(distalPairs)
		)

	write.table(nPairs, file=paste0(outPrefix, ".paralog_pairs_filtering.txt"),
		sep="\t", quote=FALSE, col.names=FALSE)

	message("INFO: Start sampling of gene pairs...")
	#===================================================================
	# Sample gene pairs
	#===================================================================
	
	#-----------------------------------------------------------------------
	# Sample pairs with equal probability from all genes
	#-----------------------------------------------------------------------
	randPairs <- lapply(1:N_RAND, function(x){getRandomPairs(nrow(paralogPairsUniq), names(tssGR))})
	Sys.sleep(3) # hack to fix problems with bplapply on MOGON
	
	# add same strand information
	randPairs <- lapply(randPairs, addSameStrand, tssGR)
	randPairs <- lapply(randPairs, addPairDist, tssGR)

	# filter for random pairs in Cis
	randPairsInCis <- lapply(randPairs, getCisPairs, tssGR)
	
	#-----------------------------------------------------------------------
	# Sample all cis pairs by distance using 90 bins 
	#-----------------------------------------------------------------------
	# Now sample from all possible gene pairs within 1 - DISTAL_MAX_DIST bp
	allCisGenePairsOVL <- getAllGenePairs(tssGR, maxDist=DISTAL_MAX_DIST, minDist=1)
	
	# filter out overlapping pairs
	allNonOVL <- nonOverlappingGenePairs(allCisGenePairsOVL, genesGR, useIDs=TRUE)
	allCisGenePairs <- allCisGenePairsOVL[allNonOVL,]

	# get sampling weights for distal pairs according to distance
	cisDistBreaks <- seq(log10(1), log10(DISTAL_MAX_DIST), length.out=91)
	cisDistProp <- weightsByBin(log10(abs(allCisPairs$dist)), log10(abs(allCisGenePairs$dist)), breaks=cisDistBreaks)
	
	# sample according to distance
	sampCisPairs <- bplapply(1:N_RAND, function(x){ 
		sampleFromAllPairsByWeight(n=nrow(allCisPairs), hitDF=allCisGenePairs, tssGR, weight=cisDistProp)
		})
	Sys.sleep(3) # hack to fix problems with bplapply on MOGON

	#===================================================================
	# Annotate sampled gene pairs
	#===================================================================

	# make GRanges objects for cis paralog pairs and random paris on same chromosome
	sampCisPairsGR = lapply(sampCisPairs, getPairAsGR, tssGR)
	Sys.sleep(3) # hack to fix problems with bplapply on MOGON


	#sampCisPairsGR = bplapply(sampCisPairsGR, addWithinSubject, TAD, "TAD")
	#Sys.sleep(3) # hack to fix problems with bplapply on MOGON

	for (tadName in names(allTADs)){
				
		sampCisPairsGR = bplapply(sampCisPairsGR, addWithinSubject, allTADs[[tadName]], paste0("TAD_", tadName))
		Sys.sleep(3) # hack to fix problems with bplapply on MOGON
		
		# separated by a TAD boundary?
		sampCisPairsGR <- bplapply(sampCisPairsGR, function(gr){
	
				separated <- countOverlaps(gr, allBoundaries[[tadName]]) >= 1
				mcols(gr)[, paste0("Boundary_", tadName)] <- separated
				
				return(gr)
			})
		Sys.sleep(3) # hack to fix problems with bplapply on MOGON		 
	}

	# assign annotation in GRanges object to gene pair data.frames
	colNames <- c(paste0("TAD_", names(allTADs)), paste0("Boundary_", names(allTADs)))
	for (i in seq(N_RAND)){
		sampCisPairs[[i]][,colNames] <- data.frame( mcols(sampCisPairsGR[[i]])[,colNames] )
	}


	#===================================================================
	# Annotate sampled gene pairs
	#===================================================================
	
	# combine all sampling replicates to one data frame
	randPairsCombined <- do.call("rbind", randPairs)
	randPairsInCisCombined <- do.call("rbind", randPairsInCis)
	sampCisPairsCombined <- do.call("rbind", sampCisPairs)

	# add Hi-C contact frequencies
	sampCisPairsCombined = addHiCfreq(sampCisPairsCombined, tssGR, HiClist)
	sampCisPairsCombined = addHiCfreq(sampCisPairsCombined, tssGR, HiClistNorm, label="HiCnorm")
        
	# add expression correlation
	for (expName in names(expCoDFlist)) {
		
		message(paste("INFO: Annotate sampled pairs with expression form:", expName))
		expDF <- expCoDFlist[[expName]]
	
		sampCisPairsCombined <- addCor(sampCisPairsCombined, expDF, colName=paste0(expName, "_expCor"))
	}
	
	#-----------------------------------------------------------------------
	# save a work image after sampling and annotation.
	#-----------------------------------------------------------------------
	save(
		paralogPairsUniq, 
		allCisPairs, 
		closePairs, 
		distalPairs,
		randPairs,  randPairsCombined,
		randPairsInCis, randPairsInCisCombined,
		sampCisPairs, sampCisPairsCombined,
		file=paste0(WORKIMAGE_FILE, ".annotated_gene_pairs.Rdata")
	)


    #-----------------------------------------------------------------------
    # build data.frame with all annotations and all pairs
    #-----------------------------------------------------------------------
    breaksCis = 10^(0:9) / 10^3
    breaksTAD = c(0, 10, 1000, Inf)

    nPairs <- nrow(allCisPairs)
    allCisPairs$g1 <- allCisPairs[,1]
    allCisPairs$g2 <- allCisPairs[,2]
    
    allPairDFs <- list(allCisPairs, sampCisPairsCombined)
    names(allPairDFs) <- c("paralogs", "sampled by distance")
    
    # get common columns in all data sets
    cC <- Reduce(intersect, lapply(allPairDFs, names))
    
    # combine all pairs by using only the common columns
    aP <- do.call("rbind", lapply(allPairDFs, function(df) df[,cC]))
    
    # add noZero HI-C 
    aP <- addNoZero(aP, c("HiC", "HiCnorm"))
    
    
    # add further columns
    aP[,"group"] <- factor(rep(c("paralog", "sampled"), c(nrow(allPairDFs[[1]]), sum(sapply(allPairDFs[2:length(allPairDFs)], nrow)))), levels=c("paralog", "sampled"))
    aP[,"sampType"] <- factor(rep(names(allPairDFs), times=sapply(allPairDFs, nrow)), levels=names(allPairDFs))
    aP[,"replicate"] <- c(rep(1, nPairs), replicate(length(allPairDFs)-1, rep(1:N_RAND, each=nPairs)))
    aP[,"dist"] <- abs(aP[,"dist"]) / 10^3
    aP[,"distGroup"] <- factor(ifelse(aP[,"dist"] <= (MAX_DIST/10^3), "close", "distal"), levels=c("close", "distal"))
    aP[,"distBin"] <- as.factor(breaksCis[.bincode(aP[,"dist"], breaksCis)])
    aP[,"distTadBin"] <- factor(breaksTAD[.bincode(aP[,"dist"], breaksTAD)], levels=breaksTAD[1:3], labels=c("<10kb", "10-1000kb", ">1000kb"))

    # save bright allDF data set with column for each source:
    write.table(aP, file=paste0(outPrefix, ".allPairs_broad.csv"),
        sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
    save(aP,file=paste0(WORKIMAGE_FILE, ".aP.Rdata"))


}else{
	message(paste("INFO: Start loading annotated gene pairs form this file:",  paste0(WORKIMAGE_FILE, ".annotated_gene_pairs.Rdata")))
	load(paste0(WORKIMAGE_FILE, ".annotated_gene_pairs.Rdata"))
	load(paste0(WORKIMAGE_FILE, ".aP.Rdata"))
	message("INFO: Finished loading.")
}

#-------------------------------------------------------------------
# make an allDF that repeat paris for different expression sources
#-------------------------------------------------------------------

breaksExpCor <- c(-1, -.5, 0, .5, 1)

aP$ID <- 1:nrow(aP)

# expand TAD to wide format
colIDsTAD <- match(c(paste0("TAD_", names(allTADs)), paste0("Boundary_", names(allTADs))), names(aP))

tmpDF <- aP %>%
	gather(v, value, colIDsTAD) %>%
	separate(v, c("type", "tad.study", "tad.tissue"), sep="_") %>%
	arrange(ID) %>%
	spread(type, value) 

# expand Expression to wide format
colIDsExp <- match(paste0(names(expCoDFlist), "_expCor"), names(tmpDF))
	
allDF <- tmpDF	%>%
	gather(v, expCor, colIDsExp) %>%
	separate(v, c("expType", "expCondition", "expTrash"), c(3,5))


# to make factors with meaning
allDF$TAD <- factor(allDF$TAD, c(TRUE, FALSE), c("Same TAD", "Not same TAD")) 
allDF$Boundary <- factor(allDF$Boundary, c(FALSE, TRUE), c("Not separated", "Cross boundary")) 
allDF$tad.source <- paste0(allDF$tad.study, "_", allDF$tad.tissue)
allDF$expTrash <- NULL
allDF$expSource <- paste0(allDF$expType, "_", allDF$expCondition)
allDF$eypCorBin <- factor(breaksExpCor[.bincode(allDF$expCor, breaksExpCor)], labels=c("high neg", "low neg", "low pos", "high pos"))
allDF$expCondition <- factor(allDF$expCondition, c("WT", "KO"))

# save bright allDF data set with column for each source:
write.table(allDF, file=paste0(outPrefix, ".allDF.csv"),
	sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
save(allDF,file=paste0(WORKIMAGE_FILE, ".allDF.Rdata"))


#~ # create data.frame with multiple rows per gene pair for the different data sources
#~ allDF <- cbind(
#~     aPflat[rep(1:nrow(aP), nExp),],
#~     data.frame(
#~         expSource = factor(rep(names(expCoDFlist), each=nrow(aP)), levels=names(expCoDFlist)),
#~         expCor = unlist(aP[,paste0(names(expCoDFlist), "_expCor")]),
#~         expCorBin=factor(breaksExpCor[.bincode( unlist(aP[,paste0(names(expCoDFlist), "_expCor")]), breaksExpCor)], labels=c("high neg", "low neg", "low pos", "high pos"))
#~     )
#~ )


########################################################################
# Analyse allDF
########################################################################

sExpRaw <- allDF$expType == "exp"
sExpRep <- allDF$expType == "rep"
sTadVit <- allDF$tad.source == "VietriRudan_liver"
sTadRao <- allDF$tad.source == "Rao_CH12-LX"

#~ sExpRep <- allDF$expSource %in% c("rep_WT", "rep_KO")

subList <- list(
	sExpRaw=sExpRaw,
	sExpRep=sExpRep
)

#-----------------------------------------------------------------------
# expression correlation in WT and KO by TAD and group
#-----------------------------------------------------------------------
for (expSubName in names(subList)){
	
	expDataSub <- subList[[expSubName]]
	#message(class(expDataSub))
	

#~ 	subDF <- subset(allDF, expDataSub & allDF$dist <= 10^3)
	subDF <- subset(allDF, expDataSub)
	
	for (tadGroup in c("TAD", "Boundary")){
	
		# calculate p-values
		pvalDF <- ddply(subDF, c("tad.source", "group", tadGroup), summarize, 
			p=wilcox.test(expCor ~ expCondition)$p.value)
				
		#~ nDF <- data.frame(table(subDF[, c("group", "expSource", "TAD")], useNA="ifany"))
		nDF <- ddply(subDF, c("tad.source", "group", tadGroup, "expCondition"), summarize, 
			n=sum(!is.na(expCor)),
			nAll=length(expCor),
			avg=mean(expCor, na.rm=TRUE),
			sd=sd(expCor, na.rm=TRUE)
			)
		
		p = ggplot(subDF, aes(x=expCondition, y=expCor)) + 
		    geom_boxplot(aes(colour = group), lwd=1.5) + 
		    facet_grid(as.formula(paste("tad.source ~ group +", tadGroup))) + 
		    scale_color_manual(values=COL, guide_legend(title = "")) +
		    theme_bw() + theme(text = element_text(size=15), axis.text.x=element_text(angle = 45, hjust = 1), legend.position="bottom") + 
		    guides(fill=guide_legend(title="")) +
		    labs(y="Expression Correlation [R]", x="") + 
			geom_text(aes(label=paste0("p=", signif(p,2)), x=1.5, y=1.1), data=pvalDF, size=5) +
		    geom_text(aes(label=paste0("n=",n), y=-1), data=nDF) +
		    geom_text(aes(label=paste0("N=",nAll), y=-1.2), data=nDF)
		    
		ggsave(p, file=paste0(outPrefix, ".expCor_WTvsKO.", expSubName, ".by_", tadGroup, "_byGroup.boxplot.pdf"), w=7, h=14)


	}

	#---------------------------------------------------------------
	# plot expression correlation by distance
	#---------------------------------------------------------------
#~ 	subDF = subset(subDF, tad.source == "Dixon_mESC" & allDF$dist <= 10^3)
	subDF = subset(subDF, tad.source == "Dixon_mESC")

	ALPHA=.1
	COL_COMB=brewer.pal(12, "Paired")[c(5,6,7,8,1,2,3,4)]

	p <- ggplot(subDF, aes(x=dist, y=expCor, fill=interaction(expCondition,TAD, group), col=interaction(expCondition,TAD,group))) + #interaction(expCondition,group))) + #, color=zvar, fill=zvar)) + 
	  geom_point(alpha=ALPHA, size=.5, shape=20)  +
	  geom_smooth(alpha=ALPHA) + 
	  guides(col=guide_legend(title=""), fill=guide_legend(title="")) + 
	  scale_color_manual(values = COL_COMB) + scale_fill_manual(values = COL_COMB) + 
	  theme_bw() + 
	  theme(legend.position="bottom") + 
	  labs(x="Distance between gene pairs [kb]", y="Expression Correlation [R]")+
	  scale_x_log10()
	ggsave(p, file=paste0(outPrefix, ".expCor_WTvsKO.", expSubName, ".by_Dixon_mESC_TAD_byGroup_byDist.fit.pdf"), w=8, h=8)

	p <- ggplot(subDF, aes(x=dist, y=expCor, fill=interaction(expCondition, Boundary, group), col=interaction(expCondition,Boundary,group))) + #interaction(expCondition,group))) + #, color=zvar, fill=zvar)) + 
	  geom_point(alpha=ALPHA, size=.5, shape=20)  +
	  geom_smooth(alpha=ALPHA) + 
	  guides(col=guide_legend(title=""), fill=guide_legend(title="")) + 
	  scale_color_manual(values = COL_COMB) + scale_fill_manual(values = COL_COMB) + 
	  theme_bw() + 
	  theme(legend.position="bottom") + 
	  labs(x="Distance between gene pairs [kb]", y="Expression Correlation [R]")+
	  scale_x_log10()
	ggsave(p, file=paste0(outPrefix, ".expCor_WTvsKO.", expSubName, ".by_Dixon_mESC_Boundary_byGroup_byDist.fit.pdf"), w=8, h=8)


}

#-----------------------------------------------------------------------
# pairs in same TAD
#-----------------------------------------------------------------------

subDF <- subset(allDF, allDF$expType == "exp" & allDF$expCondition == "WT" & distGroup=="close")
#~ subDF <- subset(allDF, allDF$expType == "exp" & allDF$expCondition == "WT")

pvalDF <- ddply(subDF, .(tad.source), summarize, 
	p=fisher.test(table(TAD, group))$p.value)

nDF <- ddply(subDF, .(tad.source, group), summarize, 
			n=sum(TAD == "Same TAD"),
			N=length(TAD),
			percent=100*n/N
			)			
p <- ggplot(nDF, aes(x=group, y=percent)) + 
	geom_bar(aes(fill=group), stat="identity", col="black") + 
	facet_grid(.~tad.source) +
	scale_fill_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=15), axis.text.x=element_blank(), legend.position="bottom") + 
    labs(y="Gene pairs in same TAD [%]", x="") + 
	geom_text(aes(label=paste0("p=", signif(p,2)), x=1.5, y=90), data=pvalDF, size=5) +
    geom_text(aes(label=paste0(round(percent,1),"%" )), vjust=1.5, size=4, data=nDF) + 
    geom_text(aes(label=paste0("n=",n)), vjust=3, size=4, data=nDF)

ggsave(p, file=paste0(outPrefix, ".paralogPairs.sameTAD_by_TADsource.barplot.pdf"), w=7, h=7)
		
		

#-----------------------------------------------------------------------
# Paris not separated by TAD boundary
#-----------------------------------------------------------------------
pvalDF <- ddply(subDF, .(tad.source), summarize, 
	p=fisher.test(table(Boundary, group))$p.value)

nDF <- ddply(subDF, .(tad.source, group), summarize, 
			n=sum(Boundary == "Not separated"),
			N=length(Boundary),
			percent=100*n/N
			)
p <- ggplot(nDF, aes(x=group, y=percent)) + 
	geom_bar(aes(fill=group), stat="identity", col="black") + 
	facet_grid(.~tad.source) +
	scale_fill_manual(values=COL, guide_legend(title = "")) +
    theme_bw() + theme(text = element_text(size=15), axis.text.x=element_blank(), legend.position="bottom") + 
    labs(y="Not sparated by TAD boundary [%]", x="") + 
	geom_text(aes(label=paste0("p=", signif(p,2)), x=1.5, y=100), data=pvalDF, size=5) +
    geom_text(aes(label=paste0(round(percent,1),"%" )), vjust=1.5, size=4, data=nDF) + 
    geom_text(aes(label=paste0("n=",n)), vjust=3, size=4, data=nDF)

ggsave(p, file=paste0(outPrefix, ".paralogPairs.Not_separated_by_TADsource.barplot.pdf"), w=7, h=7)

#===================================================================
# Run basic analysis
#===================================================================
message("Start with basic plots...")

#-----------------------------------------------------------------------
# plot number of genes in each group
#-----------------------------------------------------------------------
paralogs = unique(unlist(paralogPairs[,1:2]))
nonParalogs = setdiff(names(tssGR), paralogs)
nrGenes = sapply(list(paralogs, nonParalogs), length)    

pdf(paste0(outPrefix, ".number_of_genes.pdf"), width=3.5)    
	par(cex=1.5, mgp=c(3,1,0), mar = c(7, 4, 4, 2) + 0.1)
	barP = barplot(nrGenes, beside=TRUE, col=COL2, names.arg=NA, ylab="Number of genes", ylim=c(0, 1.25*max(nrGenes)), main=tissueName) 

	text(barP, nrGenes, nrGenes, pos=3)

	axis(1, at=barP, labels = FALSE)
	labels = c("Paralog genes", "Non-paralog genes")
	text(barP, par("usr")[3] - 0.1*max(nrGenes), srt = 45, adj = 1, labels = labels, xpd = TRUE)

dev.off()

#-----------------------------------------------------------------------
# Plot distribution of paralog group size
#-----------------------------------------------------------------------
# count paralog partners for each gene
countGenes = table(paralogPairs[,1])
countSizes = table(as.data.frame(countGenes)[,2])
countSizes = c("0"=length(nonParalogs), countSizes)
groupSizes = as.numeric(names(countSizes))

pdf(paste0(outPrefix, ".paralog_group_size.pdf"))
	par(cex=1.5, lwd=3)
	plot(groupSizes, as.vector(countSizes), xlim=c(0,40), type="o", pch=19, xlab="Number of paralogs", ylab="Number of genes", main="Paralog group size distribution", col=COL[1])
dev.off()

#-----------------------------------------------------------------------
# fraction of paralog pairs on same chromosom
#-----------------------------------------------------------------------
pdf(paste0(outPrefix, ".paralogPairs_on_same_chrom.barplot.pdf"), width=3.5)

	paraSameChrom <- !is.na(paralogPairsUniq$dist)
	sampSameChrom <- lapply(randPairs, function(gP) !is.na(gP$dist))

	paraPercent <- percentTrue(paraSameChrom)
	randPercent <- sapply(sampSameChrom, percentTrue) 
	
	height=c(paraPercent, mean(randPercent))
	
	# create contingency table and run Fisher test
	contab = rbind(
		para=table(paraSameChrom),  
		rand=table(unlist(sampSameChrom))
	)
	pval <- fisher.test(contab)$p.value
	
	yMax = 1.3*max(height) + sd(randPercent)
	par(cex=1.5, lwd=2)
	
	bp = my.barplot(height, 
		names=paste0(c("Paralogs\n (n=", "Random pairs\n (n="), c(nrow(paralogPairsUniq), paste0(N_RAND, "x",nrow(paralogPairsUniq))), c(")", ")")), 
		addValues=TRUE, yMax=yMax, col=COL_RAND,
		main=paste(tissueName, "gene pairs\n on same\n chromosome"), 
		ylab="Pairs on same chromosome [%]")
	error.bar(bp,height, c(NA, sd(randPercent)))

	# p-value matrix for only two group comparison
	pval_mat = matrix(c(NA, signif(pval, 3)), 1)
	add_pval_all_allPairs(bp, ymax=1.1*max(height)+ sd(randPercent), pval_mat, min_pval=10^-16)


dev.off()

#-----------------------------------------------------------------------
# Validate sampling by analysing distribution of distance, linked enhancer...
#-----------------------------------------------------------------------
pdf(paste0(outPrefix, ".sampling.Dist.pdf"), w=9, h=2.25)

	par(cex=1, lwd=1.5, mfrow=c(1,4))

	paraDist <- abs(allCisPairs[,"dist"] / 10^3)
	paraDist <- paraDist[paraDist <=1000]
	sampDist <- abs(sampCisPairsCombined[,"dist"] / 10^3)
	sampDist <- sampDist[sampDist <=1000]

	hist(paraDist, 20, col=COL[1],
	main="Paralogs", xlab="Distance [kb]")
	hist(sampDist, 20, col=COL[2],
	main="Sampled by distance", xlab="Distance [kb]")   
	qqplot(paraDist, sampDist, xlab="Paralogs distance [kb]", ylab="Sampled distance [kb]", main="QQ-Plot")
	abline(0,1, col="red")    
	qqplot(paraDist, sampDist, log="xy", xlab="Paralogs distance [kb]", ylab="Sampled distance [kb]", main="QQ-Plot\n(log scale)")
	abline(0,1, col="red")    

dev.off()

#-----------------------------------------------------------------------
# Show bias of equaly sampled genes to paralog genes
#-----------------------------------------------------------------------

# get paralog distances <= MAX_DIST
paraDist = abs(closePairs[,"dist"])/10^3
# get cispairs and distance form uniformaly random genes
randDist <- abs(randPairsInCisCombined[abs(randPairsInCisCombined[,"dist"])<= MAX_DIST, "dist"])/10^3
# distance form sampled pairs
sampledDist = abs(sampCisPairsCombined[abs(sampCisPairsCombined[,"dist"]) <= MAX_DIST,"dist"])/10^3

randDistPval <- wilcox.test(paraDist, randDist)$p.value
sampleDistPval <- wilcox.test(paraDist, sampledDist)$p.value
			
pdf(paste0(outPrefix, ".random_genes_distance.hist.pdf"))
	par(cex=1.5, lwd=2, mfrow=c(2,1))
	hist(paraDist, 50, col=COL_RAND[1],
	main="Distance between paralog genes", xlab="Distance (kb)")
	hist(randDist, 50, col=COL_RAND[2],
	main="Distance between random genes", xlab="Distance (kb)")    
dev.off()

pdf(paste0(outPrefix, ".samped_random_para_dist.hist.pdf"))

	par(lwd=2, mfrow=c(3,1))
	par(cex=1.5,  mar=c(3, 4.1, 1.5, 2.1))
	
	hp <- hist(paraDist, 50, col=COL_RAND[1],
	main=paste("Paralog gene pairs in", tissueName), xlab="")
	hr <- hist(randDist, 50, col=COL_RAND[2],
	main=paste0("Random gene pairs (p=", signif(randDistPval, 2), ")"), xlab="")
	hs <- hist(sampledDist, 50, col=COL[2],
	main=paste0("Sampled gene pairs (p=", signif(sampleDistPval, 2), ")"), xlab="")
	mtext("Distance (kb)", side=1, line=2, cex=1.5)    

dev.off()

#~ #-----------------------------------------------------------------------
#~ # Analyse TADs
#~ message("Start to analyse TADs...")

#~ # Plot fraction of paralog pairs within TAD
#~ # get boolean vectors
#~ paraInTAD = closePairs[, "TAD_Rao_CH12-LX"]
#~ sampInTAD = lapply(sampCisPairs, function(gP) gP[abs(gP$dist) <= MAX_DIST, "TAD_Rao_CH12-LX"])

#~ # create contingency table and run Fisher test
#~ contab = rbind(
#~ 	para=table(paraInTAD),  
#~ 	rand=table(unlist(sampInTAD))
#~ )
#~ fs.test = fisher.test(contab)

#~ # get percent values for barplot
#~ paraPercent = percentTrue(paraInTAD)
#~ sampPercent = sapply(sampInTAD, percentTrue) 
#~ heights = c(paraPercent, mean(sampPercent, na.rm=TRUE))
#~ sdRand = sd(sampPercent, na.rm=TRUE)

#~ pdf(paste0(outPrefix, ".paralogPairs.within_TAD.barplot.pdf"), width=3.5)

#~ 	par(cex=1.5, lwd=2)
#~ 	yMax = 1.3*max(heights) + sdRand
	
#~ 	bp = my.barplot(heights, addValues=TRUE, yMax=yMax,
#~ 		names=c("Paralogs", "Sampled"),
#~ 		main=paste("Pairs in", tissueName, "TAD"), ylab="Gene pairs in same TAD [%]", col=COL)

#~ 	error.bar(bp,heights, c(NA,  sdRand), lwd=2)
	
#~ 	# pvalue matrix for only two group comparision
#~ 	pval_mat = matrix(c(NA, signif(fs.test$p.value, 3)), 1)
#~ 	add_pval_all_allPairs(bp, ymax=1.1*max(heights)+sdRand, pval_mat)

#~ dev.off()
