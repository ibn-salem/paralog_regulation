#!/usr/bin/Rscript
#=======================================================================
#
#   Check for sequence conservation in human TADs
#
#=======================================================================

# load some useful libraries
require(stringr)        # for functions like paste0()
require(GenomicRanges)  # for genomic ranges
#~ require(GenomicAlignments)  # for summarizeOverlaps function
require(biomaRt)        # to retrieve human paralogs from Ensembl
require(plyr)           # count() function

# load custom scripts
source("R/functions.GRanges.R")

#=======================================================================
# parameters and input data files
#=======================================================================
VERSION_DATA_ENSEMBL="v01"

SPECIES=c("mmusculus", "cfamiliaris")

chromFile = "data/hg19/hg19.genome"
realChromFile = "data/hg19/hg19.genome.realChroms"

chromFileMouse = "data/mm10/mm10.chrom.sizes"
chromFileDog = "data/canFam3/canFam3.chrom.sizes"

outPrefixDataEnsembl=paste0("results/data/ensembl.data.", VERSION_DATA_ENSEMBL)

# make directory if not exist already
dir.create(dirname(outPrefixDataEnsembl), showWarnings = FALSE)

# TODO use makefiel like dependency instead of the USE_LOCAL_DATA_ENSEMBL option
USE_LOCAL_DATA_ENSEMBL = TRUE
#~ USE_LOCAL_DATA_ENSEMBL = FALSE

#=======================================================================
# 1.) Downloads and format data 
#=======================================================================
if ( !USE_LOCAL_DATA_ENSEMBL) {

    # define database and choose the human gene dataset
    # use last ensembl VERSION_DATA_ENSEMBL corresponding to human assembly GRCh37 (ensembl75) 
    ensemblGRCh37 = useMart(host="grch37.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
    
    #-------------------------------------------------------------------
    # get all genes with annotation:
    #-------------------------------------------------------------------
    geneAttributes = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "transcript_start", "start_position", "end_position", "strand", "status", "gene_biotype")
    geneFilters="chromosome_name"
    # read "normal" human chromosome names (without fixes and patches)
    geneValues=c(1:22, "X", "Y")
    allGenes = getBM(attributes=geneAttributes, mart=ensemblGRCh37, filters=geneFilters, values=geneValues)
    #allGenes = getBM(attributes=geneAttributes, mart=ensemblGRCh37)
    
    # save or load downloaded data 
    save(allGenes, file=paste0(outPrefixDataEnsembl, ".ENSEMBL_GRCh37.allGenes.RData"))

    #-------------------------------------------------------------------
    # Human paralog gene pairs
    #-------------------------------------------------------------------
    # get all human paralog gene pairs 
    paralogParisAttr = c("ensembl_gene_id", "hsapiens_paralog_ensembl_gene", "hsapiens_paralog_perc_id", "hsapiens_paralog_perc_id_r1", "hsapiens_paralog_dn", "hsapiens_paralog_ds")
    paralogPairsALL = getBM(attributes=paralogParisAttr, filters="with_paralog_hsap", values=TRUE, mart=ensemblGRCh37)    
    save(paralogPairsALL, file=paste0(outPrefixDataEnsembl, ".ENSEMBL_GRCh37.paralogPairsALL.RData"))

    #-------------------------------------------------------------------
    # Get for all human known coding genes the orthologs in mouse and dogs
    #-------------------------------------------------------------------
    orthologsSpeciesList = lapply(SPECIES, function(orgStr) {
    
        # atributes 
        orthologAttr = c("ensembl_gene_id",
            paste0(orgStr, c("_homolog_ensembl_gene", "_homolog_orthology_type", "_homolog_subtype", "_homolog_orthology_confidence", "_homolog_perc_id", "_homolog_perc_id_r1", "_homolog_dn", "_homolog_ds")))
            
        orthologs = getBM(attributes=orthologAttr, filters=c("status", "biotype","chromosome_name"), values=list(status="KNOWN", biotype="protein_coding", chromosome_name=c(1:22, "X", "Y")), mart=ensemblGRCh37)  
        
        # for genes that do not have orthologs set empty string to NA
        orthologs[orthologs[,2] == "", 2:4] <- NA
        
        return(orthologs)
    })
    names(orthologsSpeciesList) = SPECIES
    
    save(orthologsSpeciesList, file=paste0(outPrefixDataEnsembl, ".ENSEMBL_GRCh37.orthologsSpeciesList.RData"))
    
    #-------------------------------------------------------------------
    # get all exon regions
    #-------------------------------------------------------------------
    # get all known, coding exon IDs on chromosome 1-22, X, and Y
#~     exonIDs = getBM(attributes="ensembl_exon_id", mart=ensemblGRCh37, 
#~         filters=c("status", "biotype","chromosome_name"), 
#~         values=list(status="KNOWN", biotype="protein_coding", chromosome_name=c(1:22, "X", "Y")),
#~         )[,1]
#~ 
#~     allExonIDs = getBM(mart=ensemblGRCh37, 
#~         attributes=c("ensembl_exon_id", "gene_biotype"),
#~         filters=c("status","chromosome_name"), 
#~         values=list(status="KNOWN", chromosome_name=c(1:22, "X", "Y")),
#~         )

#~     exonAttributesStructure = c("ensembl_exon_id", "chromosome_name", "exon_chrom_start", "exon_chrom_end", "rank")
#~     codingExons = getBM(attributes=exonAttributesStructure, mart=ensemblGRCh37, filters="ensembl_exon_id", values=exonIDs)
    
#~     exonAttributesStructure = c("ensembl_exon_id", "chromosome_name", "exon_chrom_start", "exon_chrom_end", "rank", "gene_biotype")
#~     allExons = getBM(attributes=exonAttributesStructure, mart=ensemblGRCh37, filters="ensembl_exon_id", values=allExonIDs[,1])
    
    # get all exons with status KNOWN on chromosome 1-22, X, and Y
    allExons = getBM(mart=ensemblGRCh37, 
        attributes=c("ensembl_exon_id", "chromosome_name", "exon_chrom_start", "exon_chrom_end", "strand", "rank", "gene_biotype"), 
        filters=c("status","chromosome_name"), 
        values=list(status="KNOWN", chromosome_name=c(1:22, "X", "Y")))

    save(allExons, file=paste0(outPrefixDataEnsembl, ".ENSEMBL_GRCh37.allExons.RData"))

    #-------------------------------------------------------------------
    # Download mouse genes and paralog pairs
    #-------------------------------------------------------------------
    # define database and choose the mouse gene dataset
    # use last ensembl VERSION_DATA_ENSEMBL corresponding to mouse assembly GRCm38 (ensembl75) 
    ensemblMouse = useMart(host="grch37.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl")

    # get all genes with annotation:
    #-------------------------------------------------------------------
    MouseGeneAttributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name", "transcript_start", "start_position", "end_position", "strand", "status", "gene_biotype")

    # read "normal" human chromosome names (without fixes and patches)
    allMouseGenes = getBM(attributes=MouseGeneAttributes, mart=ensemblMouse, filters="chromosome_name", values=c(1:19, "X", "Y"))
    #allGenes = getBM(attributes=geneAttributes, mart=ensemblGRCh37)

    # save or load downloaded data 
    save(allMouseGenes, file=paste0(outPrefixDataEnsembl, ".ENSEMBL_GRCh37.allMouseGenes.RData"))

    # Mouse paralog gene pairs
    paralogParisMouseAttr = c("ensembl_gene_id", "mmusculus_paralog_ensembl_gene", "mmusculus_paralog_perc_id", "mmusculus_paralog_perc_id_r1", "mmusculus_paralog_dn", "mmusculus_paralog_ds")     
    paralogPairsMouseALL = getBM(attributes=paralogParisMouseAttr, filters="with_paralog_mmus", values=TRUE, mart=ensemblMouse)    
    
    save(paralogPairsMouseALL, file=paste0(outPrefixDataEnsembl, ".ENSEMBL_GRCh37.paralogPairsMouseALL.RData"))

    # mouse to human orthologs to get one (in mouse) to many (in humna)
    mouseOrthologsInHumanAttr = c("ensembl_gene_id", "hsapiens_homolog_ensembl_gene", "hsapiens_homolog_orthology_type", "hsapiens_homolog_subtype", "hsapiens_homolog_perc_id", "hsapiens_homolog_perc_id_r1", "hsapiens_homolog_dn", "hsapiens_homolog_ds")     
    mouseOrthologsInHumanALL = getBM(attributes=mouseOrthologsInHumanAttr, filters="with_homolog_hsap", values=TRUE, mart=ensemblMouse)    

    save(mouseOrthologsInHumanALL, file=paste0(outPrefixDataEnsembl, ".ENSEMBL_GRCh37.mouseOrthologsInHumanALL.RData"))

    #-------------------------------------------------------------------
    # Download DOG genes and paralog pairs
    #-------------------------------------------------------------------
    # define database and choose the dog gene dataset
    # use last ensembl VERSION_DATA_ENSEMBL corresponding to dog assembly CanFam3.1 (ensembl75) 
    ensemblDog = useMart(host="grch37.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL",dataset="cfamiliaris_gene_ensembl")

    # get all genes with annotation:
    #-------------------------------------------------------------------
    DogGeneAttributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name", "transcript_start", "start_position", "end_position", "strand", "status", "gene_biotype")

    # read "normal" human chromosome names (without fixes and patches)
    allDogGenes = getBM(attributes=DogGeneAttributes, mart=ensemblDog, filters="chromosome_name", values=c(1:38, "X", "Y"))

    # save or load downloaded data 
    save(allDogGenes, file=paste0(outPrefixDataEnsembl, ".ENSEMBL_GRCh37.allDogGenes.RData"))

    # Dog paralog gene pairs
    paralogParisDogAttr = c("ensembl_gene_id", "cfamiliaris_paralog_ensembl_gene", "cfamiliaris_paralog_perc_id", "cfamiliaris_paralog_perc_id_r1", "cfamiliaris_paralog_dn", "cfamiliaris_paralog_ds")     
    paralogPairsDogALL = getBM(attributes=paralogParisDogAttr, filters="with_paralog_cfam", values=TRUE, mart=ensemblDog)    
    
    save(paralogPairsDogALL, file=paste0(outPrefixDataEnsembl, ".ENSEMBL_GRCh37.paralogPairsDogALL.RData"))

    

# if USE_LOCAL_DATA_ENSEMBL, just load the saved data from local hard disc
}else{
    message("Start to load ENSEMBL data from device...")
    load(paste0(outPrefixDataEnsembl, ".ENSEMBL_GRCh37.allGenes.RData"))
    load(paste0(outPrefixDataEnsembl, ".ENSEMBL_GRCh37.paralogPairsALL.RData"))
    load(paste0(outPrefixDataEnsembl, ".ENSEMBL_GRCh37.orthologsSpeciesList.RData"))
    load(paste0(outPrefixDataEnsembl, ".ENSEMBL_GRCh37.allExons.RData"))
    load(paste0(outPrefixDataEnsembl, ".ENSEMBL_GRCh37.allMouseGenes.RData"))
    load(paste0(outPrefixDataEnsembl, ".ENSEMBL_GRCh37.paralogPairsMouseALL.RData"))
    load(paste0(outPrefixDataEnsembl, ".ENSEMBL_GRCh37.mouseOrthologsInHumanALL.RData"))
    load(paste0(outPrefixDataEnsembl, ".ENSEMBL_GRCh37.allDogGenes.RData"))
    load(paste0(outPrefixDataEnsembl, ".ENSEMBL_GRCh37.paralogPairsDogALL.RData"))
    message("Finished.")
}

#=======================================================================
# 2.) Convert data into GRanges 
#=======================================================================

# parse seqInfo object
seqInfo = parseSeqInfo(chromFile)
seqInfoRealChrom = parseSeqInfo(realChromFile)

#-----------------------------------------------------------------------
# subset of known coding Genes
#-----------------------------------------------------------------------
# filter for known genes only and only for protein coding genes
knownCodingGenes = allGenes[allGenes$status=="KNOWN" & allGenes$gene_biotype== "protein_coding",]

# unique gene entry by ENSG ID symbol:
genes = knownCodingGenes[!duplicated(knownCodingGenes$ensembl_gene_id),]
rownames(genes) <- genes$ensembl_gene_id

# unique gene entry by HGNC symbol:
genesHGNC = knownCodingGenes[!duplicated(knownCodingGenes$hgnc_symbol),]
rownames(genesHGNC) <- genesHGNC$hgnc_symbol

# map ENSG to HGNC symbols
gene2hgnc = as.list(knownCodingGenes$hgnc_symbol)
names(gene2hgnc) = knownCodingGenes$ensembl_gene_id

# make GRanges object for all known prot coding genes
genesGR = getGenesGR(genes, seqInfo)

# make GRanges object for all genes including non-coding
knownGenes = allGenes[allGenes$status=="KNOWN",]
allGenesGR = getGenesGR(knownGenes[!duplicated(knownGenes$ensembl_gene_id),], seqInfo)

#-----------------------------------------------------------------------
# Convert Mouse and Dog genes to GRanges
#-----------------------------------------------------------------------
# parse seqInfo object
seqInfoMouse = parseSeqInfo(chromFileMouse, genome="mm10", header=FALSE, filterForReal=TRUE)

knownCodingGenesMouse = allMouseGenes[allMouseGenes$status=="KNOWN" & allMouseGenes$gene_biotype== "protein_coding",]

# unique gene entry by ENSG ID symbol:
genesMouse = knownCodingGenesMouse[!duplicated(knownCodingGenesMouse$ensembl_gene_id),]
rownames(genesMouse) <- genesMouse$ensembl_gene_id

# make GRange object for the TSS positions
tssGRmouse = getTssGRfromENSEMBLGenes(genesMouse, seqInfoMouse, colNames=c("external_gene_name"))
tssGRmouse$gene_size = genesMouse[names(tssGRmouse), "end_position"] - genesMouse[names(tssGRmouse), "start_position"]

# same for DOG:

# parse seqInfo object
seqInfoDog = parseSeqInfo(chromFileDog, genome="canFam3", header=FALSE, filterForReal=TRUE)

knownCodingGenesDog = allDogGenes[allDogGenes$status %in% c("KNOWN", "KNOWN_BY_PROJECTION") & allDogGenes$gene_biotype== "protein_coding",]

# unique gene entry by ENSG ID symbol:
genesDog = knownCodingGenesDog[!duplicated(knownCodingGenesDog$ensembl_gene_id),]
rownames(genesDog) <- genesDog$ensembl_gene_id

# make GRange object for the TSS positions
tssGRdog = getTssGRfromENSEMBLGenes(genesDog, seqInfoDog, colNames=c("external_gene_name"))
tssGRdog$gene_size = genesDog[names(tssGRdog), "end_position"] - genesDog[names(tssGRdog), "start_position"]

#-------------------------------------------------------------------
# Filter human paralogs and get non-paralog gene sets
#-------------------------------------------------------------------
# filter out non-coding or novel paralogs:
paralogPairs = paralogPairsALL[paralogPairsALL[,1] %in% rownames(genes) & paralogPairsALL[,2] %in% rownames(genes),]

# get IDs of human paralog and non-paralog gene sets
paralogs = genes[unique(paralogPairs[,1]), "ensembl_gene_id"]
paralogsHGNC = genesHGNC[unlist(gene2hgnc[paralogs]), "hgnc_symbol"]

nonParalogs = genes[setdiff(genes$ensembl_gene_id, paralogs), "ensembl_gene_id"]
nonParalogsHGNC = genesHGNC[unlist(gene2hgnc[nonParalogs]), "hgnc_symbol"]

#-------------------------------------------------------------------
# Filter mouse and dog paralogs
#-------------------------------------------------------------------
paralogPairsMouse = paralogPairsMouseALL[paralogPairsMouseALL[,1] %in% rownames(genesMouse) & paralogPairsMouseALL[,2] %in% rownames(genesMouse),]
    
paralogPairsDog = paralogPairsDogALL[paralogPairsDogALL[,1] %in% rownames(genesDog) & paralogPairsDogALL[,2] %in% rownames(genesDog),]

#-----------------------------------------------------------------------
# convert exons to GRanges objects
#-----------------------------------------------------------------------
exonGR = getExonGR(allExons[!duplicated(allExons$ensembl_exon_id),], seqInfo)

#~ codingExonGR = getExonGR(allExons[allExons$gene_biotype == "protein_coding",], seqInfo)
codingExonGR = exonGR[exonGR$gene_biotype == "protein_coding"]

ncExonGR = exonGR[exonGR$gene_biotype != "protein_coding"]
# Note, that some non-coding exons might overlap with codig exons


#-----------------------------------------------------------------------
# set of introns and intergenci regions
#-----------------------------------------------------------------------

# get an Grange with all genic regions without strand information
unstrandedAllGenesGR = noStrand(allGenesGR)

# get the gaps between the gene ranges as unstranded intergenic ranges
#~ intergenicGR = gaps(unstrandedGenesGR)
#~ intergenicGR = intergenicGR[strand(intergenicGR) == "*"]

genomeGR = as(seqInfo,'GRanges')
intergenicGR = setdiff(genomeGR, reduce(unstrandedAllGenesGR), ignore.strand=TRUE)

# tests
stopifnot(covBases(allGenesGR)  + covBases(intergenicGR)  == covBases(genomeGR) )

# get intron regions as gaps between exons that are completely contained (within) a gene region
intronGR = subsetByOverlaps(gaps(exonGR), allGenesGR, type="within", ignore.strand=TRUE)

# this would be an more exclusive VERSION_DATA_ENSEMBL that 
#intronGR = subsetByOverlaps(filterForNoStrand(gaps(noStrand(exonGR))), unstrandedAllGenesGR, type="within")

#~ intronGR = setdiff(allGenesGR, reduce(exonGR, ignore.strand=TRUE), ignore.strand=TRUE)

#otherGenicGR = setdiff(reduce(unstrandedAllGenesGR, ignore.strand=TRUE), reduce(noStrand(c(exonGR[,NULL], intronGR))))

# tests
#stopifnot( covBases(intronGR)  + covBases(exonGR) + covBases(otherGenicGR)  == covBases(allGenesGR) )

#-----------------------------------------------------------------------
# set of non-coding RNAs
#-----------------------------------------------------------------------
ncRNAs = allGenes[allGenes$gene_biotype %in% c("lincRNA","miRNA","rRNA"),]
ncRNAsUNIQ = ncRNAs[!duplicated(ncRNAs$ensembl_gene_id),]
# convert it GRanges object
ncRNAsGR = GRanges(
        paste0("chr", ncRNAsUNIQ$chromosome_name),
        IRanges(ncRNAsUNIQ$start_position, ncRNAsUNIQ$end_position),
        strand=ifelse(ncRNAsUNIQ$strand == 1, "+", "-"),
        ncRNAsUNIQ[,c("hgnc_symbol", "status", "gene_biotype")],
        seqinfo=seqInfo
    )
names(ncRNAsGR) <- ncRNAsUNIQ$ensembl_gene_id

# do the same for only the lincRNAs
lincRNAs = allGenes[allGenes$gene_biotype == "lincRNA",]
lincRNAsUNIQ = lincRNAs[!duplicated(lincRNAs$ensembl_gene_id),]

# convert it GRanges object
lincRNAsGR = GRanges(
        paste0("chr", lincRNAsUNIQ$chromosome_name),
        IRanges(lincRNAsUNIQ$start_position, lincRNAsUNIQ$end_position),
        strand=ifelse(lincRNAsUNIQ$strand == 1, "+", "-"),
        lincRNAsUNIQ[,c("hgnc_symbol", "status", "gene_biotype")],
        seqinfo=seqInfo
    )
names(lincRNAsGR) <- lincRNAsUNIQ$ensembl_gene_id

#-----------------------------------------------------------------------
# set of all non-protein-coding genes
#-----------------------------------------------------------------------
npc = allGenes[allGenes$gene_biotype != "protein_coding",]
npcUNIQ = npc[!duplicated(npc$ensembl_gene_id),]
# convert it GRanges object
npcGR = GRanges(
        paste0("chr", npcUNIQ$chromosome_name),
        IRanges(npcUNIQ$start_position, npcUNIQ$end_position),
        strand=ifelse(npcUNIQ$strand == 1, "+", "-"),
        npcUNIQ[,c("hgnc_symbol", "status", "gene_biotype")],
        seqinfo=seqInfo
    )
names(npcGR) <- npcUNIQ$ensembl_gene_id
