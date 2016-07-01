#!/usr/bin/Rscript
#=======================================================================
#
#   Download and format genes and paralog gene pairs from ENSEMBL.org
#
#=======================================================================

# load some useful libraries
require(stringr)        # for functions like paste0()
require(GenomicRanges)  # for genomic ranges
require(biomaRt)        # to retrieve human paralogs from Ensembl
require(plyr)           # count() function

# load custom scripts
# require the following script to be loaded already:
#~ source("R/functions.GRanges.R")

#=======================================================================
# parameters and input data files
#=======================================================================
VERSION_DATA_ENSEMBL="v01"

SPECIES=c("mmusculus", "cfamiliaris")

chromFile = "data/hg19/hg19.genome"
realChromFile = "data/hg19/hg19.genome.realChroms"

chromFileMouse = "data/mm10/mm10.chrom.sizes"
chromFileDog = "data/canFam3/canFam3.chrom.sizes"

OLFACTORY_RECEPTOR_GO="GO:0004984"

outPrefixDataEnsembl=paste0("results/data/ensembl.data.", VERSION_DATA_ENSEMBL)

# make directory if not exist already
dir.create(dirname(outPrefixDataEnsembl), showWarnings = FALSE, recursive=TRUE)

USE_LOCAL_DATA_ENSEMBL = FALSE

#=======================================================================
# 1.) Downloads and format data 
#=======================================================================
if ( !USE_LOCAL_DATA_ENSEMBL) {

    # define database and choose the human gene dataset
    # use last ensembl VERSION_DATA_ENSEMBL corresponding to human assembly GRCh37 (ensembl75) 
    ensemblGRCh37 <- useMart(host="grch37.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", verbose=FALSE)
    
    #-------------------------------------------------------------------
    # get all genes with annotation:
    #-------------------------------------------------------------------
    geneAttributes = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position", "strand", "status", "gene_biotype")
    geneFilters="chromosome_name"
    # read "normal" human chromosome names (without fixes and patches)
    geneValues=c(1:22, "X", "Y")
    allGenes = getBM(attributes=geneAttributes, mart=ensemblGRCh37, filters=geneFilters, values=geneValues)
    
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
    # get RefSeq ID to ENSG mapping
    #-------------------------------------------------------------------
    refSeqAttributes <- c("refseq_mrna", "ensembl_gene_id")
    refSeqFilters <- c("chromosome_name", "with_refseq_mrna")
    refSeqValues <- list("chromosome_name"=c(1:22, "X", "Y"), with_refseq_mrna=TRUE) 
    refSeqToEnsgDF = getBM(attributes=refSeqAttributes, mart=ensemblGRCh37, filters=refSeqFilters, values=refSeqValues)
    
    # take only unique refSeq IDs (This will delete 20 out of 35945)
    refSeqToEnsgDF <- refSeqToEnsgDF[!duplicated(refSeqToEnsgDF$refseq_mrna),]
    
    # convert data.frame into a named vector
    refSeqToENSG <- refSeqToEnsgDF[,2]
    names(refSeqToENSG) <- refSeqToEnsgDF[,1]
    
    # save or load downloaded data 
    save(refSeqToENSG, file=paste0(outPrefixDataEnsembl, ".ENSEMBL_GRCh37.refSeqToENSG.RData"))

    #-------------------------------------------------------------------
    # get EntrezGene ID to ENSG mapping
    #-------------------------------------------------------------------
    entrezAttributes <- c("entrezgene", "ensembl_gene_id")
    entrezFilters <- c("chromosome_name", "with_entrezgene")
    entrezValues <- list("chromosome_name"=c(1:22, "X", "Y"), with_entrezgene=TRUE) 
    entrezToEnsgDF = getBM(attributes=entrezAttributes, mart=ensemblGRCh37, filters=entrezFilters, values=entrezValues)
    
    # take only unique refSeq IDs (This will delete 2088 out of 27791)
    entrezToEnsgDF <- entrezToEnsgDF[!duplicated(entrezToEnsgDF$entrezgene),]
    
    # convert data.frame into a named vector
    entrezToENSG <- entrezToEnsgDF[,2]
    names(entrezToENSG) <- entrezToEnsgDF[,1]
    
    # save or load downloaded data 
    save(entrezToENSG, file=paste0(outPrefixDataEnsembl, ".ENSEMBL_GRCh37.entrezToENSG.RData"))

    #-------------------------------------------------------------------
    # Download mouse genes and paralog pairs
    #-------------------------------------------------------------------
    # define database and choose the mouse gene dataset
    # use last ensembl VERSION_DATA_ENSEMBL corresponding to mouse assembly GRCm38 (ensembl75) 
    ensemblMouse = useMart(host="grch37.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl")

    # get all genes with annotation:
    #-------------------------------------------------------------------
    MouseGeneAttributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position", "end_position", "strand", "status", "gene_biotype")

    # read "normal" human chromosome names (without fixes and patches)
    allMouseGenes = getBM(attributes=MouseGeneAttributes, mart=ensemblMouse, filters="chromosome_name", values=c(1:19, "X", "Y"))
    #allGenes = getBM(attributes=geneAttributes, mart=ensemblGRCh37)

    # save or load downloaded data 
    save(allMouseGenes, file=paste0(outPrefixDataEnsembl, ".ENSEMBL_GRCh37.allMouseGenes.RData"))

    # Mouse paralog gene pairs
    paralogParisMouseAttr = c("ensembl_gene_id", "mmusculus_paralog_ensembl_gene", "mmusculus_paralog_perc_id", "mmusculus_paralog_perc_id_r1", "mmusculus_paralog_dn", "mmusculus_paralog_ds")     
    paralogPairsMouseALL = getBM(attributes=paralogParisMouseAttr, filters="with_paralog_mmus", values=TRUE, mart=ensemblMouse)    
    
    save(paralogPairsMouseALL, file=paste0(outPrefixDataEnsembl, ".ENSEMBL_GRCh37.paralogPairsMouseALL.RData"))

    #-------------------------------------------------------------------
    # Download DOG genes and paralog pairs
    #-------------------------------------------------------------------
    # define database and choose the dog gene dataset
    # use last ensembl VERSION_DATA_ENSEMBL corresponding to dog assembly CanFam3.1 (ensembl75) 
    ensemblDog = useMart(host="grch37.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL",dataset="cfamiliaris_gene_ensembl")

    # get all genes with annotation:
    #-------------------------------------------------------------------
    DogGeneAttributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name", "start_position", "end_position", "strand", "status", "gene_biotype")

    # read "normal" human chromosome names (without fixes and patches)
    allDogGenes = getBM(attributes=DogGeneAttributes, mart=ensemblDog, filters="chromosome_name", values=c(1:38, "X", "Y"))

    # save or load downloaded data 
    save(allDogGenes, file=paste0(outPrefixDataEnsembl, ".ENSEMBL_GRCh37.allDogGenes.RData"))

    # Dog paralog gene pairs
    paralogParisDogAttr = c("ensembl_gene_id", "cfamiliaris_paralog_ensembl_gene", "cfamiliaris_paralog_perc_id", "cfamiliaris_paralog_perc_id_r1", "cfamiliaris_paralog_dn", "cfamiliaris_paralog_ds")     
    paralogPairsDogALL = getBM(attributes=paralogParisDogAttr, filters="with_paralog_cfam", values=TRUE, mart=ensemblDog)    
    
    save(paralogPairsDogALL, file=paste0(outPrefixDataEnsembl, ".ENSEMBL_GRCh37.paralogPairsDogALL.RData"))

    #-------------------------------------------------------------------
    # get UniProt/SwissProt Accession to ENSG mapping
    #-------------------------------------------------------------------
    uniPSwissAttributes <- c("uniprot_swissprot_accession", "ensembl_gene_id")
    uniPSwissToEnsgDF = getBM(attributes=uniPSwissAttributes, mart=ensemblGRCh37)
    
    # filter for those with UniProt/Swiss prot Accession
    uniPSwissToEnsgDF <- uniPSwissToEnsgDF[uniPSwissToEnsgDF[,1] != "",]
    
    # save or load downloaded data 
    save(uniPSwissToEnsgDF, file=paste0(outPrefixDataEnsembl, ".ENSEMBL_GRCh37.uniPSwissToEnsgDF.RData"))
    
    # take only unique uniPSwissT IDs 
    uniPSwissToEnsgDFuniq <- uniPSwissToEnsgDF[!duplicated(uniPSwissToEnsgDF$uniprot_swissprot_accession),]
    
    # convert data.frame into a named vector
    uniPSwissToEns <- uniPSwissToEnsgDFuniq[,2]
    names(uniPSwissToEns) <- uniPSwissToEnsgDFuniq[,1]
    
    # save or load downloaded data 
    save(uniPSwissToEns, file=paste0(outPrefixDataEnsembl, ".ENSEMBL_GRCh37.uniPSwissToEns.RData"))
    
# if USE_LOCAL_DATA_ENSEMBL, just load the saved data from local hard disc
}else{
    message("Start to load ENSEMBL data from device...")
    load(paste0(outPrefixDataEnsembl, ".ENSEMBL_GRCh37.allGenes.RData"))
    load(paste0(outPrefixDataEnsembl, ".ENSEMBL_GRCh37.paralogPairsALL.RData"))
    load(paste0(outPrefixDataEnsembl, ".ENSEMBL_GRCh37.orthologsSpeciesList.RData"))
    load(paste0(outPrefixDataEnsembl, ".ENSEMBL_GRCh37.refSeqToENSG.RData"))
    load(paste0(outPrefixDataEnsembl, ".ENSEMBL_GRCh37.entrezToENSG.RData"))
    load(paste0(outPrefixDataEnsembl, ".ENSEMBL_GRCh37.allMouseGenes.RData"))
    load(paste0(outPrefixDataEnsembl, ".ENSEMBL_GRCh37.paralogPairsMouseALL.RData"))
    load(paste0(outPrefixDataEnsembl, ".ENSEMBL_GRCh37.allDogGenes.RData"))
    load(paste0(outPrefixDataEnsembl, ".ENSEMBL_GRCh37.paralogPairsDogALL.RData"))
    load(paste0(outPrefixDataEnsembl, ".ENSEMBL_GRCh37.uniPSwissToEnsgDF.RData"))
    load(paste0(outPrefixDataEnsembl, ".ENSEMBL_GRCh37.uniPSwissToEns.RData"))
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
genesGR = sort(getGenesGR(genes, seqInfo), ignore.strand=TRUE)

# make GRanges object for all genes including non-coding
knownGenes = allGenes[allGenes$status=="KNOWN",]
allGenesGR = sort(getGenesGR(knownGenes[!duplicated(knownGenes$ensembl_gene_id),], seqInfo), ignore.strand=TRUE)

# make GenomicRange objects for TSS of Genes:
tssGR = getTssGRfromENSEMBLGenes(genes, seqInfoRealChrom, colNames=c("hgnc_symbol"))
tssGR$gene_size = genes[names(tssGR), "end_position"] - genes[names(tssGR), "start_position"]
tssGR <- sort(tssGR, ignore.strand=TRUE)

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

# make GRanges object for all known prot coding genes
genesGRmouse = sort(getGenesGR(genesMouse, seqInfoMouse, annotCols=c("external_gene_name")), ignore.strand=TRUE)

#-----------------------------------------------------------------------
# same for DOG:
#-----------------------------------------------------------------------

# parse seqInfo object
seqInfoDog = parseSeqInfo(chromFileDog, genome="canFam3", header=FALSE, filterForReal=TRUE)

knownCodingGenesDog = allDogGenes[allDogGenes$status %in% c("KNOWN", "KNOWN_BY_PROJECTION") & allDogGenes$gene_biotype== "protein_coding",]

# unique gene entry by ENSG ID symbol:
genesDog = knownCodingGenesDog[!duplicated(knownCodingGenesDog$ensembl_gene_id),]
rownames(genesDog) <- genesDog$ensembl_gene_id

# make GRange object for the TSS positions
tssGRdog = getTssGRfromENSEMBLGenes(genesDog, seqInfoDog, colNames=c("external_gene_name"))
tssGRdog$gene_size = genesDog[names(tssGRdog), "end_position"] - genesDog[names(tssGRdog), "start_position"]

# make GRanges object for all known prot coding genes
genesGRdog = sort(getGenesGR(genesDog, seqInfoDog, annotCols=c("external_gene_name")), ignore.strand=TRUE)

# combine tssGR's form all species to list
speciesTssGR <- list(
    "mmusculus" = tssGRmouse,
    "cfamiliaris" = tssGRdog
)
speciesGenesGR <- list(
    "mmusculus" = genesGRmouse,
    "cfamiliaris" = genesGRdog
)

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

