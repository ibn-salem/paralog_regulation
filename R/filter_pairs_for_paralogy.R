########################################################################
#
# A script to filter a list of gene pairs (from HIPPIE, or coile-coile co-evolution) for
# paralog pairs.
#
# Approach:
# - Map UniProt/TrEMBL Accession or UniProt/SwissProt Accession to ENSG for each input pair
# - Filter against list of paralog ENSG pairs
#
########################################################################

require(biomaRt)        # to retrieve human paralogs from Ensembl
require(stringr)        # for some string functionality
require(RColorBrewer)   # for nice colors
require(ggplot2)        # for nice plots
require(BiocParallel)   # for parallel computing

#-----------------------------------------------------------------------
# Load parameters from external script
#-----------------------------------------------------------------------

VERSION="v01"
USE_LOCAL_DATA=TRUE

INFILES <- c("common_pairs_TRUE.txt", "results_all_filtered.txt")


outPrefix = paste0("results/filter_protein_pairs/", VERSION,".EnsemblGRCh37_paralog_genes")

# create directory, if not exist
dir.create(dirname(outPrefix), recursive=TRUE, showWarnings = FALSE)

# image file to save temporary and final R session with all data
WORKIMAGE_FILE = paste0(outPrefix, "workspace.Rdata")

#-----------------------------------------------------------------------
# Options for parallel computation
#-----------------------------------------------------------------------
# use all available cores but generate random number streams on each worker
RANDOM_SEED=13521
multicorParam <- MulticoreParam(RNGseed=RANDOM_SEED)   
# set options
register(multicorParam)  
# bpparam() # to print current options

#-----------------------------------------------------------------------
# load some custom functions
#-----------------------------------------------------------------------
source("R/functions.genePairs.R")

#=======================================================================
# 1.) Download paralog pairs with UniPort IDs
#=======================================================================
if ( !USE_LOCAL_DATA) {

    # define database and choose the human gene dataset
    # use last ensembl VERSION_DATA_ENSEMBL corresponding to human assembly GRCh37 (ensembl75) 
    ensemblGRCh37 <- useMart(host="grch37.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", verbose=FALSE)
    
    #-------------------------------------------------------------------
    # get UniProt/TrEMBL Accession to ENSG mapping
    #-------------------------------------------------------------------
    uniPTrEMBLAttributes <- c("uniprot_sptrembl", "ensembl_gene_id")
    uniPTrEMBLToEnsgDF = getBM(attributes=uniPTrEMBLAttributes, mart=ensemblGRCh37, 
        filters="with_uniprotsptrembl", values=TRUE)
    
    # save or load downloaded data 
    save(uniPTrEMBLToEnsgDF, file=paste0(outPrefix, ".uniPTrEMBLToEnsgDF.RData"))
    
    # take only unique uniPTrEMBL IDs
    uniPTrEMBLToEnsgDFuniq <- uniPTrEMBLToEnsgDF[!duplicated(uniPTrEMBLToEnsgDF$uniprot_sptrembl),]
    
    # convert data.frame into a named vector
    uniPTrEMBLToENSG <- uniPTrEMBLToEnsgDFuniq[,2]
    names(uniPTrEMBLToENSG) <- uniPTrEMBLToEnsgDFuniq[,1]
    
    # save or load downloaded data 
    save(uniPTrEMBLToENSG, file=paste0(outPrefix, ".uniPTrEMBLToENSG.RData"))
    
    #-------------------------------------------------------------------
    # get UniProt/SwissProt Accession to ENSG mapping
    #-------------------------------------------------------------------
    uniPSwissAttributes <- c("uniprot_swissprot_accession", "ensembl_gene_id")
    uniPSwissToEnsgDF = getBM(attributes=uniPSwissAttributes, mart=ensemblGRCh37)
    
    # filter for those with UniProt/Swiss prot Accession
    uniPSwissToEnsgDF <- uniPSwissToEnsgDF[uniPSwissToEnsgDF[,1] != "",]
    
    # save or load downloaded data 
    save(uniPSwissToEnsgDF, file=paste0(outPrefix, ".uniPSwissToEnsgDF.RData"))
    
    
    # take only unique uniPSwissT IDs 
    uniPSwissToEnsgDFuniq <- uniPSwissToEnsgDF[!duplicated(uniPSwissToEnsgDF$uniprot_swissprot_accession),]
    
    # convert data.frame into a named vector
    uniPSwissToEns <- uniPSwissToEnsgDFuniq[,2]
    names(uniPSwissToEns) <- uniPSwissToEnsgDFuniq[,1]
    
    # save or load downloaded data 
    save(uniPSwissToEns, file=paste0(outPrefix, ".uniPSwissToEns.RData"))
    
    #-------------------------------------------------------------------
    # Human paralog gene pairs
    #-------------------------------------------------------------------
    # get all human paralog gene pairs 
    paralogParisAttr = c("ensembl_gene_id", "hsapiens_paralog_ensembl_gene")
    paralogPairsALL = getBM(attributes=paralogParisAttr, filters="with_paralog_hsap", values=TRUE, mart=ensemblGRCh37)    
    save(paralogPairsALL, file=paste0(outPrefix, ".paralogPairsALL.RData"))
}else{
    message("Start to load ENSEMBL data from device...")
    load(paste0(outPrefix, ".uniPTrEMBLToEnsgDF.RData"))
    load(paste0(outPrefix, ".uniPTrEMBLToENSG.RData"))
    load(paste0(outPrefix, ".uniPSwissToEnsgDF.RData"))
    load(paste0(outPrefix, ".uniPSwissToEns.RData"))
    load(paste0(outPrefix, ".paralogPairsALL.RData"))
    message("Finished.")
}

# FILTER PARALOGS
possibleGenes <- unique(c(uniPTrEMBLToEnsgDF[,2], uniPSwissToEnsgDF[,2]))
paralogs <- paralogPairsALL[paralogPairsALL[,1] %in% possibleGenes & paralogPairsALL[,2] %in% possibleGenes, ]

# convert paralog paris to ID string consitiong of both IDs
paralogIdStr <- paste0(paralogs[,1], "_", paralogs[,2])

#=======================================================================
# Read input protein pairs
#=======================================================================

# iterate over all input files
for (inFile in INFILES) {
    
#~     inFile = INFILES[1]
    d = read.table(paste0("data/Pablo_CoilCoil/", inFile), header=FALSE, sep="|")
        
    #-----------------------------------------------------------------------
    # Test if input gene pairs are in set of paralog pairs
    #-----------------------------------------------------------------------
    isParalogPair <- function(i){
        
        if (i %% 1000 == 0){
            print(paste0("INFO: pair id: ", i))
        }
        
        p1 <- as.character(d[i,1])
        p2 <- as.character(d[i,2])
        
        # map to UniProt IDs
        g1 <- c(uniPTrEMBLToEnsgDF[uniPTrEMBLToEnsgDF[,1] == p1,2], uniPSwissToEnsgDF[uniPSwissToEnsgDF[,1] == p1,2])
        g2 <- c(uniPTrEMBLToEnsgDF[uniPTrEMBLToEnsgDF[,1] == p2,2], uniPSwissToEnsgDF[uniPSwissToEnsgDF[,1] == p2,2])
        
        # get all possible pairs between the two gene sets
        genePairs <- expand.grid(g1[!is.na(g1)], g2[!is.na(g2)])
        
        genePairIdStr <- paste0(genePairs[,1], "_", genePairs[,2])

        isParalog <- any(genePairIdStr %in% paralogIdStr)
                
        return(isParalog)
    }
    
    
    chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
    parts <- chunk2(1:nrow(d), 30)
    
    isParaList <- bplapply(parts, function(ids) sapply(ids, isParalogPair))
    
    d[,"paralog"] <- unlist(isParaList)
    
    # write to output file
    write.table(d, file=paste0(outPrefix, ".", inFile, ".annotParalog"), col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
    # write only filtered pairs
    write.table(d[d$paralog, 1:3], file=paste0(outPrefix, inFile, ".paralogs"), col.names=FALSE, row.names=FALSE, quote=FALSE, sep="|")

}


#=======================================================================
# save workspace image
#=======================================================================
save.image(WORKIMAGE_FILE)
