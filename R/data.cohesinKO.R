#!/usr/bin/Rscript
#=======================================================================
#
#   parse cohesin KO expression data
#
#=======================================================================

require(reshape2)	# for melt() and cast()


#=======================================================================
# parameters and input data files
#=======================================================================

COHESIN_KO_FILE="data/Dimitirs/RNAseqWTRad21KOMacrophages.RData"

#=======================================================================
# Parse data and reformat
#=======================================================================

# load R file:
cohesinEnv <- new.env()
load(COHESIN_KO_FILE, envi=cohesinEnv )

expDF <- cohesinEnv$FPKM.macrophages

# drop non-unique ENSG rows
expDF <- expDF[!duplicated( expDF[,"ensembl_gene_id"]),]

# set ENSG as row names and remove gene name and ID columns
row.names(expDF) <- expDF[,"ensembl_gene_id"]
rawExpDF <- expDF[,-(1:2)]

cd <- cohesinEnv$colData

expKO <- subset(rawExpDF, select=row.names(cd)[cd$Genotype == "Rad21KO"])
expWT <- subset(rawExpDF, select=row.names(cd)[cd$Genotype == "WT"])

#=======================================================================
# combine replicates
#=======================================================================

# melt data to only one value and one variable column
lexpDF <- melt(expDF, id=2, measure.vars=3:ncol(expDF))

# add metadata from colData
llexpDF <- cbind(lexpDF, cd[as.character(lexpDF$variable),])

# compute mean of replicates

repKO <- dcast(subset(llexpDF, Genotype=="Rad21KO"), ensembl_gene_id ~ Treatment + Time_hrs, mean, value.var="value", na.rm=TRUE)
repWT <- dcast(subset(llexpDF, Genotype=="WT"), ensembl_gene_id ~ Treatment + Time_hrs, mean, value.var="value", na.rm=TRUE)

# use ENSG as row.names and remove it from columns
row.names(repKO) <- repKO$ensembl_gene_id
row.names(repWT) <- repWT$ensembl_gene_id

repKO$ensembl_gene_id <- NULL
repWT$ensembl_gene_id <- NULL
	
expCoDFlist <- list("expWT"=expWT,
        "expKO"=expKO,
        "repWT"=repWT,
        "repKO"=repKO
)

# test that mean works 
#~ expKO["ENSMUSG00000000028", c("BM1_FL8", "BM2_FL8", "BM3_FL8")]
#~ sum(expKO["ENSMUSG00000000028", c("BM1_FL8", "BM2_FL8", "BM3_FL8")]) / 3
#~ names(repKO)
#~ cohesinEnv$colData
#~ repKO["ENSMUSG00000000028", ]
