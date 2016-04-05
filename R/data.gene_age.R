#!/usr/bin/Rscript
#=======================================================================
#
#    gene pair age data from Pablo Mier Munoz 
#
#=======================================================================
require(stringr)

#=======================================================================
# parameters and input data files
#=======================================================================
VERSION_DATA_GENE_AGE="v01"

# gene pair age file
PAIR_AGE_FILE="data/Pablo/paralog_pairs_age_filtered.txt"

# singel gene file
GENE_AGE__FILE="data/Pablo/results.txt"

#=======================================================================
# parse data and save it as R object
#=======================================================================

# parse organism names and divergence age
age.org.str <- read.table(GENE_AGE__FILE, nrow=1, stringsAsFactors=FALSE)
age.org <- data.frame(str_split_fixed(age.org.str, "_", 2))
colnames(age.org) <- c("organism", "myr")

# read age data from file
pair.age.df <- read.table(PAIR_AGE_FILE, stringsAsFactors=FALSE)

# split up first column to get two vectors with swiss prot accession IDs
pair.ac <- data.frame(str_split_fixed(pair.age.df[,1], "\\+", 2))

# add age column to the new data.frame
pair.ac$age <- pair.age.df[,2] 

#map.df <- uniPSwissToEnsgDF

#=======================================================================
# adds the age information (from pair.ac) using a mapping of IDs (from uniPSwissToEnsgDF)
#=======================================================================
addAge <- function(genePair, pair.ac, map.df){
    
    ac1 <- map.df[match(genePair[,1], map.df[,2]), 1]
    ac2 <- map.df[match(genePair[,2], map.df[,2]), 1]

    # sort ac IDs from the input gene pairs to be unique
    acComb = sapply(apply(cbind(ac1, ac2), 1, sort), paste, collapse="+")    
    
    # also sort IDs of the data paralog pairs
    pair.ac.comb <- apply(apply(pair.ac[,1:2], 1, sort), 2, paste, collapse="+")
    
    # get age fro matching IDs
    age <- pair.ac[match(acComb, pair.ac.comb), 3]
    
    # add age column to the input gene pairs
    genePair[,"age"] <- age
    
    return(genePair)
}

