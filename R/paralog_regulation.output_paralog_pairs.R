################################################################################
#
# A script to create a table with paralog pairs from the allDF data frame
#
################################################################################

require(GenomicRanges)  # for genomic intervals and overlap calculation
require(tidyverse)

#-----------------------------------------------------------------------
# Load parameters from external script
#-----------------------------------------------------------------------

PARAM_SCRIPT="R/paralog_regulation.param.R"
source(PARAM_SCRIPT)

# #-----------------------------------------------------------------------
# # load some custom functions
# #-----------------------------------------------------------------------
# source("R/functions.plot.R")
# source("R/functions.GRanges.R")
# source("R/functions.genePairs.R")
# source("R/functions.genePairs.randomization.R")
# source("R/functions.genePairs.paralog_analysis.R")
# source("R/parseHiC.R")



#=======================================================================
# 1.) Load data from exported data.frames
#=======================================================================

#~ INallDF <- read.delim(paste0(outPrefix, ".allDF.csv"), header=TRUE, sep="\t")

# load aP and allDF
message(paste("INFO: Start loading data from file:", paste0(WORKIMAGE_FILE, ".aP_allDF.Rdata")))
load(paste0(WORKIMAGE_FILE, ".aP_allDF.Rdata"))

# define some often used subsetting vectors for the allDF
sSamp <- allDF$sampType %in% c("paralogs", "sampled by distance")
sSampEh <- allDF$sampType %in% c("paralogs", "sampled by distance\nand enhancer")
sTAD <- allDF$tadSource == "Rao_IMR90"
sExp <- allDF$expSource == "GTEx"
sSpec <- allDF$species == "mouse"

#=======================================================================
# 2.) Define subset for output
#=======================================================================
subDF <- allDF %>% 
  as.tibble() %>% 
  filter(
    group == "paralog",
    sTAD,
    sExp,
    sSpec
  ) %>% 
  select(g1, g2, HGNC_g1, HGNC_g2, everything())

# write to .tsv file
write_tsv(subDF, paste0(outPrefix, ".paralog_cis_pairs.with_annotation.tsv"))

