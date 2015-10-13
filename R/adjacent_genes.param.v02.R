########################################################################
#
# This script contains parameters for the adjacent_genes.R analysis.
#
########################################################################

#-----------------------------------------------------------------------
# External data file paths
#-----------------------------------------------------------------------
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

# Mouse and dog TADs from Rudan et al. 2015
RudanFile = "data/Rudan2015/mmc2.xlsx"


#-----------------------------------------------------------------------
# Colors used for plotting
#-----------------------------------------------------------------------
COL_TAD = brewer.pal(8, "Set1")[c(3,5)]
COL_ORTHO = c(brewer.pal(12, "Paired")[5], brewer.pal(8, "Pastel2")[8])
COL_SPECIES = brewer.pal(8, "Accent")[1:2]

#-----------------------------------------------------------------------
# Parameters for data loading
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
# Parameters for randomizations
#-----------------------------------------------------------------------
RANDOM_SEED=13521

#-----------------------------------------------------------------------
# Version of the analysis and output file paths
#-----------------------------------------------------------------------

VERSION="v02"

outDataPrefix = "results/paralog_regulation/EnsemblGRCh37_paralog_genes"
outPrefix = paste0("results/adjacent_genes/", VERSION, "/", VERSION, ".EnsemblGRCh37_adjacent_genes")

# create directory, if not exist
dir.create(dirname(outPrefix), recursive=TRUE, showWarnings = FALSE)

# image file to save temporary and final R session with all data
WORKIMAGE_FILE = paste0(outPrefix, "workspace.Rdata")

