########################################################################
#
# This script contains parameters for the paralog_regulation.R analysis.
#
########################################################################

require(RColorBrewer)   # for nice colors
require(colorspace)     # for some more colors

#-----------------------------------------------------------------------
# External data file paths
#-----------------------------------------------------------------------
REGMAP_FILE_FANTOM5 <- "data/Andersson2014/enhancer_tss_associations.bed"
EH_FILE_FANTOM5 = "data/Andersson2014/permissive_enhancers.bed"


RaoDomainFiles <- c(
    Rao_HeLa="data/Rao2014/GSE63525_HeLa_Arrowhead_domainlist.txt",
    Rao_HUVEC="data/Rao2014/GSE63525_HUVEC_Arrowhead_domainlist.txt",
    Rao_K562="data/Rao2014/GSE63525_K562_Arrowhead_domainlist.txt",
    Rao_KBM7="data/Rao2014/GSE63525_KBM7_Arrowhead_domainlist.txt",
    Rao_NHEK="data/Rao2014/GSE63525_NHEK_Arrowhead_domainlist.txt",
    Rao_IMR90="data/Rao2014/GSE63525_IMR90_Arrowhead_domainlist.txt",
    Rao_GM12878="data/Rao2014/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt"
)

RaoSubcompartmentFile <- "data/Rao2014/GSE63525_GM12878_subcompartments.bed"

DixonDomainFiles = c(
    Dixon_hESC="data/Dixon2012/hESC.hg18.bed.hg19.bed",
    Dixon_IMR90="data/Dixon2012/IMR90.hg18.bed.hg19.bed"
)

# Mouse and dog TADs from Rudan et al. 2015
RudanFile <- "data/Rudan2015/mmc2.xlsx"


subTADfigPaths <- paste0("data/figs/gene_pairs_", c(
    "no_TAD",
    "diff_TAD",
    "diff_subTAD",
    "same_subTAD"
), ".png")

sub3TADfigPaths <- paste0("data/figs/gene_pairs_", c(
    "not_same_TAD",
    "diff_subTAD",
    "same_subTAD"
), ".png")

#-----------------------------------------------------------------------
# Rao et al. 2014 Hi-C sample
#-----------------------------------------------------------------------
CELL <- "IMR90"
#~ HIC_RESOLUTION <- 50*10^3 # 50kb
HIC_RESOLUTION <- 5*10^3 # 5kb
HIC_DATA_DIR <- "data/Rao2014"

#-----------------------------------------------------------------------
# Colors used for plotting
#-----------------------------------------------------------------------
# two colors for paralog and non-paralog genes
COL=brewer.pal(9, "Set1")[c(1,9)]   # for paralog vs. sampled genes
COL_RAND=c(COL[1], brewer.pal(8, "Dark2")[8])   # random genes
COL2=brewer.pal(9, "Set1")[c(1,2)]  # for paralog vs. non-paralog
COL_DOMAIN=brewer.pal(12, "Set3")
COL_FAMILY=brewer.pal(8, "Dark2")
COL_EH = brewer.pal(9, "Set1")[5]
COL_TAD = brewer.pal(8, "Set1")[c(3,5)]
COL_ORTHO = c(brewer.pal(12, "Paired")[5], brewer.pal(8, "Pastel2")[8])
COL_SPECIES = brewer.pal(8, "Accent")[1:2]
COL_STRAND = brewer.pal(8, "Pastel1")[2:1]
COL_AGE = brewer.pal(9, "YlOrRd")[c(6,9)]
COL_EH_POS=brewer.pal(9, "Set3")[c(5,9,6)]
COL_COMP=c("#92AC51","#F2ED98","#B6534D", "#DFDFDF") # AA, BB, AB, <NA>

COL_COMP_REG=brewer.pal(12, "Paired")[c(8,2,9,7,1,11)]
names(COL_COMP_REG) <- c("A/A_REG", "B/B_REG", "DIF_REG", "AA", "BB", "A/B")
#-----------------------------------------------------------------------
# Parameters for data loading from image files
#-----------------------------------------------------------------------

# use local data or downlaod data from ensemble
USE_LOCAL_HIC = TRUE
LOAD_INPUT_DATA=TRUE
LOAD_PAIRS=TRUE

#-----------------------------------------------------------------------
# Parameters critical to the analysis
#-----------------------------------------------------------------------

# number of random permutations of whole data sets
N_RAND=10            
SELECT_OLD_PAIRS=FALSE

# number of bins for distance density estimation as parameter to adjust sampling resolution
RANDOM_SEED=13521

MAX_DIST=10^6
DISTAL_MAX_DIST=10^9
DISTAL_MIN_DIST=10^6

DIST_TH=c(10^6, 10^5)

HIPPIE_MEDIUM=0.63
orgStr2Name = c(mmusculus="mouse", cfamiliaris="dog")

#-----------------------------------------------------------------------
# Version of the analysis and output file paths
#-----------------------------------------------------------------------

VERSION="v16"

outDataPrefix = "results/paralog_regulation/EnsemblGRCh37_paralog_genes"

# if SELECT_OLD_PAIRS ture, indicate this in the prefix!
if (SELECT_OLD_PAIRS){
    outPrefix = paste0("results/paralog_regulation/", VERSION, "_maxDist_", MAX_DIST, "_nrand_", N_RAND, "_HiC-res_", HIC_RESOLUTION, "_SelectOldPairs_", SELECT_OLD_PAIRS, "/EnsemblGRCh37_paralog_genes")
}else{
    outPrefix = paste0("results/paralog_regulation/", VERSION, "_maxDist_", MAX_DIST, "_nrand_", N_RAND, "_HiC-res_", HIC_RESOLUTION, "/EnsemblGRCh37_paralog_genes")
}

# create directory, if not exist
dir.create(dirname(outPrefix), recursive=TRUE, showWarnings = FALSE)

# image file to save temporary and final R session with all data
WORKIMAGE_FILE = paste0(outPrefix, "workspace.Rdata")

