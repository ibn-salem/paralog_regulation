########################################################################
#
# This script contains parameters for the paralog_regulation.R analysis.
#
########################################################################

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

#-----------------------------------------------------------------------
# Rao et al. 2014 Hi-C sample
#-----------------------------------------------------------------------
CELL <- "IMR90"
HIC_RESOLUTION <- 50*10^3 # 50kb
#~ HIC_RESOLUTION <- 5*10^3 # 5kb
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
#~ COL_COMP=c(AA="#92AC51", AB="#F2ED98", BB="#B6534D")
COL_COMP=c("#92AC51","#F2ED98","#B6534D", "#DFDFDF") # AA, AB, BB, <NA>
#-----------------------------------------------------------------------
# Parameters for data loading
#-----------------------------------------------------------------------

# use local data or downlaod data from ensemble
USE_LOCAL = FALSE
USE_LOCAL_HIC_CONTACTS = FALSE
#~ N_CPUS=20

#-----------------------------------------------------------------------
# Parameters critical to the analysis
#-----------------------------------------------------------------------

# number of random permutations of whole data sets
N_RAND=10            

# number of bins for distance density estimation as parameter to adjust sampling resolution
N_SAMPLING_BIN <- 20
RANDOM_SEED=13521

MAX_DIST=10^6
DISTAL_MAX_DIST=10^9
DISTAL_MIN_DIST=10^6

DIST_TH=c(10^6, 10^5)

#-----------------------------------------------------------------------
# Version of the analysis and output file paths
#-----------------------------------------------------------------------

VERSION="v15"

outDataPrefix = "results/paralog_regulation/EnsemblGRCh37_paralog_genes"
outPrefix = paste0("results/paralog_regulation/", VERSION, "_maxDist_", MAX_DIST, "_nrand_", N_RAND, "_HiC-res_", HIC_RESOLUTION, "/EnsemblGRCh37_paralog_genes")

# create directory, if not exist
dir.create(dirname(outPrefix), recursive=TRUE, showWarnings = FALSE)

# image file to save temporary and final R session with all data
WORKIMAGE_FILE = paste0(outPrefix, "workspace.Rdata")

