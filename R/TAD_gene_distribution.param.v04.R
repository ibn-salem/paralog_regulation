#!/usr/bin/Rscript
#=======================================================================
#
#   Parameters for TAD_gene_distribution.R 
# 
#=======================================================================


#=======================================================================
# parameters and input data files
#=======================================================================
VERSION="v04"

chromFile = "data/hg19/hg19.genome"

RaoDomainFiles = c(
    Rao_HeLa="data/Rao2014/GSE63525_HeLa_Arrowhead_domainlist.txt",
    Rao_HUVEC="data/Rao2014/GSE63525_HUVEC_Arrowhead_domainlist.txt",
    Rao_K562="data/Rao2014/GSE63525_K562_Arrowhead_domainlist.txt",
    Rao_KBM7="data/Rao2014/GSE63525_KBM7_Arrowhead_domainlist.txt",
    Rao_NHEK="data/Rao2014/GSE63525_NHEK_Arrowhead_domainlist.txt",
    Rao_IMR90="data/Rao2014/GSE63525_IMR90_Arrowhead_domainlist.txt",
    Rao_GM12878="data/Rao2014/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt"
)

RaoSubcompartmentFile="data/Rao2014/GSE63525_GM12878_subcompartments.bed"

DixonDomainFiles = c(
    Dixon_hESC="data/Dixon2012/hESC.hg18.bed.hg19.bed",
    Dixon_IMR90="data/Dixon2012/IMR90.hg18.bed.hg19.bed"
)

# pseudogenes data
PseudoGenesFiles = c(
    "Psg_proc_TS"="data/PseudoGenes/Processed_Psg_Transcribed_L1.txt",
    "Psg_proc_non-TS"="data/PseudoGenes/Processed_Psg_non-Transcribed_PL1_2_3.txt",
    "Psg_dup_TS"="data/PseudoGenes/Duplicated_Psg_Transcribed_L1.txt",
    "Psg_dup_non-TS"="data/PseudoGenes/Duplicated_Psg_non-Transcribed_PL1_2_3.txt"
)

# essential genes
OGEEdbFile="data/OGEEdb/9606_dataset348.txt"


COL_REGIONS = brewer.pal(8, "Accent")
COL_PAIRED = brewer.pal(12, "Paired")
COL_PAIRED_1 = COL_PAIRED[seq(1,12,2)]
COL_PAIRED_2 = COL_PAIRED[seq(2,12,2)]

COL_TEST_CTL = brewer.pal(9, "Set1")[c(1,9)]   # for test vs. control comp.
COL_DARK2 = brewer.pal(8, "Dark2")
COL_STRAND = brewer.pal(8, "Pastel1")[2:1]
COL_DOMAIN=brewer.pal(12, "Set3")

N_BINS=4
N_BIN_LIST=c(4, 40)

RANDOM_SEED=13521

WINDOW_SIZE = 400*10^3

# take the center of a gene for overlap calculations
GENE_COORDINATE_CENTER=TRUE

outPrefix=paste0("results/TAD_gene_distribution/", VERSION, ".TAD_gene_distribution.nonOverlppingTAD")
WORKIMAGE_FILE=paste0(outPrefix, ".workspace.Rdata")

# make directory if not exist already
dir.create(dirname(outPrefix), showWarnings = FALSE)

