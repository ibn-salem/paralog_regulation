#!/usr/bin/Rscript
#=======================================================================
#
#   Parameters for TAD_gene_distribution.R 
# 
#=======================================================================


#=======================================================================
# parameters and input data files
#=======================================================================
VERSION="v02"

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

# pseudogenes data
PseudoGenesFiles = c(
    "Psg_proc_TS"="data/PseudoGenes/Processed_Psg_Transcribed_L1.txt",
    "Psg_proc_non-TS"="data/PseudoGenes/Processed_Psg_non-Transcribed_PL1_2_3.txt",
    "Psg_dup_TS"="data/PseudoGenes/Duplicated_Psg_Transcribed_L1.txt",
    "Psg_dup_non-TS"="data/PseudoGenes/Duplicated_Psg_non-Transcribed_PL1_2_3.txt"
)


outPrefix=paste0("results/TAD_gene_distribution/", VERSION, ".TAD_gene_distribution")
WORKIMAGE_FILE=paste0(outPrefix, ".workspace.Rdata")

# make directory if not exist already
dir.create(dirname(outPrefix), showWarnings = FALSE)

COL_REGIONS = brewer.pal(8, "Accent")
COL_PAIRED = brewer.pal(12, "Paired")
COL_PAIRED_1 = COL_PAIRED[seq(1,12,2)]
COL_PAIRED_2 = COL_PAIRED[seq(2,12,2)]

COL_TEST_CTL = brewer.pal(9, "Set1")[c(1,9)]   # for test vs. control comp.
COL_DARK2 = brewer.pal(8, "Dark2")

N_BINS=4
N_BIN_LIST=c(4, 40)

RANDOM_SEED=13521


WINDOW_SIZE = 400*10^3
