#!/bin/bash

#=======================================================================
# this script is suppost do document all downlaods in the data folder
#=======================================================================

# set some variables here:
BIN=../bin
mkdir -p ${BIN}

#-----------------------------------------------------------------------
# General genome assembly based data and tools from UCSC:
#-----------------------------------------------------------------------

# UCSC chrom size tables
mkdir -p hg19
#~ wget -P hg19 http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes 
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from hg19.chromInfo"  > hg19/hg19.genome
# filter for real chromosomes:
head -n 25  hg19/hg19.genome > hg19/hg19.genome.realChroms

# get hg19 reference genome sequence
wget -P hg19 http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit
# convert .2bit to .fasta format
wget -P ${BIN} http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa
chmod u+x ${BIN}/twoBitToFa
${BIN}/twoBitToFa hg19/hg19.2bit hg19/hg19.fa

# mouse genome chromosome sizes
mkdir -p mm10
wget -P mm10 http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes

# dog genome chromosome sizes
mkdir -p canFam3
wget -P canFam3 http://hgdownload.cse.ucsc.edu/goldenPath/canFam3/bigZips/canFam3.chrom.sizes

# UCSC liftover chains
mkdir -p UCSC
wget -P UCSC http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg18.over.chain.gz
wget -P UCSC http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz
gunzip UCSC/*.gz

# download liftOver tool from UCSC:
wget -P ${BIN} http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
chmod u+x ${BIN}/liftOver

#-----------------------------------------------------------------------
# specific data sets
#-----------------------------------------------------------------------

# Anderson2014 interaction map based on correlation of CAGE data
mkdir -P Andersson2014
wget -P Andersson2014 http://enhancer.binf.ku.dk/presets/enhancer_tss_associations.bed
wget -P Andersson2014 http://enhancer.binf.ku.dk/presets/robust_enhancers.bed
wget -P Andersson2014 http://enhancer.binf.ku.dk/presets/permissive_enhancers.bed



# alternative topological domains from Filippova et al. 2014
#~ mkdir -p Filippova2014
#~ wget -P Filippova2014 http://www.cs.cmu.edu/~ckingsf/software/armatus/IMR90.consensus.tar
#~ wget -P Filippova2014 http://www.cs.cmu.edu/~ckingsf/software/armatus/mESC.consensus.tar
#~ cd Filippova2014/
#~ tar xf IMR90.consensus.tar
#~ tar xf mESC.consensus.tar
#~ cd ..

# Association between genes and tissues:
#~ mkdir -p TISSUES
#~ wget -P TISSUES http://download.jensenlab.org/human_tissue_experiments_full.tsv
#~ wget -P TISSUES http://download.jensenlab.org/human_tissue_experiments_filtered.tsv


# manual downlaods: of 1000 Genome deletion eQTL loci 
# Email by Sascha Meyers from 30.01.15
# save the file here: /home/ibnsalem/PhD/data/1kG/final.bed
# remove header : /home/ibnsalem/PhD/data/1kG/final.bed.noHeader
# use only the firts three columns (deletion coordinates): /home/ibnsalem/PhD/data/1kG/final.bed.noHeader.del

# reformat 1kG eQTL deletions to interaction file format
python ../python/bed2LongRangeFormat.py -i 1kG/final.bed.noHeader -e Ensembl/GRCh37_human_genes_ENSID_chr_txStart.txt.noHeader.uniqENSG.bed -t -o 1kG/final.bed.noHeader.LongRangeTextFormat


#=======================================================================
# ENSEMBL data for ENSEMBL (GRCh37):
#=======================================================================

# installing biomat-perl to fetch data automatically:
# Following this instructions: http://www.biomart.org/other/install-overview.html
    
#=======================================================================
# ENSEMBL manual downloads:
#=======================================================================

# 30.01.15
# manualy download ENSEMBL Genes with ENSEMBLE GENE ID, chr, TSS
# save file here: /home/ibnsalem/PhD/data/Ensembl/GRCh38_human_genes_ENSID_chr_TSS.txt
# rm header:
#~ tail -n +2 GRCh38_human_genes_ENSID_chr_TSS.txt > GRCh38_human_genes_ENSID_chr_TSS.txt.noHeader
# reformat into BED
#~ awk '{print "chr"$2"\t"$3"\t"$3+1"\t"$1}' GRCh38_human_genes_ENSID_chr_TSS.txt.noHeader > GRCh38_human_genes_ENSID_chr_TSS.txt.noHeader.bed

# take only one unique ENS Gene ID
#~ sort -u  -k 4 GRCh38_human_genes_ENSID_chr_TSS.txt.noHeader.bed > GRCh38_human_genes_ENSID_chr_TSS.txt.noHeader.uniqENSG.bed

# manualy download ENSEMBL Genes (Ensembl GRCh37 release 78) with ENSEMBLE GENE ID, chr, TSS
# save file here: /home/ibnsalem/PhD/data/Ensembl/GRCh37_human_genes_ENSID_chr_txStart.txt
# rm header:
#~ tail -n +2 Ensembl/GRCh37_human_genes_ENSID_chr_txStart.txt > Ensembl/GRCh37_human_genes_ENSID_chr_txStart.txt.noHeader
# reformat into BED
#~ awk '{print "chr"$2"\t"$3"\t"$3+1"\t"$1}' Ensembl/GRCh37_human_genes_ENSID_chr_txStart.txt.noHeader > Ensembl/GRCh37_human_genes_ENSID_chr_txStart.txt.noHeader.bed
# unique by ENS Gene ID
#~ sort -u  -k 4 Ensembl/GRCh37_human_genes_ENSID_chr_txStart.txt.noHeader.bed > Ensembl/GRCh37_human_genes_ENSID_chr_txStart.txt.noHeader.uniqENSG.bed



#=======================================================================
# Hi-C data from Rao et al 2014 Cell
#=======================================================================
mkdir -p Rao2014
wget -P Rao2014 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525%5FIMR90%5FREADME%2Ertf

# IMR90 combined matrices (~3.8 GB in .gz)
wget -P Rao2014 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_IMR90_intrachromosomal_contact_matrices.tar.gz

# GM12878 combined matrices (~35 GB in .gz)
wget -P Rao2014 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_GM12878_combined_intrachromosomal_contact_matrices.tar.gz

# mouse CH12-LX 
wget -P Rao2014 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525%5FCH12%2DLX%5FArrowhead%5Fdomainlist%2Etxt%2Egz

# unpack matrices
tar xvfz Rao2014/*.tar.gz -C Rao2014

# unzip all
gunzip Rao2014/*.gz


RAO_CELLS="GM12878_primary+replicate HMEC HUVEC HeLa IMR90 K562 KBM7 NHEK"

for CELL in ${RAO_CELLS} ; do
    # download
    wget -P Rao2014 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_${CELL}_Arrowhead_domainlist.txt.gz
    wget -P Rao2014 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_${CELL}_HiCCUPS_looplist.txt.gz
    wget -P Rao2014 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_${CELL}_HiCCUPS_looplist_with_motifs.txt.gz

    # unzip 
    gunzip Rao2014/GSE63525_${CELL}_Arrowhead_domainlist.txt.gz
    gunzip Rao2014/GSE63525_${CELL}_HiCCUPS_looplist.txt.gz
    gunzip Rao2014/GSE63525_${CELL}_HiCCUPS_looplist_with_motifs.txt.gz
        
    # re-format TADs into bed file
    tail -n +2 Rao2014/GSE63525_${CELL}_Arrowhead_domainlist.txt |cut -f 1-3 | sed -e 's/^/chr/' > Rao2014/GSE63525_${CELL}_Arrowhead_domainlist.txt.bed 
done

# make zip for Sweta:
#zip Rao2014/Rao2014_HiC_domains.zip Rao2014/*domainlist*

#=======================================================================
# HIPPIE and expression data
#=======================================================================
mkdir -P HIPPIE
# E-Mail by Miguel from 9.03.15 "gene expression in human tissues"
# save file as 'HIPPIE/final_expression_matrix_v2.txt'


#=======================================================================
# Hi-C data form Dixon et al. 2012
#=======================================================================
mkdir -p Dixon2012


# copy the toplogical domains from Hi-C data in hESCs (Dixon et al. 2012)
mkdir -p Dixon2012

DIXON_CELLS="hESC IMR90"

for C in ${DIXON_CELLS} ; do
    
    # download
    wget -P Dixon2012 http:///chromosome.sdsc.edu/mouse/hi-c/${C}.domain.tar.gz
    
    # extract and rename
    tar xvfz Dixon2012/${C}.domain.tar.gz -C Dixon2012
    cp Dixon2012/${C}/combined/total.combined.domain Dixon2012/${C}.hg18.bed
    
    # liftover to hg19
	${BIN}/liftOver \
        Dixon2012/${C}.hg18.bed \
        UCSC/hg18ToHg19.over.chain \
        Dixon2012/${C}.hg19.bed \
        Dixon2012/${C}.hg19.bed_unmapped.bed

done


#=======================================================================
# Hi-C data form Dixon et al. 2015
#=======================================================================

mkdir -p Dixon2015 
#~ wget -P Dixon2015  ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE52nnn/GSE52457/suppl/GSE52457%5FBED%5Ffiles%2Etar%2Egz

#~ tar xvfz Dixon2015/GSE52457_BED_files.tar.gz

#=======================================================================
# Enhancer promoter interaction predictions from  He et al. 2014
#=======================================================================
mkdir -p He2014
wget -P He2014 www.healthcare.uiowa.edu/labs/tan/EP_predictions.xlsx


#=======================================================================
# Conservation data in human genome
#=======================================================================
mkdir -p ConservedRegions

# copy manually from Sweta's project
# scp jibnsale@mogon.zdv.uni-mainz.de:/project/jgu-cbdm/andradeLab/miscellaneous/project/hg19_conserved_regions.bed ConservedRegions/hg19_conserved_regions.bed

#=======================================================================
# EBI Expression Atlas 
#=======================================================================
mkdir -p ExpressionAtlas
# downloaded at 31.03.2015
#~ wget -P ExpressionAtlas ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/atlas-latest-data.tar.gz

# extract data
#tar xvfz ExpressionAtlas/atlas-latest-data.tar.gz

#-------------------
# Download some specific data sets:
#-------------------
# FANTOM5 tissue table (Tissues - 68 RIKEN FANTOM5 project (CAGE))
wget "http://www.ebi.ac.uk/gxa/experiments/E-MTAB-3358.tsv?accessKey=&geneQuery=&cutoff=-0.1" -O ExpressionAtlas/E-MTAB-3358.tsv

# Cell Lines - ENCODE
wget "http://www.ebi.ac.uk/gxa/experiments/E-GEOD-26284.tsv?accessKey=&geneQuery=&cutoff=-0.1" -O ExpressionAtlas/E-GEOD-26284.tsv

# Tissues - Illumina Body Map
wget "http://www.ebi.ac.uk/gxa/experiments/E-MTAB-513.tsv?accessKey=&geneQuery=&cutoff=-0.1" -O ExpressionAtlas/E-MTAB-513.tsv


# Tissues - GTEx
wget "http://www.ebi.ac.uk/gxa/experiments/E-MTAB-2919.tsv?accessKey=&geneQuery=&cutoff=-0.1" -O ExpressionAtlas/E-MTAB-2919.tsv

#~ # Cell Lines - NCI-60 cancer (CCLE) (1437 assays)
#~ wget "http://www.ebi.ac.uk/gxa/experiments/E-MTAB-2980.tsv?accessKey=&geneQuery=&cutoff=-0.1" -O ExpressionAtlas/E-MTAB-2980.tsv

# Tissues - 32 Uhlen's Lab (200 assays)
wget "http://www.ebi.ac.uk/gxa/experiments/E-MTAB-2836.tsv?accessKey=&geneQuery=&cutoff=-0.1" -O ExpressionAtlas/E-MTAB-2836.tsv

#wget -P  ExpressionAtlas http://www.ebi.ac.uk/gxa/experiments/E-GEOD-26284.tsv
#wget -P  ExpressionAtlas http://www.ebi.ac.uk/gxa/experiments/E-MTAB-513.tsv

#=======================================================================
# Pseudogenes on MOGON from Sweta (28.04.15)
#=======================================================================
mkdir -p PseudoGenes

# copy manually from Sweta's project
#scp jibnsale@mogon.zdv.uni-mainz.de:/project/jgu-cbdm/andradeLab/miscellaneous/project/Transcribed_Psg_data/* PseudoGenes/

#=======================================================================
# Recombination map from Kong et al. 2014
#=======================================================================
mkdir -p Kong2014
#~ wget -P Kong2014 http://www.nature.com/ng/journal/v46/n1/extref/ng.2833-S2.zip
#~ unzip Kong2014/ng.2833-S2.zip -d Kong2014
#~ 
#~ RECOMB_FILES="Kong2014/DECODE_Recombination_Events_Maternal_01OCT2013.txt Kong2014/DECODE_Recombination_Events_Paternal_01OCT2013.txt"
#~ for F in ${RECOMB_FILES} ; do
    #~ 
    #~ # remove header
    #~ tail -n +2 ${F} > ${F}.hg18.bed
    #~ 
    #~ # liftover to hg19
    #~ ${BIN}/liftOver \
        #~ ${F}.hg18.bed \
        #~ UCSC/hg18ToHg19.over.chain \
        #~ ${F}.hg19.bed \
        #~ ${F}.hg19.bed_unmapped.bed
#~ done

# combine maternal and paternal recombination events:
#~ cat Kong2014/DECODE_Recombination_Events_Maternal_01OCT2013.txt.hg19.bed Kong2014/DECODE_Recombination_Events_Paternal_01OCT2013.txt.hg19.bed > Kong2014/DECODE_Recombination_Events.maternal_and_paternal.hg19.bed

#=======================================================================
# Recombination map from Kong et al. 2010 (See:http://www.decode.com/addendum/ )
#=======================================================================
#~ mkdir -p Kong2010
#~ wget -P Kong2010 http://www.decode.com/additional/sex-averaged.rmap
#~ 
#~ # refomat file: remove header and duplicate positon column
#~ 
#~ tail -n +2 Kong2010/sex-averaged.rmap \
    #~ | awk '{ print $1, $2, $2+1, $3, $4}' \
    #~ > Kong2010/sex-averaged.rmap.hg18.bed
#~ 
#~ # liftover to hg19
#~ F=Kong2010/sex-averaged.rmap
#~ ${BIN}/liftOver \
    #~ ${F}.hg18.bed \
    #~ UCSC/hg18ToHg19.over.chain \
    #~ ${F}.hg19.bed \
    #~ ${F}.hg19.bed_unmapped.bed

#=======================================================================
# TADs in mouse and dog from Rudan et al 2015
#=======================================================================
mkdir -p Rudan2015
wget -P Rudan2015 http://www.cell.com/cms/attachment/2026643989/2045457290/mmc2.xlsx

#=======================================================================
# Hi-C interactions (within 2Mb)  from Rudan et al 2015
#=======================================================================
wget -P Rudan2015 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65126/suppl/GSE65126%5FHiCseq%5FSummary%2Etxt%2Egz

wget -P Rudan2015 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65126/suppl/GSE65126%5FHiC%5Fdog%5Fliver%5Fmerged%5F50000%2Etxt%2Egz

wget -P Rudan2015 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65126/suppl/GSE65126%5FHiC%5Fmouse%5Fliver%5Fmerged%5F50000%2Etxt%2Egz

#~ wget -P Rudan2015 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65126/suppl/GSE65126%5FHiC%5Fmacaque%5Fliver%5Fmerged%5F50000%2Etxt%2Egz
#~ wget -P Rudan2015 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65126/suppl/GSE65126%5FHiC%5Frabbit%5Fliver%5Fmerged%5F50000%2Etxt%2Egz
gunzip Rudan2015/*.gz

#=======================================================================
# JASPAR vertebrate core data base (Downloaded 26.05.15)
#=======================================================================
#~ mkdir -p JASPAR
#~ wget -P JASPAR http://jaspar.genereg.net/html/DOWNLOAD/JASPAR_CORE/pfm/nonredundant/pfm_vertebrates.txt
#~ # convert matrix to transfac format using RASAT convert-matrix
#~ convert-matrix -i JASPAR/pfm_vertebrates.txt -o JASPAR/pfm_vertebrates.txt.tf -from jaspar -to transfac
#~ 

#=======================================================================
# Capture Hi-C data from Mifsud2015
#=======================================================================
mkdir -p Mifsud2015
wget -P Mifsud2015 http://www.ebi.ac.uk/arrayexpress/files/E-MTAB-2323/E-MTAB-2323.additional.1.zip
unzip Mifsud2015/E-MTAB-2323.additional.1.zip -d Mifsud2015


#=======================================================================
# Evolutionary breakpoints from Lemaitre et al 2009 BMC Genomics
#=======================================================================
mkdir -p Lemaitre2009
wget -P Lemaitre2009 http://www.biomedcentral.com/content/supplementary/1471-2164-10-335-s2.txt

#=======================================================================
# Pseudo genes
#=======================================================================
# copy pseudo gene data manually to data/PseudoGenes

#=======================================================================
# OGEEdb (essential genes)
#=======================================================================
mdir -p OGEEdb
wget -P OGEEdb http://ogeedb.embl.de/downloads/9606_dataset348.txt.gz
gunzip OGEEdb/9606_dataset348.txt.gz

#=======================================================================
# ChIP-exo data for CTCF by Rhee and Pugh 2012 (mapped and by Starick and Ibn-Salem et al. 2015)
#=======================================================================
#mkdir -p Rhee2012 
#wget -P Rhee2012 http://trg.molgen.mpg.de/Data/chip-seq-ibnsalem/bam/hg19/CTCF/all_rep.bam.sorted.bam
# manually copy files from MPI 
#mkdir -p Starick2015/CTCF
#scp -r ibnsalem@geniux.molgen.mpg.de:/project/trg-chip/Jonas_IS/results/GSM325897_HeLa-CTCF-bsites.txt.hg19/matrix-scan/MA0139.1 Starick2015/CTCF/

#=======================================================================
# ENCODE data from UCSC Genome Browser
#=======================================================================
# see https://genome.ucsc.edu/ENCODE/downloads.html
# see https://www.encodeproject.org/search/

#~ ENCODE_LINKS="
#~ http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibTfbs/wgEncodeHaibTfbsGm12878Rad21V0416101RawRep1.bigWig
#~ http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs/wgEncodeSydhTfbsGm12878Znf143166181apStdSig.bigWig
#~ http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs/wgEncodeSydhTfbsGm12878Stat3IggmusSig.bigWig
#~ "
#~ ALL_TF_LINKS="
#~ wgEncodeSydhTfbsGm12878Bhlhe40cIggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878Brca1a300IggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878Cdpsc6327IggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878CfosStdSig.bigWig
#~ wgEncodeSydhTfbsGm12878Chd1a301218aIggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878Chd2ab68301IggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878Corestsc30189IggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878Ctcfsc15914c20StdSig.bigWig
#~ wgEncodeSydhTfbsGm12878E2f4IggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878Ebf1sc137065StdSig.bigWig
#~ wgEncodeSydhTfbsGm12878Elk112771IggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878ErraIggrabSig.bigWig
#~ wgEncodeSydhTfbsGm12878Gcn5StdSig.bigWig
#~ wgEncodeSydhTfbsGm12878Ikzf1iknuclaStdSig.bigWig
#~ wgEncodeSydhTfbsGm12878InputIggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878InputIggrabSig.bigWig
#~ wgEncodeSydhTfbsGm12878InputStdSig.bigWig
#~ wgEncodeSydhTfbsGm12878InputTnfaIggrabSig.bigWig
#~ wgEncodeSydhTfbsGm12878Irf3IggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878JundIggrabSig.bigWig
#~ wgEncodeSydhTfbsGm12878JundStdSig.bigWig
#~ wgEncodeSydhTfbsGm12878MafkIggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878MaxIggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878MaxStdSig.bigWig
#~ wgEncodeSydhTfbsGm12878Mazab85725IggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878Mxi1IggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878Nfe2sc22827StdSig.bigWig
#~ wgEncodeSydhTfbsGm12878NfkbTnfaIggrabSig.bigWig
#~ wgEncodeSydhTfbsGm12878NfyaIggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878NfybIggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878Nrf1IggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878P300IggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878P300bStdSig.bigWig
#~ wgEncodeSydhTfbsGm12878P300sc584IggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878Pol2IggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878Pol2StdSig.bigWig
#~ wgEncodeSydhTfbsGm12878Pol2s2IggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878Pol3StdSig.bigWig
#~ wgEncodeSydhTfbsGm12878Rad21IggrabSig.bigWig
#~ wgEncodeSydhTfbsGm12878Rfx5200401194IggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878Sin3anb6001263IggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878Smc3ab9263IggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878Spt20StdSig.bigWig
#~ wgEncodeSydhTfbsGm12878Srebp1IggrabSig.bigWig
#~ wgEncodeSydhTfbsGm12878Srebp2IggrabSig.bigWig
#~ wgEncodeSydhTfbsGm12878Stat1StdSig.bigWig
#~ wgEncodeSydhTfbsGm12878Stat3IggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878Tblr1ab24550IggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878TbpIggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878Tr4StdSig.bigWig
#~ wgEncodeSydhTfbsGm12878Usf2IggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878WhipIggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878Yy1StdSig.bigWig
#~ wgEncodeSydhTfbsGm12878Znf143166181apStdSig.bigWig
#~ wgEncodeSydhTfbsGm12878Znf274StdSig.bigWig
#~ wgEncodeSydhTfbsGm12878Znf384hpa004051IggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878Zzz3StdSig.bigWig
#~ "
#~ 
#~ SELECTED_ENCODE="
#~ wgEncodeSydhTfbsGm12878Stat1StdSig.bigWig
#~ wgEncodeSydhTfbsGm12878Stat3IggmusSig.bigWig
#~ wgEncodeSydhTfbsGm12878Yy1StdSig.bigWig
#~ wgEncodeSydhTfbsGm12878Znf143166181apStdSig.bigWig
#~ wgEncodeSydhTfbsGm12878Znf274StdSig.bigWig
#~ wgEncodeSydhTfbsGm12878Znf384hpa004051IggmusSig.bigWig
#~ "
#~ 
#~ mkdir -p ENCODE
#~ for F in $SELECTED_ENCODE; do
    #~ wget -P ENCODE http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs/${F} 
#~ done
#~ 

