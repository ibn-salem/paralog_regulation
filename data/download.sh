#!/bin/bash

#=======================================================================
# this script is suppost do document all downlaods in the data folder
#=======================================================================

# set some variables here:
BIN=../bin
mkdir -p ${BIN}

#=======================================================================
# General genome assembly based data and tools from UCSC:
#=======================================================================

#-----------------------------------------------------------------------
# UCSC chrom size tables
#-----------------------------------------------------------------------

# human hg19
mkdir -p hg19
#~ wget -P hg19 http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes 
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e "select chrom, size from hg19.chromInfo"  > hg19/hg19.genome
# filter for real chromosomes:
head -n 25  hg19/hg19.genome > hg19/hg19.genome.realChroms

# mouse genome chromosome sizes
mkdir -p mm10
wget -P mm10 http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes

# dog genome chromosome sizes
mkdir -p canFam3
wget -P canFam3 http://hgdownload.cse.ucsc.edu/goldenPath/canFam3/bigZips/canFam3.chrom.sizes

#-----------------------------------------------------------------------
# UCSC liftover data and tool
#-----------------------------------------------------------------------

# UCSC liftover chains
mkdir -p UCSC
wget -P UCSC http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg18.over.chain.gz
wget -P UCSC http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz
wget -P UCSC http://hgdownload.cse.ucsc.edu/goldenPath/mm9/liftOver/mm9ToMm10.over.chain.gz
gunzip UCSC/*.gz

# download liftOver tool from UCSC:
wget -P ${BIN} http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
chmod u+x ${BIN}/liftOver

#=======================================================================
# specific data sets
#=======================================================================

#-----------------------------------------------------------------------
# Anderson2014 interaction map based on correlation of CAGE data
#-----------------------------------------------------------------------
mkdir -P Andersson2014
wget -P Andersson2014 http://enhancer.binf.ku.dk/presets/enhancer_tss_associations.bed
wget -P Andersson2014 http://enhancer.binf.ku.dk/presets/robust_enhancers.bed
wget -P Andersson2014 http://enhancer.binf.ku.dk/presets/permissive_enhancers.bed

#=======================================================================
# Hi-C data from Rao et al 2014 Cell
#=======================================================================
mkdir -p Rao2014

# IMR90 combined matrices (~3.8 GB in .gz)
wget -P Rao2014 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_IMR90_intrachromosomal_contact_matrices.tar.gz

# unpack matrices
tar xvfz Rao2014/*.tar.gz -C Rao2014

# subcompartments in GM12878 
wget -P Rao2014 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525%5FGM12878%5Fsubcompartments%2Ebed%2Egz

# unzip all
gunzip Rao2014/*bed.gz


#~ RAO_CELLS="GM12878_primary+replicate HMEC HUVEC HeLa IMR90 K562 KBM7 NHEK"
RAO_CELLS="GM12878_primary+replicate HMEC HUVEC HeLa IMR90 K562 KBM7 NHEK CH12-LX"


for CELL in ${RAO_CELLS} ; do
    # download
    wget -P Rao2014 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_${CELL}_Arrowhead_domainlist.txt.gz

    # unzip 
    gunzip Rao2014/GSE63525_${CELL}_Arrowhead_domainlist.txt.gz
        
    # re-format TADs into bed file
	if [ "$CELL" == "CH12-LX" ] 
	then		
	    tail -n +2 Rao2014/GSE63525_${CELL}_Arrowhead_domainlist.txt \
			|cut -f 1-3 \
			> Rao2014/GSE63525_${CELL}_Arrowhead_domainlist.txt.bed 
	else
	    tail -n +2 Rao2014/GSE63525_${CELL}_Arrowhead_domainlist.txt \
			|cut -f 1-3 \
			| sed -e 's/^/chr/' \
			> Rao2014/GSE63525_${CELL}_Arrowhead_domainlist.txt.bed 
	fi
	
done

# convert mouse domains from mm9 to mm10 assembly
${BIN}/liftOver \
	Rao2014/GSE63525_CH12-LX_Arrowhead_domainlist.txt.bed  \
	UCSC/mm9ToMm10.over.chain \
	Rao2014/GSE63525_CH12-LX_Arrowhead_domainlist.txt.bed.mm9.bed.mm10.bed \
	Rao2014/GSE63525_CH12-LX_Arrowhead_domainlist.txt.bed.mm9.bed_unmapped.bed


#=======================================================================
# HIPPIE and expression data
#=======================================================================
mkdir -p HIPPIE

wget -P HIPPIE http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/hippie_current.txt


#=======================================================================
# Hi-C data form Dixon et al. 2012
#=======================================================================
mkdir -p Dixon2012

DIXON_CELLS="hESC IMR90"

for C in ${DIXON_CELLS} ; do
    
    # download
    #~ wget -P Dixon2012 http:///chromosome.sdsc.edu/mouse/hi-c/${C}.domain.tar.gz
    wget -P Dixon2012 http://132.239.201.216/mouse/hi-c/${C}.domain.tar.gz    
    # extract and rename
    tar xvfz Dixon2012/${C}.domain.tar.gz -C Dixon2012
    cp Dixon2012/${C}/combined/total.combined.domain Dixon2012/${C}.hg18.bed
    
    # liftover to hg19
	${BIN}/liftOver \
        Dixon2012/${C}.hg18.bed \
        UCSC/hg18ToHg19.over.chain \
        Dixon2012/${C}.hg18.bed.hg19.bed \
        Dixon2012/${C}.hg18.bed.hg19.bed_unmapped.bed
done

#http://132.239.201.216/mouse/hi-c/mESC.domain.tar.gz
# http://132.239.201.216/mouse/hi-c/cortex.domain.tar.gz
DIXON_MOUSE="mESC cortex"
for C in ${DIXON_MOUSE} ; do

    # download
    wget -P Dixon2012 http://132.239.201.216/mouse/hi-c/${C}.domain.tar.gz    
    # extract and rename
    tar xvfz Dixon2012/${C}.domain.tar.gz -C Dixon2012
    cp Dixon2012/${C}/*combined/total.*combined.domain Dixon2012/mouse.${C}.mm9.bed

    # liftover to mm10
	${BIN}/liftOver \
        Dixon2012/mouse.${C}.mm9.bed \
        UCSC/mm9ToMm10.over.chain \
        Dixon2012/mouse.${C}.mm9.bed.mm10.bed \
        Dixon2012/mouse.${C}.mm9.bed_unmapped.bed


done



#=======================================================================
# EBI Expression Atlas 
#=======================================================================
mkdir -p ExpressionAtlas

# downloaded at 31.03.2015
# Download some specific data sets:

# FANTOM5 tissue table (Tissues - 68 RIKEN FANTOM5 project (CAGE))
wget "http://www.ebi.ac.uk/gxa/experiments/E-MTAB-3358.tsv?accessKey=&geneQuery=&cutoff=-0.1" -O ExpressionAtlas/E-MTAB-3358.tsv

# Cell Lines - ENCODE
wget "http://www.ebi.ac.uk/gxa/experiments/E-GEOD-26284.tsv?accessKey=&geneQuery=&cutoff=-0.1" -O ExpressionAtlas/E-GEOD-26284.tsv

# Tissues - Illumina Body Map
wget "http://www.ebi.ac.uk/gxa/experiments/E-MTAB-513.tsv?accessKey=&geneQuery=&cutoff=-0.1" -O ExpressionAtlas/E-MTAB-513.tsv

# Tissues - GTEx
wget "http://www.ebi.ac.uk/gxa/experiments/E-MTAB-2919.tsv?accessKey=&geneQuery=&cutoff=-0.1" -O ExpressionAtlas/E-MTAB-2919.tsv


#=======================================================================
# House keeping genes from Eisenberg 2013 et al. (http://www.ncbi.nlm.nih.gov/pubmed/23810203)
#=======================================================================
mkdir -p Eisenberg2013
wget -P Eisenberg2013 http://www.tau.ac.il/~elieis/HKG/HK_genes.txt


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

gunzip Rudan2015/*.gz

#=======================================================================
# Capture Hi-C data from Mifsud2015
#=======================================================================
mkdir -p Mifsud2015
wget -P Mifsud2015 http://www.ebi.ac.uk/arrayexpress/files/E-MTAB-2323/E-MTAB-2323.additional.1.zip
unzip Mifsud2015/E-MTAB-2323.additional.1.zip -d Mifsud2015

#=======================================================================
# Cohesin KO expression data from Dimitris Polychronopoulos <d.polychronopoulos@imperial.ac.uk>
#=======================================================================
# copied file manually to 
# Dimitirs/RNAseqWTRad21KOMacrophages.RData


