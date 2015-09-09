#!/bin/bash

#=======================================================================
# this script is suppost do document all downlaods in the data folder
#=======================================================================


#=======================================================================
# Hi-C data from Rao et al 2014 Cell
# (see: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525)
#=======================================================================

# IMR90 contact matrices
mkdir -p Rao2014 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525%5FIMR90%5Fcontact%5Fmatrices%2Etar%2Egz
# IMR90 domain list:
#wget -P Rao2014 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525%5FIMR90%5FHiCCUPS%5Flooplist%2Etxt%2Egz

# GM12878 combined contact matrices:
wget -P Rao2014 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525%5FGM12878%5Fcombined%5Fcontact%5Fmatrices%2Etar%2Egz

# unzip all files
gunzip Rao2014/*.gz



