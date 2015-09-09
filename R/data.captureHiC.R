#!/usr/bin/Rscript
#=======================================================================
#
#    parse capture Hi-C data between promoters from Mifsud et al. 2015
#
#=======================================================================

# load some useful libraries
require(stringr)        # for functions like paste0()

# depend on the following scripts
# load costume scripts 
#source("R/parseHiC.R")
#source("R/data.ensembl.R")

#=======================================================================
# parameters and input data files
#=======================================================================
VERSION_DATA_CAPTUREHIC="v01"

#promoter-promoter interaction from Capture Hi-C (Mifsud2015a)
CAPTURC_FILE="data/Mifsud2015/TS5_GM12878_promoter-promoter_significant_interactions.txt"

outPrefixDataCaptureC=paste0("results/data/captureC.Mifsud2015.", VERSION_DATA_CAPTUREHIC)

USE_LOCAL_DATA_CAPTUREHIC = TRUE


# make directory if not exist already
dir.create(dirname(outPrefixDataCaptureC), showWarnings = FALSE)


#=======================================================================
# parse captureHiC data as matrix of raw counts and obs/exp values:
#=======================================================================



if ( !USE_LOCAL_DATA_CAPTUREHIC) {    

    # parse data and save it as R object

    # genes and seqInfoRealChrom should be loaded by data.ensembl.R script
    tssGR = getTssGRfromENSEMBLGenes(genes, seqInfoRealChrom, colNames=c("hgnc_symbol"))
    
    captureHiC <- parseCaptureHiC(inFile=CAPTURC_FILE, tssGR)    

#~     # due to sparse matrix data structure non available pairs will get 0 counts
#~     # since no 0 count pair is in the original data, we can replace all 0 with NA
#~     zeroSub = captureHiC[["raw"]] == 0    
#~     captureHiC[["raw"]][zeroSub] = NA
#~     captureHiC[["obsExp"]][zeroSub] = NA
    
    save(captureHiC, file=paste0(outPrefixDataCaptureC, ".captureHiC.Rdata"))

}else{
    # just load the R object
    load(paste0(outPrefixDataCaptureC, ".captureHiC.Rdata"))
}



