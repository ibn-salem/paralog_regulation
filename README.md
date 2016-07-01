# paralog_regulation

This repository contains the source code for all analysis in the study Ibn-Salem et al. 2016, "Co-regulation of paralog genes in the three-dimensional chromatin architecture" (in revision).

The following documentation guides through all steps including downloading of external source data, retrieving paralog genes from ENSEMBL, filtering, annotation, and running all analysis.


### Requirements:

The analysis is mainly implemented in R by using several additional packages and depend on:

 - R version 3.2.2 (2015-08-14)
 - R packages:
     - biomaRt
     - BiocParallel
     - Biostrings
     - BSgenome.Hsapiens.UCSC.hg19
     - colorspace
     - data.table
     - entropy
     - gdata
     - GenomicAlignments
     - GenomicRanges
     - ggplot2
     - gplots
     - grid
     - gridExtra
     - gtable
     - HiTC
     - igraph
     - Matrix
     - minerva
     - plyr
     - png
     - RColorBrewer
     - rPython
     - rtracklayer
     - scales
     - stringr
     - TFBSTools
     - xtable


 - Python 2.7 
 - Python packages: 
 	- networkx (https://networkx.github.io/)


## Retrieving public data

Downloading of all external data is documented in the file [`data/download.sh`](data/download.sh). It is a bash script which can be executed in the `data` folder to download all required external data.

The script can be executed as follows:
```bash
cd data
./download.sh
```

## Downloading of ENSEMBL data

The next step after retiring external annotation data is to download paralog gene pairs from ENSEMBL database from within R using the `biomaRt` package. This is implemented in the R script [`R/data.ensembl.R`](R/data.ensembl.R)


Similar the parsing of expression and Capture Hi-C data is implemented in the scripts [`R/data.expression.R`](R/data.expression.R) and [`R/data.captureHiC.R`](R/data.captureHiC.R), respectively. 

These scripts are executed automatically within the next step.


## Filtering sampling and annotation


From the project root folder start an R session and run
```
source("R/paralog_regulation.sampling_annotation.R")
```

This will download and parse data, filter gene pairs, sample random background control gene pairs, and annotate all pairs. 



## Main analysis

The main analyis after filtering, sampling and annotation is implemented in the R scripts [`R/paralog_regulation.analysis.R`](R/paralog_regulation.analysis.R) and [`R/paralog_regulation.analysis_allDF.R`](R/paralog_regulation.analysis_allDF.R).

From within R execute the following commands:

```
source("R/paralog_regulation.analysis.R")
source("R/paralog_regulation.analysis_allDF.R")
```

This will create a folder `results` with several plots in PDF format.

## Parameters and temporary files

The script [`R/paralog_regulation.param.v16.R`](R/paralog_regulation.param.v16.R) defines several parameters like file paths, colors, distance thresholds, and number sampling replications, which can be modified if needed.


