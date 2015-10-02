#!/usr/bin/Rscript
#=======================================================================
#
#   parse baseline expression data sets from EBI Expression atlas
#
#=======================================================================

# load some useful libraries
require(stringr)        # for functions like paste0()

#=======================================================================
# parameters and input data files
#=======================================================================
VERSION_DATA_EXPRESSION="v02"

# Human gene expression data sets from Expression Atlas:
EXPRESSION_FILE_LIST = list(
    "ENCODE_cell_lines"="data/ExpressionAtlas/E-GEOD-26284.tsv",
    "GTEx"="data/ExpressionAtlas/E-MTAB-2919.tsv",
    "FANTOM5_tissues"="data/ExpressionAtlas/E-MTAB-3358.tsv",
    "Illumina_Body_Map"="data/ExpressionAtlas/E-MTAB-513.tsv")

#~     "Uhlen_Lab_tissues"="data/ExpressionAtlas/E-MTAB-2836.tsv",

outPrefixDataExpression=paste0("results/data/expression.ebi_expression_atlas.", VERSION_DATA_EXPRESSION)

# make directory if not exist already
dir.create(dirname(outPrefixDataExpression), showWarnings = FALSE)

#=======================================================================

#-----------------------------------------------------------------------
# parse the base line expression data sets from EBI's Expression Atlas files
#-----------------------------------------------------------------------
parseExpressionAtlas = function(inFile){
    ge = read.delim(inFile, comment.char='#')
    names(ge) = gsub("\\.", "_", names(ge))
    row.names(ge) = ge[,1]
    tissues = names(ge)[3:ncol(ge)]

    # get only the columns with expression values of all tissues
    expDF = ge[3:ncol(ge)]
}


#=======================================================================
# parse baseline gene expression data set:
#=======================================================================

# parse expression data as list of data.frame
expDFlist = lapply(EXPRESSION_FILE_LIST, parseExpressionAtlas)

