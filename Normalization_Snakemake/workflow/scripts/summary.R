library(dplyr)
library(tidyr)
library(optparse)
option.list <- list(
    make_option("--summarySheets", type="character"),
    make_option("--samples", type="character"),
    make_option("--outputFile", type="character"),
    make_option("--ATAC", type="logical")
)
opt <- parse_args(OptionParser(option_list=option.list))
sheets= opt$summarySheets %>% strsplit(split=",") %>% unlist()
samples = opt$samples %>% strsplit(split=",") %>% unlist()

#######################################################
getSample <- function(path){
        pathlist=unlist(strsplit(path, "/"))
        sample=pathlist[length(pathlist)-3]
        return(sample)
}
#######################################################

if (!opt$ATAC){
    output = as.data.frame(matrix(nrow=length(sheets),ncol=2))
    rownames(output) = samples
    colnames(output)=c("cellnumber_before_filtering", "cellnumber_after_filtering")
    for (sheet in sheets){
	    df = read.table(sheet, sep="\t", header=T, row.names=1, stringsAsFactors=F)
	    sample = getSample(sheet)
	    output[sample, "cellnumber_before_filtering"]=df[sample, "cellnumber_before_filtering"]
	    output[sample, "cellnumber_after_filtering"]=df[sample, "cellnumber_after_filtering"]
    }
} else {
    output = as.data.frame(matrix(nrow=length(sheets),ncol=3))
    rownames(output) = samples
    colnames(output)=c("cellnumber_before_filtering", "multiomic_cellnumber_before_filtering", "cellnumber_after_filtering")
    for (sheet in sheets){
	    print(sheet)
	    df = read.table(sheet, sep="\t", header=T, row.names=1, stringsAsFactors=F)
	    sample = getSample(sheet)
	    output[sample, "cellnumber_before_filtering"]=df[sample, "cellnumber_before_filtering"]
	    if ("cellnumber_before_filtering_after_ATAC_filtering" %in% colnames(df)) {
		output[sample, "multiomic_cellnumber_before_filtering"]=df[sample, "cellnumber_before_filtering_after_ATAC_filtering"]
            }
	    output[sample, "cellnumber_after_filtering"]=df[sample, "cellnumber_after_filtering"]
    }
}

write.table(output, opt$outputFile, sep="\t", row.names=T, col.names=NA, quote=FALSE)

