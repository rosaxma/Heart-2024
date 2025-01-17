library(tidyverse)
library(optparse)
option.list <- list(
    make_option("--summarySheets", type="character"),
    make_option("--samples", type="character"),
    make_option("--outputFile", type="character")
)
opt <- parse_args(OptionParser(option_list=option.list))
sheets= opt$summarySheets %>% strsplit(split=",") %>% unlist()
samples = opt$samples %>% strsplit(split=",") %>% unlist()
output = as.data.frame(matrix(nrow=length(sheets),ncol=5))
rownames(output) = samples
colnames(output)=c("Number_of_reads", "SequencingSaturation", "EstimatedNumberOfCells","MeanReadsPerCell","MeanGenePerCell")

getSample <- function(path){
	pathlist=unlist(strsplit(path, "/"))
	sample=pathlist[length(pathlist)-3]
	return(sample)
}
for (sheet in sheets){
	df = read.table(sheet, sep=",", header=F, row.names=1, stringsAsFactors=F)
	sample = getSample(sheet)
	print(sample)
	output[sample, "Number_of_reads"]=df["Number of Reads",1]
	output[sample, "SequencingSaturation"]=df["Sequencing Saturation",1]
	output[sample, "EstimatedNumberOfCells"]=df["Estimated Number of Cells",1]
	output[sample, "MeanReadsPerCell"]=df["Mean Reads per Cell",1]
	output[sample, "MeanGenePerCell"]=df["Mean GeneFull_Ex50pAS per Cell",1]
	output[sample, "ReadsWithValidBarcodes"]=df["Reads With Valid Barcodes",1]
	output[sample, "ReadsMapped2GeneFull_Ex50pAS"]=df["Reads Mapped to GeneFull_Ex50pAS: Unique GeneFull_Ex50pAS",1]
}

write.table(output, opt$outputFile, sep="\t", row.names=T, col.names=NA, quote=FALSE)

