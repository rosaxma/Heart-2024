library(tidyverse)
library(optparse)
option.list <- list(
    make_option("--barcodeSheets", type="character"),
    make_option("--markdupSheets", type="character"),
    make_option("--samples", type="character"),
    make_option("--outputFile", type="character")
)


opt <- parse_args(OptionParser(option_list=option.list))
print(opt)
barcodesheets= opt$barcodeSheets %>% strsplit(split=",") %>% unlist()
markdupsheets= opt$markdupSheets %>% strsplit(split=",") %>% unlist()
samples = opt$samples %>% strsplit(split=",") %>% unlist()
output = as.data.frame(matrix(nrow=length(samples),ncol=5))
rownames(output) = samples
colnames(output)=c("NUMBER_OF_READS", "NUMBER_OF_READS_with_valid_barcodes","READ_PAIRS_EXAMINED", "ESTIMATED_LIBRARY_SIZE","PERCENT_DUPLICATION")

getMarkDupSample <- function(path){
        pathlist=unlist(strsplit(path, "/"))
        sample=pathlist[length(pathlist)-1]
        return(sample)
}

getBarcodeSample <- function(path){
	pathlist=unlist(strsplit(path, "/"))
	sample=pathlist[length(pathlist)-2]
	return(sample)
}
print(barcodesheets)
for (sheet in barcodesheets){
        sample = getBarcodeSample(sheet)
	str_total_reads = read.table(sheet, sep=' ', header=F, nrows=1, stringsAsFactors=F)	
	output[sample, "NUMBER_OF_READS"]=as.numeric(unlist(strsplit(str_total_reads[1,"V1"],"/"))[2])
	output[sample, "NUMBER_OF_READS_with_valid_barcodes"]=as.numeric(unlist(strsplit(str_total_reads[1,"V1"],"/"))[1])
}

print(markdupsheets)
for (sheet in markdupsheets){
	sample = getMarkDupSample(sheet)
	temp =readLines(sheet)
	start=grep('ESTIMATED_LIBRARY_SIZE', temp)
	end=start+1
	markdup_file=read.table(text=paste0(trimws(temp[start:end]), collapse="\n"),sep="\t", header=T)
	output[sample, "READ_PAIRS_EXAMINED"]=as.numeric(markdup_file[1,"READ_PAIRS_EXAMINED"])
	output[sample, "ESTIMATED_LIBRARY_SIZE"]=as.numeric(markdup_file[1,"ESTIMATED_LIBRARY_SIZE"])
	output[sample, "PERCENT_DUPLICATION"]=as.numeric(markdup_file[1,"PERCENT_DUPLICATION"])
}
#print(output)
#print(opt$outputFile)
write.table(output, opt$outputFile, sep="\t", row.names=T, col.names=NA, quote=FALSE)
