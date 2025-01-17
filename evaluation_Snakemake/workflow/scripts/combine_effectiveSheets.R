suppressPackageStartupMessages(library(optparse))

option.list <- list(
    make_option("--outdir", type="character"), 
    make_option("--allEffectiveSheets", type="character")
    )

opt <- parse_args(OptionParser(option_list=option.list))
library(SeuratObject)
library(Seurat)
library(Matrix)
suppressPackageStartupMessages(library(tidyverse))
set.seed(42)

sheets=opt$allEffectiveSheets %>% strsplit(split=",") %>% unlist()

sheet_df_list <- vector(mode = "list", length = length(sheets))

for (i in 1:length(sheet_df_list)){
	df=read.table(sheets[i], sep="\t", header=T, stringsAsFactors=F) %>% select(cbc,Multiplets,EffectiveCells,orig.ident,source)
	sheet_df_list[[i]]=df
}
combined_df <- do.call("rbind", sheet_df_list)
write.table(combined_df, file.path(opt$outdir,"combined_effective_sheets.tsv"),sep="\t", quote=F, row.names=F)


