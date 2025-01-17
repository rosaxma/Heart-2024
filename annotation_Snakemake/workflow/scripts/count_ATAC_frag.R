library(tidyverse)
suppressPackageStartupMessages(library(optparse))
set.seed(42)
option.list <- list(
	make_option("--cluster_barcodes", type="character", default=""),
	make_option("--ATAC_QC_sheet ", type="character", default=""),
	make_option("--outDir",type="character", default="")
)
opt = parse_args(OptionParser(option_list=option.list))
#################################################################################################
cluster_barcodes_table <- read.table(opt$cluster_barcodes, header=T, sep="\t", stringsAsFactors=F)
ATAC_QC_sheet <- read.table(opt$ATAC_QC_sheet, header=T, sep="\t", stringsAsFactors=F) %>% mutate(GEX_Barcodes=paste0(Sample,"_",GEX_Barcodes))
cluster_barcodes_table <- cluster_barcodes_table %>% left_join(ATAC_QC_sheet %>% select(GEX_Barcodes,nFrags), by=c("cbc"="GEX_Barcodes")) %>% group_by(cluster_manual) %>% summarize(number_of_cells=n(), total_Frags=sum(nFrags)) %>% ungroup() %>% mutate(ABC_eligible=ifelse(total_Frags > 3000000, "T", "F"))
write.table(cluster_barcodes_table, file.path(opt$outDir, "ATAC_fragments_per_cluster.tsv"), sep="\t", row.names=F, quote=F)
