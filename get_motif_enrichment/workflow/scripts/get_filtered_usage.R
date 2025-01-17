library(tidyverse)
library(optparse)
option.list <- list(
    make_option("--usage_mtx", type="character"),
    make_option("--atac_cbc", type="character"),
    make_option("--outdir", type="character"),
    make_option("--cell_info", type="character"),
    make_option("--low_exp_gep", type="character"),
    make_option("--contamination_gep", type="character"),
    make_option("--sample", type="character")
)

opt <- parse_args(OptionParser(option_list=option.list))
usage_mtx <- read.table(opt$usage_mtx, header=T, row.names=1, sep="\t", stringsAsFactors=F)
cell_info <- read.table(opt$cell_info, header=T,sep="\t", stringsAsFactors=F)
atac_cbc <- read.table(opt$atac_cbc, header=F, sep="\t", comment.char = "", stringsAsFactors=F)
low_exp_gep <- read.table(opt$low_exp_gep, header=F,sep="\t", stringsAsFactors=F)
contamination_gep <- read.table(opt$contamination_gep, header=F,sep="\t", stringsAsFactors=F)
gep_remove=unique(c(low_exp_gep$V1,contamination_gep$V1))
#######################################################################################
atac_cbc <- atac_cbc %>% mutate(V2=str_replace(V1, paste0(opt$sample, "#"), "")) %>% rename("sample_atac_cbc"=V1)%>% rename("atac_cbc"=V2) %>% mutate(sample_atac_cbc=str_replace(sample_atac_cbc, paste0("#"), "."))
cell_info <- cell_info %>% select(sample_orig, GEX_Barcode, ATAC_barcode) %>% mutate(gex_cbc=paste0(sample_orig, "_", GEX_Barcode)) %>% mutate(atac_cbc=paste0(sample_orig, "_", ATAC_barcode)) %>% select(-c(GEX_Barcode, ATAC_barcode))
cell_info <- cell_info %>% inner_join(atac_cbc, by=c("atac_cbc"))
usage_mtx <- usage_mtx[cell_info$gex_cbc,]
colnames(usage_mtx) <- paste0(opt$sample,"_", colnames(usage_mtx))
usage_mtx <- usage_mtx %>% select(-any_of(gep_remove))
tmp_mtx <- usage_mtx %>% rownames_to_column(var="gex_cbc") %>% left_join(cell_info %>% select(c(gex_cbc,sample_atac_cbc))) 
dup=duplicated(tmp_mtx$sample_atac_cbc)

usage_mtx <- usage_mtx %>% rownames_to_column(var="gex_cbc") %>% left_join(cell_info %>% select(c(gex_cbc,sample_atac_cbc))) %>% select(-gex_cbc) %>% column_to_rownames(var="sample_atac_cbc")
write.table(usage_mtx, file.path(opt$outdir, paste0("filtered_usage_table_", opt$sample, ".tsv")), sep="\t", row.names=T, col.names=NA, quote=F)


