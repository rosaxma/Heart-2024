suppressPackageStartupMessages(library(optparse))
option.list <- list(
    make_option("--outdir", type="character", default="", help="Output directory"), 
    make_option("--seurat", type="character", default="", help="Path to the Seurat Object"),
    make_option("--sampleIDColName", type="character", default="")
)
opt <- parse_args(OptionParser(option_list=option.list))

library(SeuratObject)
library(Seurat)
suppressPackageStartupMessages(library(tidyverse))

s <- readRDS(opt$seurat)
print(length(colnames(s)))
set.seed(seed=42)
cellIDs <- sample(colnames(s))

metadata <- s@meta.data
sample_ids <- metadata %>% select(!!sym(opt$sampleIDColName)) %>% distinct() %>% pull(!!sym(opt$sampleIDColName))
sample_id_list <-  metadata %>% select(!!sym(opt$sampleIDColName)) %>% pull(!!sym(opt$sampleIDColName))
print(table(sample_id_list))
print(length(sample_id_list))
rep_id_list <- c()
counter=0
for (sample in sample(sample_ids)){
    counter = counter + 1
    print(sample)
    count <- sum(sample_id_list == sample)
    rep_id_list <- c(rep_id_list, rep(paste0("control_rep_", as.character(counter)), count))
}
print(length(rep_id_list))
print(rep_id_list[1:10])
newID_df <- data.frame(cellIDs, rep_id_list)
colnames(newID_df) <- c("cellbarcodes", "Rep_ID")
print(head(newID_df))
metadata <- s@meta.data
metadata$cellbarcodes <- rownames(metadata)
print(head(metadata))
metadata <- metadata %>% left_join(newID_df, by="cellbarcodes")
rownames(metadata) <- metadata$cellbarcodes
metadata <- metadata %>% select(-cellbarcodes)
print(table(metadata$Rep_ID))
s@meta.data <- metadata
print(str(s))
saveRDS(s, file.path(opt$outdir, "control.RDS"))

