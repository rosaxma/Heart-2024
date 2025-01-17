suppressPackageStartupMessages(library(optparse))
option.list <- list(
    make_option("--outdir", type="character"), 
    make_option("--seurat", type="character"),
    make_option("--res", type="numeric")
)

opt <- parse_args(OptionParser(option_list=option.list))
library(Seurat)
library(SeuratObject)
library(utils)
suppressPackageStartupMessages(library(tidyverse))

s <- readRDS(opt$seurat)

cluster_col=paste0("res_", as.character(opt$res))
Idents(s) <- cluster_col
cluster_list <- s@meta.data %>% select(!!sym(cluster_col)) %>% distinct() %>% pull(!!sym(cluster_col))
Idents(s) <- cluster_col
library(Matrix)
counter=0
all_counter=0

all_df <- 0
for (cluster in cluster_list){
        all_counter=all_counter+1
        gene_matrix <- GetAssayData(object = subset(s, idents=cluster), slot = "counts")
        gene_sum <- Matrix::rowSums(gene_matrix)
        gex_df <- data.frame(genes=names(gene_sum), nUMI=gene_sum)
        norm_gex_df <- gex_df %>% mutate(nUMI=nUMI/sum(gene_sum)*1000000) %>% rename(!!sym(paste0(cluster,"_transcript_per_million")):="nUMI")
        if (all_counter ==1) {
                all_df <- norm_gex_df
        } else {
                all_df <- all_df %>% left_join(norm_gex_df, by="genes")
        }

}
write.table(all_df, file.path(opt$outdir,paste0("GEX_TPM.tsv")), sep="\t", quote=FALSE, row.names=FALSE)
