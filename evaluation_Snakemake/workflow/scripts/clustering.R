suppressPackageStartupMessages(library(optparse))

option.list <- list(
    make_option("--outdir", type="character"), 
    make_option("--inputSeuratObject", type="character"),
    make_option("--res_range", type="character", default="0,1"),
    make_option("--increment", type="numeric", default="0.1")
)

opt <- parse_args(OptionParser(option_list=option.list))
library(SeuratObject)
library(Seurat)
library(ggplot2)
suppressPackageStartupMessages(library(tidyverse))
set.seed(42)

#############################################################################
s <- readRDS(opt$inputSeuratObject)

meta <- s@meta.data
meta <- meta %>% select(-starts_with("RNA_snn")) %>% select(-starts_with("res_"))
s@meta.data <- meta
#####################################
s <- FindVariableFeatures(object = s)
s <- ScaleData(s)
s <- RunPCA(s,features = VariableFeatures(object = s), verbose = FALSE)
pdf(file.path(opt$outdir,"elbow_plot.pdf"))
p <- ElbowPlot(s)
print(p)
dev.off()
s <- FindNeighbors(s, dims = 1:15)
s=DietSeurat(s,scale.data = FALSE, dimreducs="pca",graphs=c("RNA_nn","RNA_snn"))
#############################################################################
print(opt$res_range)
start=as.numeric(unlist(strsplit(opt$res_range, split=","))[1])
print(start)
end=as.numeric(unlist(strsplit(opt$res_range, split=","))[2])
print(end)
res_list <- seq(start, end, by=as.numeric(opt$increment))
print(res_list)
cluster_count <- c()
for (res in res_list){
    s <- FindClusters(s, resolution=as.numeric(res))
    print(str(s))
    metadata <- s@meta.data
    rownames <- rownames(metadata)
    colname <- paste0("res_", as.character(res))
    count <- length(unique(metadata$seurat_clusters))
    metadata <- metadata %>% rename(!!colname :="seurat_clusters")
    rownames(metadata) <- rownames
    s@meta.data <- metadata
    cluster_count <- c(cluster_count,count)
}
s <- RunUMAP(s, dims=1:15)
saveRDS(s,gsub(".RDS", ".cluster.pick.res.RDS", opt$inputSeuratObject))
print(head(s@meta.data))
#############################################################################
res_count=length(res_list)

pdf(file.path(opt$outdir, "clustree.pdf"), height=20, width=30)
library(clustree)
print(str(s))
if (res_count>1) {
	clustree(s, prefix = "res_", layout = "sugiyama")
}
dev.off()

df <- data.frame(res_list,cluster_count)
write.table(df, file.path(opt$outdir,"cluster_count_res.tsv"), quote=FALSE, sep = "\t", row.names = FALSE)
