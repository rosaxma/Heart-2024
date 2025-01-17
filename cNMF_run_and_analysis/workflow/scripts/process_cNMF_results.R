library(SeuratObject)
library(Seurat)
library(Matrix)
library(optparse)
library(dplyr)
library(tidyr)
library(ggpubr)
option.list <- list(
    make_option("--SeuratObject", type="character"),
    make_option("--input_dir", type="character"),
    make_option("--name", type="character"),
    make_option("--numK", type="character"),
    make_option("--density", type="character"),
    make_option("--outdir", type="character")
)
opt <- parse_args(OptionParser(option_list=option.list))

s <- readRDS(opt$SeuratObject)

k=opt$numK
density=opt$density
beta_filename=paste0(opt$name, ".gene_spectra_score.k_", k, ".dt_", density, ".txt")
gene_beta <- read.table(file.path(opt$input_dir,beta_filename), header=T, row.names=1, stringsAsFactors=F, sep="\t") %>% t()
colnames(gene_beta) <- paste0("GEP_", c(1:as.numeric(k)), "_beta_score")


usage_filename=paste0(opt$name, ".usages.k_", k, ".dt_", density, ".consensus.txt")
gep_usage <- read.table(file.path(opt$input_dir,usage_filename), header=T, row.names=1, stringsAsFactors=F, sep="\t")%>% apply(1, function(x) x*100/sum(x)) %>% t()
colnames(gep_usage) <- paste0("GEP_", gsub("X", "", colnames(gep_usage)))
s <- AddMetaData(s, metadata=gep_usage, col.name=paste0("k_", k, "_", colnames(gep_usage)))

write.table(gene_beta, file.path(opt$outdir,"tables", paste0("beta_score_K_", k, "_df_", density, ".tsv")), sep="\t", row.names=T, quote=F, col.names=NA)
write.table(gep_usage, file.path(opt$outdir,"tables", paste0("gep_usages_K_", k, "_df_", density, ".tsv")), sep="\t", row.names=T, quote=F, col.names=NA)

#plot UMAP
plot_list <- vector(mode="list", length=length(colnames(gep_usage)))
ncol=4
if ((as.numeric(k) %% ncol) != 0){
    nrow=as.numeric(k)/ncol + 1
 } else {
    nrow=as.numeric(k)/ncol
}
feature_list= col.name=paste0("k_", k, "_", colnames(gep_usage))
for (n in c(1: length(colnames(gep_usage)))){
    gep <- feature_list[n]
    print(gep)
    p <- FeaturePlot(s, raster=T, reduction="umap", features=gep, order=TRUE, pt.size=2)
    plot_list[[n]] <- p 
}
pdf(file.path(opt$outdir,"figures", "GEP.expression.UMAP.pdf"), height=nrow*2, width=ncol*3)
print(ggarrange(plotlist=plot_list, ncol=ncol, nrow=nrow))
dev.off()


