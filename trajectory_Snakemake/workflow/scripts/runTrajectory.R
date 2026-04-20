suppressPackageStartupMessages(library(optparse))
library(tidyverse)
library(slingshot)
library(viridis)
library(Seurat)
library(ggrastr)
library(ggplot2)
library(scales)
set.seed(42)
option.list <- list(
    make_option("--inputSeuratObject", type="character", default=""),
    make_option("--outdir", type="character", default=""),
    make_option("--sample", type="character", default=NA),
    make_option("--sampleInfo", type="character", default=NA)
)
opt <- parse_args(OptionParser(option_list=option.list))
object <- readRDS(opt$inputSeuratObject)

if (length(colnames(object) < 10000)) {
    approx_points_input = NULL
} else {
    approx_points_input = 200
}
sce<- as.SingleCellExperiment(object)
sce <- slingshot(sce, clusterLabels="RNA_snn_res.0.5", reducedDim="PCA")
saveRDS(sce, file.path(opt$outdir, paste0(opt$sample, "_slingshot.RDS")))
################################################################################################
colors <- colorRampPalette(viridis(10))(100)
pseudo.paths <- slingPseudotime(sce)
embedded <- embedCurves(sce, "UMAP")
embedded <- slingCurves(embedded)
pdf(file.path(opt$outdir, paste0(opt$sample,"_UMAP.trajectory.pdf")))
cluster_df <- data.frame(reducedDims(sce)$UMAP)
cluster_df$Cluster <- colData(sce)$RNA_snn_res.0.5

cluster_plot <- ggplot(cluster_df, aes(x=UMAP_1, y=UMAP_2, color=Cluster)) +
  rasterize(geom_point(size=0.8), dpi=300) +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(size=3))) +
  ggtitle(paste(opt$sample, "- Cluster Projection (res.0.5)"))
  print(cluster_plot)
for (path in colnames(pseudo.paths)){
  pseudo.path.subset=as.data.frame(pseudo.paths) %>% select(!!sym(path))
  umap_df=merge(data.frame(reducedDims(sce)$UMAP),pseudo.path.subset, by="row.names", all=TRUE)  %>% arrange(!is.na(!!sym(path)), !!sym(path))
  umap <- ggplot(umap_df, aes(x=UMAP_1, y=UMAP_2)) + rasterize(geom_point(mapping=aes(x=UMAP_1, y=UMAP_2, color=rescale(umap_df %>% pull(!!sym(path))))), dpi=300)+ scale_color_viridis(na.value="lightgray") + theme_bw()
  index=as.numeric(gsub("Lineage", "", path))
  path_object=embedded[[index]]
  embedded_path <- data.frame(path_object$s[path_object$ord,])
  umap <- umap + geom_path(data=embedded_path, mapping=aes(x=UMAP_1, y=UMAP_2), size=1.2)
  print(umap)
}
dev.off()
################################################################################################
print(head(pseudo.paths))
object <- readRDS(opt$inputSeuratObject)
metadata <- object@meta.data %>% rownames_to_column("cbc") %>% left_join(pseudo.paths %>% as.data.frame() %>% rownames_to_column("cbc"), by=c("cbc"))
if (!("Rounded.PCW" %in% colnames(metadata))){
    info <- read.table(opt$sampleInfo, header=T, sep="\t", stringsAsFactors=F) %>% mutate(cbc=paste0(sample_orig, "_", ATAC_barcode)) %>% select(cbc, Rounded.PCW)
    metadata <- metadata %>% left_join(info, by=c("cbc"))
}
write.table(pseudo.paths, file.path(opt$outdir, paste0(opt$sample, "_slingshot.tsv")), sep="\t", quote=F, row.names=T)

pdf(file.path(opt$outdir, paste0(opt$sample, "_boxplot.pseudotime.annot.pdf")))
for (path in colnames(pseudo.paths)){
    lineage=metadata %>% filter(!is.na(!!sym(path)))
    lineage_path=ggplot(metadata, aes(x=as.factor(Rounded.PCW), y=!!sym(path)))+geom_violin(width=1.4) +
    geom_boxplot(width=0.1, color="grey", alpha=0.2)  + ggtitle(path)
    print(lineage_path)
}
dev.off()
################################################################################################
sink(file = file.path(opt$outdir, "slingshot_summary.txt"))
print(SlingshotDataSet(sce))
sink(file = NULL)
