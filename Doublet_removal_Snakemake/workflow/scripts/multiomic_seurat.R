suppressPackageStartupMessages(library(optparse))
option.list <- list(
    make_option("--inputSeuratObject", type="character"),
    make_option("--summary_sheet",type="character"),
    make_option("--outdir",type="character"),
    make_option("--sample", type="character"),
    make_option("--atacInfoSheet", type="character"),
    make_option("--GEXSaturation", type="character"),
    make_option("--min_cell", type="numeric")
    )

opt <- parse_args(OptionParser(option_list=option.list))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggrastr))
set.seed(42)

input_doublets <- read.table(opt$summary_sheet, header=TRUE, sep = '\t', stringsAsFactors=FALSE)
atacInfoSheet <- read.table(opt$atacInfoSheet, header=TRUE, sep = '\t', stringsAsFactors=FALSE)
s <- readRDS(opt$inputSeuratObject)
##############################################################################
cells_to_keep=input_doublets %>% filter(EffectiveCells=="True") %>% pull(cbc)
s <- subset(s, cells=cells_to_keep)
tokeep <- which(s@assays$RNA@counts %>% apply(1, function(x) (( as.numeric(x) > 0) %>% as.numeric %>% sum > opt$min_cell)))
s <- s[tokeep,]
metadata <- s@meta.data
metadata$orig.ident=opt$sample
s@meta.data=metadata
print(colnames(s@meta.data))
##############################################################################
metadata <- s@meta.data
df <- metadata %>% summarise(mean_UMI = mean(nUMI), median_UMI=median(nUMI),mean_gene_counts=mean(nGene), median_gene_counts=median(nGene), mean_ratio_mt=mean(ratio.mt), mean_nUMI_non_MT=mean(nUMI_non_MT), median_nUMI_non_MT=median(nUMI_non_MT), mean_ratio_mt=mean(ratio.mt), gex_cell_number=n()) %>% as.data.frame() %>% mutate(Sample=opt$sample) %>% relocate(Sample)
write.table(df,file.path(opt$outdir, paste0(opt$sample,"_gex_data_QC_post_doublet_removal_multiomic.tsv")), quote=FALSE, sep="\t", row.names=FALSE)
s$GEX_saturation=rep(as.numeric(opt$GEXSaturation, rep=length(colnames(s))))
saveRDS(s, file.path(opt$outdir, paste0(opt$sample, "_SeuratObject_doublet_removed_multiomic_UMAP.RDS")))
rm(s)
##############################################################################
#add UMAP
s <- readRDS(opt$inputSeuratObject)
cells_to_keep=input_doublets %>% filter(InATAC=="True") %>% pull(cbc)
atac_wo_doublets=input_doublets %>% filter(EffectiveCells=="True") %>% pull(cbc)
input_doublets <- input_doublets %>% left_join(atacInfoSheet, by=c("cbc"="GEX_Barcodes"))%>% mutate(log10_nFrags=log10(nFrags))
# Keep all the "cells" in both GEX and ATAC
s <- subset(s, cells=cells_to_keep)
s <- FindVariableFeatures(object = s)
all.genes <- rownames(s)
s <- ScaleData(s)
s <- RunPCA(s,features = VariableFeatures(object = s), verbose = FALSE)
pdf(file.path(opt$outdir,"elbow_plot.pdf"))
p <- ElbowPlot(s)
print(p)
dev.off()
if (length(colnames(s)) < 3000){
    UMAP_res=0.3
} else (
    UMAP_res = ((length(colnames(s))-3000) %/% 4000 + 1) * 0.3 + 0.4
)
s <- FindNeighbors(s, dims = 1:15)
s <- FindClusters(s, resolution = UMAP_res)
s <- RunUMAP(s, dims = 1:15)
s$Seurat_cluster <- Idents(s)
metadata <- s@meta.data
metadata <- metadata %>% left_join(input_doublets, by=c("cells"="cbc")) 
rownames(metadata) <- metadata$cells
s@meta.data=metadata

pdf(file.path(opt$outdir,paste0(opt$sample, "_nfrag_nUMI.pdf")), width=10, height=10)
p1 <- metadata %>% ggplot(aes(x=nUMI, y=nFrags)) + rasterize(geom_point(), dpi=300) + scale_x_log10() + scale_y_log10() + theme_classic() + ggtitle(paste0(opt$sample, "_nUMI vs. nFrags_all_cells"))
p2 <- ggplot(metadata, aes(x=nFrags/nUMI)) + geom_density(size=1) + theme_classic()+ ggtitle(paste0(opt$sample, "nFrags/nUMI_density_all_cells"))
singlet_metadata <- metadata %>% filter(cells %in% atac_wo_doublets)
p3 <- singlet_metadata %>% ggplot(aes(x=nUMI, y=nFrags)) + rasterize(geom_point(), dpi=300)+ scale_x_log10() + scale_y_log10() + theme_classic()+ ggtitle(paste0(opt$sample, "_nUMI vs. nFrags_doublet_removed"))
p4 <- ggplot(singlet_metadata, aes(x=nFrags/nUMI)) + geom_density(size=1) + theme_classic()+ ggtitle(paste0(opt$sample, "nFrags/nUMI_density_doublet_removed"))
#print(p1)
#print(p2)
#print(p3)
#print(p4)
print(ggarrange(p1, p2, p3, p4, ncol=2, nrow=2))
dev.off()

pdf(file.path(opt$outdir, paste0(opt$sample, "_doublet_removed_multiomic_UMAP.pdf")), width=15, height=10)
Idents(s)="Multiplets"
p1 <- DimPlot(s, split.by="ident", pt.size=0.5, cols = c('True' = 'red', 'False' = 'blue'),raster=TRUE) + ggtitle("All multiplets")
Idents(s)="scrublet_doublets"
p2 <- DimPlot(s,split.by="ident", pt.size=0.5, cols = c('True' = 'red', 'False' = 'blue'),raster=TRUE)  + ggtitle("Scrublet doublets")
Idents(s)="Amulet_multiplets"
p3 <- DimPlot(s,split.by="ident", pt.size=0.5, cols = c('True' = 'red', 'False' = 'blue'),raster=TRUE) + ggtitle("Amulet multiplets")

if ("demux_doublets" %in% colnames(input_doublets)){
    Idents(s)="demux_doublets"
    p4 <- DimPlot(s,split.by="ident", pt.size=0.5, cols = c('True' = 'red', 'False' = 'blue'),raster=TRUE) + ggtitle("Demux doublets")
    if ("ArchRDoublets" %in% colnames(input_doublets)){
        Idents(s)="ArchRDoublets"
        p5 <- DimPlot(s,split.by="ident", pt.size=0.5, cols = c('True' = 'red', 'False' = 'blue'),raster=TRUE) + ggtitle("ArchR doublets")
	plot <- ggarrange(p1, p2, p3, p4, p5, ncol=3, nrow=2)
	plot <- annotate_figure(plot, top=text_grob(opt$sample))
	print(plot)
    } else{
	plot <- ggarrange(p1, p2, p3, p4, ncol=2, nrow=2)
	plot <- annotate_figure(plot, top=text_grob(opt$sample))
        print(plot)
    } 
} else if ("ArchRDoublets" %in% colnames(input_doublets)) {
    p5 <- DimPlot(s,split.by="ident", pt.size=0.5, cols = c('True' = 'red', 'False' = 'blue'),raster=TRUE) + ggtitle("ArchR doublets")
    plot <- ggarrange(p1, p2, p3, p5, ncol=2, nrow=2)
    plot <- annotate_figure(plot, top=text_grob(opt$sample))
    print(plot)
} else{
    plot <- ggarrange(p1, p2, p3,  ncol=2, nrow=2)
    plot <- annotate_figure(plot, top=text_grob(opt$sample))
    print(plot)
}
dev.off()

atac_qc_df <- input_doublets %>% filter(EffectiveCells=="True") %>%
summarise(mean_ATAC_frags = round(mean(nFrags), digits=2), median_ATAC_frags=round(median(nFrags),digits=2),mean_ATAC_TSS_enrichment=round(mean(TSSEnrichment), digits=2), ATAC_cell_number=n())  %>% as.data.frame() %>% mutate(sample=opt$sample) %>% relocate(sample)
write.table(atac_qc_df, file.path(opt$outdir, paste0(opt$sample,"_ATAC_QC_post_doublet_removal_multiomic.tsv")), quote=FALSE, sep="\t", row.names=FALSE)

pdf(file.path(opt$outdir, paste0(opt$sample, "_doublet_removed_multiomic_ATAC.pdf")), width=10, height=10)
cells_to_keep=input_doublets %>% filter(EffectiveCells=="True") %>% pull(cbc)
s <- subset(s, cells=cells_to_keep)
Idents(s)="nFrags"
p6 <- FeaturePlot(s,features="log10_nFrags", order = TRUE, cols=c("light gray", "red"), pt.size=0.5, raster=TRUE)+ ggtitle(paste0(opt$sample, "log10(nUniqueFragments)"))
p7 <- ggplot(input_doublets, aes(x=log10_nFrags)) + geom_density()+ ggtitle(paste0(opt$sample, "log10(nUniqueFragments)"))
Idents(s)="TSSEnrichment"
p8 <- FeaturePlot(s,features="TSSEnrichment", order = TRUE, cols=c("light gray", "red"), pt.size=0.5,raster=TRUE)+ ggtitle(paste0(opt$sample, "TSSEnrichment"))
p9 <- ggplot(input_doublets, aes(x=TSSEnrichment)) + geom_density()+ ggtitle(paste0(opt$sample, "TSSEnrichment"))
plot <- ggarrange(p6, p7, p8, p9, ncol=2, nrow=2)
plot <- annotate_figure(plot, top=text_grob(opt$sample))
print(plot)
dev.off()
##############################################################################
