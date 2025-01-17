library(Seurat)
library(tidyverse)
#library(ComplexHeatmap)
#library(circlize)
suppressPackageStartupMessages(library(optparse))
set.seed(42)
option.list <- list(
	make_option("--SeuratPath_full", type="character", default=""),
	make_option("--SeuratPath_test", type="character", default=""),
	make_option("--experiment_info_table", type="character", default=""),
	make_option("--Sample", type="character", default=""),
	make_option("--res", type="character", default=""),
	make_option("--outdir",type="character", default=""),
	make_option("--key_feature_table",type="character", default=""),
	make_option("--category",type="character", default=""),
	make_option("--clusters2remove", type="character", default="NA"),
	make_option("--GEX_TPM_table", type="character", default="NA"),
	make_option("--exclusive", type="logical", default=NA)
)
opt = parse_args(OptionParser(option_list=option.list))
cluster_col=paste0("res_", opt$res)
print(cluster_col)
cluster_col=gsub("\\.0","",cluster_col)
#################################################################################################
s_test <- readRDS(opt$SeuratPath_test)
s_full <- readRDS(opt$SeuratPath_full)


metadata <- s_full@meta.data
og_full_metadata <- metadata %>% group_by(orig.ident) %>% mutate(nUMI_zscore=scale(nUMI))

print(colnames(s_test@meta.data))
print(colnames(s_full@meta.data))
print(cluster_col)

metadata_test <- s_test@meta.data
metadata_full <- s_full@meta.data
exp_info_table <- read.table(opt$experiment_info_table, header=T, sep="\t",stringsAsFactors=F)
for (col in colnames(exp_info_table)){
	if (col != "orig.ident") {
		if (col %in% colnames(metadata_test)){
			metadata_test <- metadata_test %>% select(-!!sym(col))
		}
		if (col %in% colnames(metadata_full)){
			metadata_full <- metadata_full %>% select(-!!sym(col))
		}
	}

}
metadata_test <- metadata_test %>% left_join(exp_info_table, by="orig.ident")
rownames(metadata_test)=metadata_test$cbc
s_test@meta.data=metadata_test

metadata_full <- metadata_full %>% left_join(exp_info_table, by="orig.ident")
rownames(metadata_full)=metadata_full$cbc
s_full@meta.data=metadata_full


if (opt$clusters2remove!="NA"){
	clusters2remove <- read.table(opt$clusters2remove, header=F, sep="\t", stringsAsFactors=F) %>% filter(grepl(opt$Sample,V1)) %>% mutate(V1=str_replace(V1,paste0(opt$Sample, "_"), "")) %>% mutate(V1=as.character(V1))
	print(clusters2remove)
	Idents(s_test) <- cluster_col
	cluster_ids=unique(metadata_test %>% mutate(!!sym(cluster_col):=as.character(!!sym(cluster_col))) %>% pull(!!sym(cluster_col)))
	print(cluster_ids)
		if (length(intersect(cluster_ids, clusters2remove$V1)) > 0 ){
			s_test=subset(x = s_test, idents=clusters2remove$V1, invert = TRUE)
			Idents(s_full) <- cluster_col
			s_full=subset(x = s_full, idents=clusters2remove$V1, invert = TRUE)
			s_full_filename=basename(opt$SeuratPath_full)
			saveRDS(s_full, file.path(opt$outdir, gsub(".RDS", ".doublet.cluster.removed.RDS", s_full_filename)))
		} else {
		Idents(s_test) <- cluster_col
		Idents(s_full) <- cluster_col
	}
} else {
	Idents(s_test) <- cluster_col
	Idents(s_full) <- cluster_col
}
my_levels <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)
Idents(s_test) <- factor(Idents(s_test), levels=my_levels)
Idents(s_full) <- factor(Idents(s_full), levels=my_levels)
markers <- FindAllMarkers(s_test)
top <- markers %>% group_by(cluster) %>% filter(p_val_adj < 0.05,avg_log2FC>0) %>% arrange(p_val_adj, -avg_log2FC)
bottom <- markers %>% group_by(cluster) %>% filter(p_val_adj < 0.05,avg_log2FC<0) %>% arrange(p_val_adj, avg_log2FC)
combined_deg=rbind(top,bottom) %>% arrange(cluster)
write.table(combined_deg,file.path(opt$outdir, "top_bottom_DEGs.tsv"), row.names=F, sep="\t",quote=F)
top_10 <- markers %>% group_by(cluster) %>% filter(p_val_adj < 0.05,avg_log2FC>0) %>% slice_head(n = 15) %>% ungroup()
s_test <- FindVariableFeatures(object = s_test)
s_test <- ScaleData(s_test)
p <- DoHeatmap(s_test, features = top_10$gene) + NoLegend()
nCluster=length(unique(markers$cluster))
pdf(file.path(opt$outdir, "top_10_DEG_heatmap.pdf"),height=nCluster*3, width=nCluster*2 )
print(p)
dev.off()
#################################################################################################
## Dot plot
key_feature_table <- read.table(opt$key_feature_table, header=T, sep="\t", stringsAsFactors=F) %>% column_to_rownames(var ="Category")
print(head(key_feature_table))
feature_list_string=key_feature_table[opt$category,"Genes"]
print(feature_list_string)
feature_list=feature_list_string%>% strsplit(split=",") %>% unlist()
feature_vector=c()
for (feature in feature_list){
	feature_vector=c(feature_vector,trimws(feature))
}
feature_vector <- unique(feature_vector)
print(feature_vector)
pdf(file.path(opt$outdir, "Key_gene_dot_plot.pdf"), width=unique(markers$cluster)*2, height=length(feature_vector)*0.5)
p <- DotPlot(s_full, features = feature_vector)+coord_flip()
print(p)
dev.off()
#################################################################################################
distribution <- function(metadata, reference, target, outdir){
    count_table <- metadata %>% select(!!sym(reference), !!sym(target)) %>% group_by(!!sym(reference), !!sym(target)) %>% summarize(n=n()) %>% ungroup() %>% distinct() %>% pivot_wider(names_from=!!sym(target), values_from=n) %>% as.data.frame()
    write.table(count_table, file.path(outdir, paste0(reference, "_vs_",target,".tsv")), sep="\t", quote=F, row.names=F)
}
average <- function(metadata, reference, target, outdir){
    count_table <- metadata %>% select(c(!!sym(reference), !!sym(target), orig.ident)) %>% group_by(!!sym(reference), orig.ident) %>% summarize(!!sym(paste0("mean_", target)) := mean(!!sym(target))) %>% ungroup() %>% distinct() %>% as.data.frame()
    write.table(count_table, file.path(outdir, paste0("average_", reference, "_vs_",target, ".tsv")), sep="\t", quote=F, row.names=F)
}
#################################################################################################
metadata <- s_full@meta.data
target_features=c("Rounded.PCW", "Sex", "orig.ident", "Phase", "source", "experiment")
for (target in target_features) {
	distribution(metadata, cluster_col,target,opt$outdir)
}
target_features=c("ratio.mt","ratio.ribo","nUMI_zscore")
for (target in target_features) {
	average(og_full_metadata, cluster_col, target, opt$outdir)
}
#################################################################################################
z_score_table <- read.table(file.path(opt$outdir, paste0("average_", cluster_col, "_vs_",target, ".tsv")), header=T, sep="\t") %>% mutate(!!sym(cluster_col):=as.factor(!!sym(cluster_col)))
summary_z_score_table <- z_score_table %>% group_by(!!sym(cluster_col)) %>% mutate(total_sample=n()) %>% filter(mean_nUMI_zscore > 1) %>% mutate(sample_large_zscore=n()) %>% mutate(ratio=sample_large_zscore/total_sample) %>% ungroup() %>% select(!!sym(cluster_col), total_sample, sample_large_zscore, ratio) %>% distinct()
write.table(summary_z_score_table, file.path(opt$outdir, paste0("average_", cluster_col, "_vs_zscore_summary.tsv")), sep="\t", quote=F, row.names=F)

p <- ggplot(z_score_table, aes(x=mean_nUMI_zscore, color=!!sym(cluster_col))) + geom_density()+ theme_minimal()
pdf(file.path(opt$outdir, paste0("average_", cluster_col, "_vs_zscore_distribution.pdf")), width=15, height=5)
print(p)
dev.off()
#################################################################################################
# plot UMAP
target_features=c("Rounded.PCW", "Sex", "orig.ident", "Phase", "source", "experiment")
pdf(file.path(opt$outdir, "UMAP_colored_by_sample_feature.pdf"))
raster_bool=length(colnames(s_full)) > 100000
p <- DimPlot(s_full, group.by=cluster_col, shuffle=TRUE, label=TRUE, raster=raster_bool)
print(p)
for (target_feature in target_features){
	p <- DimPlot(s_full, group.by=target_feature, shuffle=TRUE, raster=raster_bool)
	print(p)		
}
p <- FeaturePlot(s_full, features="ratio.mt",order=TRUE, raster=raster_bool)
print(p)
p <- FeaturePlot(s_full, features="ratio.ribo",order=TRUE, raster=raster_bool)
print(p)
dev.off()
#################################################################################################
p <- VlnPlot(s_full, features=feature_vector, split.plot = TRUE, split.by="Rounded.PCW", group.by=cluster_col, pt.size=0, ncol=1)
pdf(file.path(opt$outdir, "violin_plot.pdf"), width=unique(markers$cluster)*2, height=length(feature_vector)*6+3)
print(p)
dev.off()
pdf(file.path(opt$outdir, "UMAP_colored_by_source.pdf"), width=8, height=5)
p <- DimPlot(s_full, group.by="source", shuffle=TRUE, raster=raster_bool)
print(p)
dev.off()
#################################################################################################
if (opt$GEX_TPM_table != "NA"){
	all_df <- read.table(opt$GEX_TPM_table, header=T, sep="\t", stringsAsFactors=F)
}
barcodes_w_annotation<- metadata %>% select(cbc,!!sym(cluster_col))
if (opt$clusters2remove != "NA"){
	clusters2remove <- read.table(opt$clusters2remove, header=F, sep="\t", stringsAsFactors=F)%>% filter(grepl(opt$Sample,V1)) %>% mutate(V2=str_replace(V1,paste0(opt$Sample, "_"), ""))
	Idents(s_test) <- cluster_col
        cluster_ids=unique(metadata_test %>% pull(cluster_col))
        if (length(intersect(cluster_ids, clusters2remove)) > 0 ){
		barcodes_w_annotation <- barcodes_w_annotation %>% filter(!(!!sym(cluster_col) %in% clusters2remove$V2))
		for (cluster in clusters2remove$V1){
        		col=paste0(cluster,"_transcript_per_million")
				if (opt$GEX_TPM_table != "NA"){
        			all_df <- all_df %>% select(-!!sym(col))
				}
		}
	}
}
barcodes_annotation <- barcodes_w_annotation %>% mutate(coarse_cluster=opt$Sample) %>% mutate(!!sym(cluster_col):=paste0(opt$Sample, "_", !!sym(cluster_col)))
write.table(barcodes_annotation, file.path(opt$outdir, "barcode_w_annotations.tsv"), sep="\t", quote=F, row.names=F)



if (opt$GEX_TPM_table != "NA"){
	if (!opt$exclusive){
		all_df <- all_df %>% select(starts_with(opt$Sample),"genes")
	}
	library(ComplexHeatmap)
	gene_df <- all_df %>% filter(genes %in% feature_list)
	rownames(gene_df) <- gene_df$genes
	gene_df <- gene_df[order(row.names(gene_df)), ]
	mat <- gene_df %>% select(-genes) %>% drop_na() %>% as.matrix()
	colnames(mat) <- sub("_transcript_per_million", "", colnames(mat))
	colnames(mat) <- sub("HEALTHY_EARLY_","",colnames(mat))
	mat <- t(scale(t(mat)))
	mat <- mat[!rowSums(is.na(mat)),]
	Seurat::BlueAndRed()
	quantile_table=quantile(mat, c(0.1, 0.95), na.rm=TRUE)
	quantile_low=as.numeric(quantile_table["10%"])
	quantile_high=as.numeric(quantile_table["95%"])
	col_fun = circlize::colorRamp2(c(quantile_low, 0, quantile_high), c("#313695", "white", "#A50026"))
	cluster_col="cluster_manual"
	p <- Heatmap(mat, name = "Transcript_per_million",
			cluster_columns = TRUE,
			show_column_dend = TRUE,
			show_row_dend = TRUE,
			cluster_rows=TRUE,
			col = col_fun,
			row_names_gp = gpar(fontsize = 20),
			column_names_gp = gpar(fontsize=15),
			column_title_gp = gpar(fontsize = 2),
			show_column_names = TRUE,
			use_raster = TRUE,
			column_title_side = "top",
			raster_quality = 4)
	print(dim(mat))
	pdf(file.path(opt$outdir,paste0("GEX_heatmap_",cluster_col,".pdf")), height=nrow(mat)*0.8, width=ncol(mat)*0.5+2)
	print(p)
	dev.off()
	write.table(gene_df, file.path(opt$outdir, paste0("GEX_heatmap_tpm_", cluster_col, ".tsv")), row.names=F, sep="\t",quote=F)
}
