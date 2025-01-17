library(Seurat)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
suppressPackageStartupMessages(library(optparse))
set.seed(42)
option.list <- list(
	make_option("--SeuratPath", type="character", default=""),
	make_option("--scSHC_clustering", type="character", default=NA),
	make_option("--scSHC_train_clustering", type="character", default=NA),
	make_option("--subclustering_results", type="character", default=""),
	make_option("--scSHC_clustering_batch", type="character", default=NA),
	make_option("--scSHC_train_clustering_batch", type="character", default=NA),
	make_option("--outputPDF",type="character", default=""),
	make_option("--annotatedSeurat",type="character", default=""),
	make_option("--outputCellDistributionPDF",type="character", default="")
)
opt = parse_args(OptionParser(option_list=option.list))

s <- readRDS(opt$SeuratPath)

#################################################################################################
subclustering_table <- read.table(opt$subclustering_results, row.names=1, header=T, sep="\t", stringsAsFactors=F) %>% rownames_to_column("cbc")
annotation_table <- subclustering_table

if (!is.na(opt$scSHC_clustering)) {
	full_table <- read.table(opt$scSHC_clustering, row.names=1, header=T, sep="\t", stringsAsFactors=F) %>% rownames_to_column("cbc")
    colnames(full_table) <- c("cbc","scSHC_full_count_clusters")
	annotation_table <- annotation_table %>% left_join(full_table , by="cbc")

}

if (!is.na(opt$scSHC_train_clustering)) {
	train_table <- read.table(opt$scSHC_train_clustering, row.names=1, header=T, sep="\t", stringsAsFactors=F) %>% rownames_to_column("cbc")
	colnames(train_table) <- c("cbc","scSHC_train_count_clusters")
	annotation_table <- annotation_table %>% left_join(train_table , by="cbc")
}


batches=c()
if (!is.na(opt$scSHC_clustering_batch)){
	batch_annot=opt$scSHC_clustering_batch %>% strsplit(split=",") %>% unlist()
	for (annot in batch_annot){
		batch_annot_table=read.table(annot, row.names=1, header=T, sep="\t", stringsAsFactors=F) %>% rownames_to_column("cbc")
		batch_name=basename(annot)
		batch_name=str_replace(batch_name,"clustering_batch_", "")
		batch_name=str_replace(batch_name,".tsv","")
		colnames(batch_annot_table)=c("cbc",paste0("scSHC_full_count_clusters_batch_", batch_name))
		annotation_table=annotation_table%>% left_join(batch_annot_table, by="cbc")
		batches=c(batches,batch_name)
	}
}

if (!is.na(opt$scSHC_train_clustering_batch)){
	batch_annot=opt$scSHC_train_clustering_batch %>% strsplit(split=",") %>% unlist()
	for (annot in batch_annot){
		batch_annot_table=read.table(annot, row.names=1, header=T, sep="\t", stringsAsFactors=F) %>% rownames_to_column("cbc")
		batch_name=basename(annot)
		batch_name=str_replace(batch_name,"train_clustering_batch_", "")
		batch_name=str_replace(batch_name,".tsv","")
		colnames(batch_annot_table)=c("cbc",paste0("scSHC_train_count_clusters_batch_", batch_name))
		annotation_table=annotation_table%>% left_join(batch_annot_table, by="cbc")
	}
}

#################################################################################################
metadata <- s@meta.data
metadata <- metadata %>% select(-starts_with("RNA_snn")) %>% select(-starts_with("res_")) %>% select(-starts_with("scSHC_full_count_clusters")) %>% select(-starts_with("scSHC_train_count_clusters")) %>% left_join(annotation_table, by="cbc")
rownames(metadata) <- metadata$cbc
s@meta.data=metadata
saveRDS(s, opt$annotatedSeurat)

pdf(opt$outputPDF)
if (length(colnames(s)) > 100000){
	p1 <- DimPlot(s, group.by="Sex",raster=TRUE,shuffle=TRUE)
	p2 <- DimPlot(s, group.by="stage",raster=TRUE,shuffle=TRUE)
	print(p1)
	print(p2)
	for (annot in colnames(annotation_table)){
		if (annot != "cbc") {
			p <- DimPlot(s, group.by=annot,raster=TRUE) + ggtitle(annot)
			print(p)
		}
	}
	} else {
    p1 <- DimPlot(s, group.by="Sex",shuffle=TRUE)
    p2 <- DimPlot(s, group.by="stage",shuffle=TRUE)
	print(p1)
	print(p2)
	for (annot in colnames(annotation_table)){
		if (annot != "cbc"){
			p <- DimPlot(s, group.by=annot,shuffle=TRUE) + ggtitle(annot)
			print(p)
		}
	}
	}
dev.off()

#################################################################################################
plotCellDistribution <- function(metadata, reference, target){
	percentage_full<- metadata %>% select(!!sym(reference),!!sym(target)) %>% group_by(!!sym(reference),!!sym(target)) %>% summarize(n=n()) %>% ungroup() %>% as.data.frame() %>% group_by(!!sym(reference)) %>% mutate(cell_counts=sum(n)) %>% ungroup() %>% as.data.frame() %>% mutate(ratio=n/cell_counts) %>% mutate(ratio=round(ratio,digits=3))

	percentage_ratio_mtx <- percentage_full %>% select(!!sym(reference),!!sym(target),ratio) %>% mutate(!!sym(reference):=paste0("reference_",!!sym(reference))) %>% mutate(!!sym(target):=paste0("target_",!!sym(target))) %>% pivot_wider(names_from=!!sym(target), values_from=ratio, values_fill=0) %>% as.data.frame() %>% column_to_rownames(var=reference) %>% as.matrix %>% t()

	f2 = colorRamp2(breaks=c(0, 1), colors=c("#EEEEEE", "blue"), space = "RGB")
	p <- Heatmap(percentage_ratio_mtx,cluster_columns = TRUE, col=f2, show_column_dend = TRUE,show_row_dend = FALSE,cluster_rows=FALSE,use_raster=TRUE,rect_gp = gpar(col = "white", lwd = 2),column_title=paste0("reference_", reference, "_vs_target_", target),cell_fun = function(j, i, x, y, width, height, fill) {grid.text(sprintf("%.2f", percentage_ratio_mtx[i, j]), x, y, gp = gpar(fontsize = 10))})

	print(p)
}


