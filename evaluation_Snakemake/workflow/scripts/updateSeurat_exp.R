suppressPackageStartupMessages(library(optparse))
option.list <- list(
    make_option("--outdir", type="character"), 
    make_option("--seurat", type="character"),
    make_option("--effectiveSheet", type="character"),
    make_option("--genelist", type="character"),
    make_option("--annotations", type="character")
)

opt <- parse_args(OptionParser(option_list=option.list))

library(Seurat)
suppressPackageStartupMessages(library(tidyverse))

effective_sheet <- read.table(opt$effectiveSheet, header=T, sep="\t", stringsAsFactors=F) %>% mutate(cbc=paste0(orig.ident,"_",cbc))
annotations=opt$annotations %>% strsplit(split=",") %>% unlist()
counter=0
annot_df=""
print(annotations)
for (annot in annotations){
        counter=counter+1
        id=basename(annot)
        id=sub(".tsv", "", id)
        annotation <- read.table(annot, row.names=1, header=T, sep=",", stringsAsFactors=F) %>% rownames_to_column("cbc") %>% select(cbc, majority_voting) %>% rename(!!sym(id):=majority_voting)
        print(head(annotation))
        if (counter==1){
                annot_df=annotation
        } else {
                annot_df <- annot_df %>% left_join(annotation, by="cbc")
        }
}
print(head(effective_sheet))
print(head(annot_df))
info_sheet <- effective_sheet %>% inner_join(annot_df, by="cbc")

s <- readRDS(opt$seurat)
s <- subset(s, cells=effective_sheet$cbc)

metadata <- s@meta.data
metadata$cbc <- rownames(metadata)
print(head(metadata))
print(head(info_sheet))
for (col in colnames(info_sheet)) {
	if ((col %in% colnames(metadata))& !(col %in% c("cbc", "orig.ident"))){
		info_sheet <- info_sheet %>% select(-sym(col))
	}
}
metadata <- metadata %>% left_join(info_sheet , by=c("cbc", "orig.ident"))
rownames(metadata) <- metadata$cbc
s@meta.data <- metadata
write.table(metadata, file.path(opt$outdir, "metadata.tsv"), row.names=T, quote=F, col.names = NA, sep="\t")

file.name=basename(opt$seurat)
saveRDS(s,file.path(opt$outdir, gsub(".RDS", ".source.RDS",file.name)))

Idents(s) <- "source"
library(Matrix)
all_counter=0
all_df <- 0
gene_list <- read.table(opt$genelist, sep="\t", header=FALSE, stringsAsFactors=FALSE)
print(head(metadata))
source_list <- unique(metadata$source)
print(source_list)
for (source in source_list){
    all_counter=all_counter+1
    gene_matrix <- GetAssayData(object = subset(s, idents=source), slot = "counts")
    print(gene_matrix[1:10,1:1])
    print(rownames(gene_matrix)[1:10])
    print(gene_list[1:10,1])
    existing_genes=intersect(gene_list[,1], rownames(gene_matrix))
    if (is.null(dim(gene_matrix[existing_genes,]))){
        next
        } else {
        gene_sum <- Matrix::rowSums(gene_matrix[existing_genes,])
         gex_df <- data.frame(genes=names(gene_sum), nUMI=gene_sum)
         print(source)
         norm_gex_df <- gex_df %>% mutate(nUMI=nUMI/sum(gene_sum)*1000000) %>% rename(!!sym(paste0(source,"_transcript_per_million")):="nUMI")
    if (all_counter ==1) {
            all_df <- norm_gex_df
    } else {
            all_df <- all_df %>% left_join(norm_gex_df, by="genes")
     }
   }
}

gene_df <- all_df %>% filter(genes %in% gene_list[,1])
write.table(gene_df, file.path(opt$outdir, "GEX_TPM_sex.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
library(ComplexHeatmap)
rownames(gene_df) <- gene_df$genes
gene_df <- gene_df[order(row.names(gene_df)), ]
mat <- gene_df %>% select(-genes) %>% as.matrix()

colnames(mat) <- sub("_transcript_per_million", "", colnames(mat))
mat <- t(scale(t(mat)))
Seurat::BlueAndRed()
quantile_table=quantile(mat, c(0.1, 0.95), na.rm=TRUE)
quantile_low=as.numeric(quantile_table["10%"])
quantile_high=as.numeric(quantile_table["95%"])
col_fun = circlize::colorRamp2(c(quantile_low, 0, quantile_high), c("#313695", "white", "#A50026"))
mat=na.omit(mat)
p <- Heatmap(mat, name = "Transcript_per_million",
        cluster_columns = TRUE,
        show_column_dend = TRUE,
        show_row_dend = TRUE,
        cluster_rows=TRUE,
        col = col_fun,
        row_names_gp = gpar(fontsize = 4),
        column_title_gp = gpar(fontsize = 2),
        show_column_names = TRUE,
        use_raster = TRUE,
        column_title_side = "top", raster_quality = 4)
print(dim(mat))
pdf(file.path(opt$outdir, "GEX_TPM_sex.pdf"))
print(p)
dev.off()

