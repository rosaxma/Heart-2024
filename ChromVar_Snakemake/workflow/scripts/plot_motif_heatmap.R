suppressPackageStartupMessages(library(optparse))
option.list <- list(
    make_option("--GenomeAssembly", type="character"),
    make_option("--ArchR", type="character"),
    make_option("--Outdir", type="character"),
    make_option("--SampleName", type="character"),
    make_option("--Threads", type="numeric"),
    make_option("--TF_regulators", type="character"),
    make_option("--clustering_info", type="character"),
    make_option("--GeneExpTPM", type="character")
)

opt <- parse_args(OptionParser(option_list=option.list))
suppressPackageStartupMessages(library(ArchR))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(Seurat))
set.seed(42)

addArchRGenome(opt$GenomeAssembly)
addArchRThreads(threads=opt$Threads)

#######################################################################
proj <- loadArchRProject(opt$ArchR) 
#####################################################################
motif_mtx <- getMatrixFromProject(ArchRProj=proj, useMatrix="MotifMatrix")
motif_z_mtx <- assays(motif_mtx)$z
tf_regulators <- read.table(opt$TF_regulators, header=T, sep="\t", stringsAsFactors=F)
tf <- tf_regulators %>% filter(TFRegulator=="YES") %>% pull(MotifMatrix_name)
motif_tf_z_mtx=motif_z_mtx[tf, ]
info <- read.table(opt$clustering_info, header=T, sep="\t", stringsAsFactors=F)
#####################################################################
info <- info %>% mutate(atac_cbc=paste0(opt$SampleName, "#", sample_orig, "_", ATAC_barcode))
info <- info %>% filter(atac_cbc %in% colnames(motif_tf_z_mtx)) %>% select(atac_cbc, celltype)
#####################################################################
motif_tf_z_df <- as.data.frame(as.matrix(motif_tf_z_mtx)) %>% rownames_to_column(var="motif") %>% pivot_longer(cols=starts_with(opt$SampleName), names_to="atac_cbc", values_to="chromvar_z_score")
motif_tf_z_df <- motif_tf_z_df %>% left_join(info) %>% select(-atac_cbc)
motif_tf_z_mtx_2 <- motif_tf_z_df %>% group_by(motif, celltype) %>% summarize(avg_chromvar_z_score=mean(chromvar_z_score))
motif_tf_z_mtx_2 <- motif_tf_z_mtx_2 %>% pivot_wider(names_from=celltype, values_from=avg_chromvar_z_score) %>% column_to_rownames(var="motif")
#####################################################################
suppressPackageStartupMessages(library(ComplexHeatmap))
mat <- t(as.matrix(motif_tf_z_mtx_2))
colnames(mat) <- sub("Motif.", "", colnames(mat))
#####################################################################
p <- Heatmap(mat, name = "ChromVar_z_score",
        cluster_columns = TRUE,
        show_column_dend = TRUE,
        show_row_dend = TRUE,
        cluster_rows=TRUE,
        row_names_gp = gpar(fontsize = 4),
        column_title_gp = gpar(fontsize = 2),
        show_column_names = TRUE,
        use_raster = TRUE,
        column_title_side = "top",
        raster_quality = 4)
pdf(file.path(opt$Outdir, "motif_heatmap.pdf"))
print(p)
dev.off()
#####################################################################
tf_gene <- tf_regulators %>% filter(TFRegulator=="YES") %>% pull(GeneScoreMatrix_name)
print(tf_gene[1:10])
all_df <- read.table(opt$GeneExpTPM, header=T, sep="\t", stringsAsFactors=F)
ct2keep <- unique(info$celltype)
tf_gene_df <- all_df %>% filter(genes %in% tf_gene) %>% column_to_rownames(var="genes")
colnames(tf_gene_df) <- gsub("_transcript_per_million", "", colnames(tf_gene_df))
print(dim(tf_gene_df))
print(colnames(tf_gene_df))
print(ct2keep)
tf_gene_df=tf_gene_df[ct2keep]
mat <- as.matrix(tf_gene_df)
colnames(mat) <- sub("GEX.", "", colnames(mat))
mat <- t(scale(t(mat)))
mat <- mat[!rowSums(is.na(mat)),]
mat <- t(mat)
#Seurat::BlueAndRed()
#quantile_table=quantile(mat, c(0.1, 0.95), na.rm=TRUE)
#quantile_low=as.numeric(quantile_table["10%"])
#quantile_high=as.numeric(quantile_table["95%"])
#col_fun = circlize::colorRamp2(c(quantile_low, 0, quantile_high), c("#313695", "white", "#A50026"))
print(head(mat))
p <- Heatmap(mat, name = "TPM",
        cluster_columns = TRUE,
        show_column_dend = TRUE,
        show_row_dend = TRUE,
        cluster_rows=TRUE,
        row_names_gp = gpar(fontsize = 4),
        column_title_gp = gpar(fontsize = 2),
        show_column_names = TRUE,
        use_raster = TRUE,
        column_title_side = "top",
        raster_quality = 4)
print(dim(mat))
pdf(file.path(opt$Outdir, "tf_heatmap.pdf"))
print(p)
dev.off()