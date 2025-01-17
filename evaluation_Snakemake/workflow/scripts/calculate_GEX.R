suppressPackageStartupMessages(library(optparse))
option.list <- list(
    make_option("--outdir", type="character"), 
    make_option("--seurat", type="character"),
    make_option("--res_range", type="character"),
    make_option("--genelist", type="character"),
    make_option("--increment", type="character")
)

opt <- parse_args(OptionParser(option_list=option.list))

library(Seurat)
library(SeuratObject)
library(utils)
suppressPackageStartupMessages(library(tidyverse))

start=as.numeric(unlist(strsplit(opt$res_range, split=","))[1])
print(start)
end=as.numeric(unlist(strsplit(opt$res_range, split=","))[2])
print(end)
res_list <- seq(start, end, by=as.numeric(opt$increment))
print(res_list)

s <- readRDS(opt$seurat)
metadata <- s@meta.data
for (res in res_list){
        cluster_col <- paste0("res_", res)
        cluster_list <- metadata %>% select(!!sym(cluster_col)) %>% distinct() %>% pull(!!sym(cluster_col))
        Idents(s) <- cluster_col
        library(Matrix)
        counter=0
        all_counter=0

        all_df <- 0
        #all_df_var <- 0
        gene_list <- read.table(opt$genelist, sep="\t", header=FALSE, stringsAsFactors=FALSE)
        print(cluster_list)
        for (cluster in cluster_list){
                all_counter=all_counter+1
                gene_matrix <- GetAssayData(object = subset(s, idents=cluster), slot = "counts")
                print(object.size(gene_matrix), units = "auto", standard = "IEC")
                gene_sum <- Matrix::rowSums(gene_matrix)
                #gene_var <- apply(gene_matrix, 1, var)
                #gex_var_df <- data.frame(genes=names(gene_var), var=gene_var)
                gex_df <- data.frame(genes=names(gene_sum), nUMI=gene_sum)
                print(cluster)
                norm_gex_df <- gex_df %>% mutate(nUMI=nUMI/sum(gene_sum)*1000000) %>% rename(!!sym(paste0(cluster,"_transcript_per_million")):="nUMI")
                #gex_var_df <- gex_var_df %>% rename(!!sym(paste0(cluster, "_var")) := "var")
                if (all_counter ==1) {
                        all_df <- norm_gex_df
                        #all_df_var <- gex_var_df
                } else {
                        all_df <- all_df %>% left_join(norm_gex_df, by="genes")
                        #all_df_var <- all_df_var %>% left_join(gex_var_df, by="genes")
                }
                print(object.size(all_df), units = "auto", standard = "IEC")
                print(head(all_df))
        
        }
        write.table(all_df, file.path(opt$outdir,paste0("GEX_TPM_res", res, ".tsv")), sep="\t", quote=FALSE, row.names=FALSE)
        #write.table(all_df_var, file.path(opt$outdir,paste0("GEX_variance_res", res, ".tsv")), sep="\t", quote=FALSE, row.names=FALSE)
        annotation_IDs=annotation <-  colnames(s@meta.data)[str_detect(colnames(s@meta.data),"_annotation")]
        pdf(file.path(opt$outdir,paste0("GEX_UMAP_annotation.pdf")))
        for (i in annotation_IDs){
                p <- DimPlot(s,raster=TRUE, group.by=i) + ggtitle(i)
                print(p)
        }
        dev.off()

        if (length(cluster_list)>1){
                library(ComplexHeatmap)
                gene_df <- all_df %>% filter(genes %in% gene_list[,1])
                rownames(gene_df) <- gene_df$genes
                gene_df <- gene_df[order(row.names(gene_df)), ]
                mat <- gene_df %>% select(-genes) %>% drop_na() %>% as.matrix()
                colnames(mat) <- sub("_transcript_per_million", "", colnames(mat))
                mat <- t(scale(t(mat)))
                mat <- mat[!rowSums(is.na(mat)),]
                Seurat::BlueAndRed()
                quantile_table=quantile(mat, c(0.1, 0.95), na.rm=TRUE)
                quantile_low=as.numeric(quantile_table["10%"])
                quantile_high=as.numeric(quantile_table["95%"])
                col_fun = circlize::colorRamp2(c(quantile_low, 0, quantile_high), c("#313695", "white", "#A50026"))

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
                        column_title_side = "top",
                        raster_quality = 4)
                print(dim(mat))
                pdf(file.path(opt$outdir,paste0("GEX_heatmap_res", res, ".pdf")), height=nrow(mat)*0.2, width=ncol(mat)*0.3+2)
                print(p)      
                dev.off()
                pdf(file.path(opt$outdir,paste0("GEX_UMAP_res", res, ".pdf")))
                print(DimPlot(s,raster=TRUE,group.by=paste0("res_",res)))
                      
                dev.off()

        } else {
                pdf(file.path(opt$outdir,paste0("GEX_heatmap_res", res, ".pdf")))
                dev.off()
        }
}
