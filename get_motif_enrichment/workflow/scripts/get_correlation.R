library(tidyverse)
library(optparse)
library(data.table)
option.list <- list(
    make_option("--filtered_gep", type="character"),
    make_option("--motif_mtx", type="character"),
    make_option("--outdir", type="character"),
    make_option("--GEPBeta", type="character")
)

opt <- parse_args(OptionParser(option_list=option.list))
usage <- read.table(opt$filtered_gep, header=T, row.names=1, sep="\t",stringsAsFactors=F)
motif_mtx <- fread(opt$motif_mtx, header=T, sep="\t",stringsAsFactors=F) %>% column_to_rownames(var="V1")
rownames(motif_mtx)=sub("_[^_]+$", "",rownames(motif_mtx))
motif_mtx=as.data.frame(t(as.matrix(motif_mtx)))
motif_list=colnames(motif_mtx)
gene_beta_metadata <- fread(opt$GEPBeta, header=T, sep="\t", stringsAsFactors=F) 
gene_beta_metadata$Gene=str_replace_all(gene_beta_metadata$Gene, pattern = "\\.", replacement = "")
result_df <- data.frame(matrix(NA, nrow=0, ncol=7))
colnames(result_df) <- c("GEP", "Motif", "CellNumber", "P", "Estimate", "TF_rank_in_GEP","P.adjust_bonferroni")

for (gep in colnames(usage)){
    gene_beta_metadata_subset <- gene_beta_metadata %>% select(c(Gene, !!sym(gep))) %>% filter(!is.na(!!sym(gep)))%>% mutate(rank=dense_rank(-!!sym(gep)))
    subset_usage=usage %>% filter(!!sym(gep)!=0) %>% rownames_to_column(var="cbc") %>% rename(gep_usage=!!sym(gep)) %>% select(cbc, gep_usage)
    gep_result_df <- data.frame(matrix(NA, nrow=0, ncol=6))
    colnames(gep_result_df) <- c("GEP", "Motif","CellNumber", "P", "Estimate", "TF_rank_in_GEP")
    for (motif in colnames(motif_mtx)){
        motif_scores=motif_mtx %>% select(!!sym(motif)) %>% rownames_to_column(var="cbc") %>% filter(cbc %in% subset_usage$cbc) %>% rename(motif_z_score=!!sym(motif))
        suppressMessages(subset_usage_motif<-subset_usage %>% inner_join(motif_scores))
        correlation=cor.test(subset_usage_motif$gep_usage, subset_usage_motif$motif_z_score, method="spearman")
        ncell=nrow(subset_usage_motif)
        if (nrow(gene_beta_metadata_subset %>% filter(Gene==motif)!=0)){
                tf_rank=as.character(gene_beta_metadata_subset %>% filter(Gene==motif) %>% pull(rank))
        } else {
                tf_rank=NA
        }
        gep_result_df[nrow(gep_result_df) + 1,] = c(gep, motif,ncell, correlation$p.value, correlation$estimate, tf_rank)
    }
    gep_result_df$P.adjust_bonferroni <- p.adjust(gep_result_df$P, method="bonferroni")
    result_df=rbind(result_df, gep_result_df)
}
write.table(result_df,gzfile(file.path(opt$outdir, "enhancer_motif_correlation.tsv.gz")), sep="\t", quote=F, row.names=F)
result_df_filtered <- result_df %>% filter(P.adjust_bonferroni<0.05) %>% arrange(GEP)
write.table(result_df_filtered %>% group_by(GEP) %>% arrange(desc(Estimate), TF_rank_in_GEP, .by_group = TRUE),file.path(opt$outdir, "filtered_enhancer_motif_correlation.tsv"), sep="\t", quote=F, row.names=F)

