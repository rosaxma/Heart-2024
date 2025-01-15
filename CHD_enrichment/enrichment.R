library(tidyverse)
library(data.table)
set.seed(42)
#######################################################
greplany <- function(patterns, v) {
  match <- rep(FALSE, length(v))
  for (pattern in patterns) {
    match <- match | grepl(pattern, v)
  }
  return(match)
}
#######################################################
celltype2include=read.table("config/assignment.tsv", header=T, sep="\t", fill=TRUE, comment.char = "") %>% pull(Celltype)%>% unique()

ref_table <- read.table("config/CollapsedGeneBounds.hg38.TSS500bp.bed", header=F, sep="\t", stringsAsFactors=F)
ref_table <- ref_table %>% filter(!greplany(c("^LINC","-AS","^MIR","RNU","^LOC","^RPS", "^RPL"),V4))

causal_genes <- read.table("config/combined_disease_causal_genes.tsv", header=T, sep="\t", stringsAsFactors=F) %>% filter(Gene %in% ref_table$V4) 
combined_CHD_genes=causal_genes %>% mutate(Trait="AllCHDGenesCombined") 
causal_genes=rbind(causal_genes,combined_CHD_genes)

gene_table_tpm <- read.table("config/GEX_TPM.tsv.gz",header=T, sep="\t", stringsAsFactors=F)
colnames_tpm=colnames(gene_table_tpm)
colnames_tpm_new=str_replace(colnames_tpm, "_transcript_per_million", "")
colnames(gene_table_tpm) <- colnames_tpm_new
#######################################################
gene_table_tpm <- gene_table_tpm %>% select(all_of(c("genes", celltype2include))) %>% filter(genes %in% ref_table$V4) %>% column_to_rownames(var="genes")
# filter for genes that have TPM>1 in at least 1 cell
gene_table_tpm=gene_table_tpm[rowSums(gene_table_tpm > 1) >0, ] %>% rownames_to_column(var="genes")
# z score normalization of TPM
gene_table_zscore=t(scale(t(gene_table_tpm %>% column_to_rownames(var="genes")))) %>% as.data.frame() %>% rownames_to_column(var="genes") %>% pivot_longer(!genes, names_to = "celltype", values_to = "Z_score") %>% dplyr::rename(Gene=genes) %>% filter(Gene %in% ref_table$V4) %>% mutate(celltype=str_replace(celltype, "_transcript_per_million", ""))

# Initialize the output table
enrichment_recall_table=data.frame(matrix(NA, nrow=0, ncol=15))
colnames(enrichment_recall_table) <- c("Category", "Disease", "Celltype", "in_Celltype_true_disease", "in_Celltype_false_disease", "not_in_Celltype_true_disease", "not_in_Celltype_false_disease", "odds_ratio", "CI_lower", "CI_upper", "P_value", "NumberOfTrueDiseaseGenes", "recall", "disease_genes_in_Celltype", "minZscore")

gene_table_long=gene_table_tpm %>% pivot_longer(!genes, names_to = "celltype", values_to = "TPM") %>% mutate(celltype=str_replace(celltype, "_transcript_per_million", "")) %>% dplyr::rename(Gene=genes) %>% filter(Gene %in% ref_table$V4) %>% filter(TPM>1)

cardiac_genes <- unique(gene_table_zscore$Gene)
tpm_zscore_df <- gene_table_long %>% left_join(gene_table_zscore %>% select(Gene, celltype, Z_score))
output_df <- tpm_zscore_df %>% select(celltype, Gene, TPM, Z_score) 
# only rank genes with Z score > 1 in each cell type
output_df_tmp <- output_df %>% filter(Z_score > 1) %>% group_by(celltype) %>% mutate(TPMRank=rank(desc(TPM), ties.method="random")) %>% ungroup()
output_df <- output_df %>% left_join(output_df_tmp)
write.table(output_df, "celltype_gene_table.tsv", sep="\t", quote=F, row.names=F)
########################################################################################
for (cat in unique(causal_genes$category)){
    for (disease in unique(causal_genes %>% filter(category==cat) %>% pull(Trait))) {
        print(disease)
        true_disease_genes=intersect(causal_genes %>% filter(Trait==disease) %>% pull(Gene) %>% unique(),cardiac_genes)
        if (length(true_disease_genes) >= 10)  {
                for (ct in unique(gene_table_zscore$celltype)) {
                    celltype_tpm_zscore_df <- tpm_zscore_df %>% filter(celltype==ct) %>% filter(Z_score > 1) %>% mutate(TPMRank=rank(desc(TPM), ties.method="random")) %>% filter(TPMRank < 301)
                    celltype_genes=unique(celltype_tpm_zscore_df$Gene)

                    in_Celltype_true_disease=length(intersect(celltype_genes, true_disease_genes))
                    in_Celltype_false_disease=length(setdiff(celltype_genes, true_disease_genes))
                    not_in_Celltype_true_disease=length(setdiff(true_disease_genes, celltype_genes))
                    bg_genes=setdiff(cardiac_genes, celltype_genes)
                    not_in_Celltype_false_disease=length(setdiff(bg_genes, true_disease_genes))

                    dat=data.frame(c(in_Celltype_true_disease, not_in_Celltype_true_disease), c(in_Celltype_false_disease,not_in_Celltype_false_disease))
                    colnames(dat)=c("true_disease", "false_disease")
                    rownames(dat)=c("in_Celltype", "not_in_Celltype")
                    row_index=nrow(enrichment_recall_table)+1
                    enrichment_recall_table[row_index, "Category"]=cat
                    enrichment_recall_table[row_index, "Disease"]=disease
                    enrichment_recall_table[row_index, "Celltype"]=ct
                    enrichment_recall_table[row_index, "in_Celltype_true_disease"]=in_Celltype_true_disease
                    enrichment_recall_table[row_index, "in_Celltype_false_disease"]=in_Celltype_false_disease
                    enrichment_recall_table[row_index, "not_in_Celltype_true_disease"]=not_in_Celltype_true_disease
                    enrichment_recall_table[row_index, "not_in_Celltype_false_disease"]=not_in_Celltype_false_disease
                    enrichment_recall_table[row_index, "NumberOfTrueDiseaseGenes"]=length(true_disease_genes)
                    enrichment_recall_table[row_index, "recall"]=round(in_Celltype_true_disease/length(true_disease_genes), digits=3)
                    enrichment_recall_table[row_index, "disease_genes_in_Celltype"]=paste(intersect(celltype_genes, true_disease_genes), collapse="|")
		                enrichment_recall_table[row_index, "minZscore"]=min(celltype_tpm_zscore_df$Z_score)
                    if (!(any(is.na(dat)))){
                            overall_enrichment_p_value <- fisher.test(dat, conf.int, alternative="greater")$p.value
                            dat_test <- fisher.test(dat, conf.int, alternative="greater")
                            enrichment_recall_table[row_index, "odds_ratio"] <- dat_test$estimate
                            enrichment_recall_table[row_index, "CI_lower"] <- dat_test$conf.int[1]
                            enrichment_recall_table[row_index, "CI_upper"] <- dat_test$conf.int[2]
                            enrichment_recall_table[row_index, "P_value"] <- dat_test$p.value
                    } else{
                    enrichment_recall_table[row_index, "odds_ratio"] <- NA
                    enrichment_recall_table[row_index, "CI_lower"] <- NA
                    enrichment_recall_table[row_index, "CI_upper"] <- NA
                    enrichment_recall_table[row_index, "P_value"] <- NA
                    }
            }
        }
    }
}

enrichment_recall_table  <- enrichment_recall_table  %>% filter(!is.na(odds_ratio)) %>% group_by(Category) %>% mutate(P_fdr=p.adjust(P_value)) %>%mutate(P_bonferroni=round(p.adjust(P_value, method="bonferroni"), digits=3)) %>%  ungroup() 
## rename
rename_df <- read.table("config/assignment.tsv", header=T, sep="\t", fill=TRUE, comment.char = "")
enrichment_recall_table <- enrichment_recall_table %>% left_join(rename_df, by=c("Celltype")) %>% select(-Celltype) %>% rename(Celltype=New_Name) %>% relocate(Celltype, .after=Disease)
write.table(enrichment_recall_table, file.path("disease_enrichment_recall.tsv"), sep="\t", quote=F, row.names=F)
write.table(enrichment_recall_table %>% filter(P_fdr < 0.05) %>% filter(CI_lower>1), file.path("fdr_filtered_disease_enrichment_recall.tsv"), sep="\t", quote=F, row.names=F)
