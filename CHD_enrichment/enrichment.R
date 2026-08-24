library(tidyverse)
library(data.table)
set.seed(42)

#######################################################
# Helper Functions
#######################################################
greplany <- function(patterns, v) {
  match <- rep(FALSE, length(v))
  for (pattern in patterns) {
    match <- match | grepl(pattern, v)
  }
  return(match)
}

#######################################################
# Initial Data Loading & Pre-processing
#######################################################
# Read cell types and rename mapping ONCE at the start
rename_df <- read.table("config/assignment.tsv", header=T, sep="\t", fill=TRUE, comment.char = "")
celltype2include <- rename_df %>% pull(Celltype) %>% unique()

# Load and filter reference table
ref_table <- read.table("config/CollapsedGeneBounds.hg38.TSS500bp.bed", header=F, sep="\t", stringsAsFactors=F)
ref_table <- ref_table %>% filter(!greplany(c("^LINC","-AS","^MIR","RNU","^LOC","^RPS", "^RPL"), V4))

# Load causal genes and add combined traits
causal_genes <- read.table("config/combined_disease_causal_genes.tsv", header=T, sep="\t", stringsAsFactors=F) %>% 
  filter(Gene %in% ref_table$V4) 
combined_CHD_genes <- causal_genes %>% mutate(Trait="AllCHDGenesCombined") 
causal_genes <- rbind(causal_genes, combined_CHD_genes)

# Load and format TPM table
gene_table_tpm <- read.table("config/GEX_TPM.tsv.gz", header=T, sep="\t", stringsAsFactors=F)
colnames_tpm_new <- str_replace(colnames(gene_table_tpm), "_transcript_per_million", "")
colnames(gene_table_tpm) <- colnames_tpm_new

#######################################################
# Gene Filtering and Z-score Normalization
#######################################################
gene_table_tpm <- gene_table_tpm %>% 
  select(all_of(c("genes", celltype2include))) %>% 
  filter(genes %in% ref_table$V4) %>% 
  column_to_rownames(var="genes")

# Filter for genes that have TPM > 1 in at least 1 cell type
gene_table_tpm <- gene_table_tpm[rowSums(gene_table_tpm > 1) > 0, ] %>% 
  rownames_to_column(var="genes")

# Z-score normalization of TPM
gene_table_zscore <- t(scale(t(gene_table_tpm %>% column_to_rownames(var="genes")))) %>% 
  as.data.frame() %>% 
  rownames_to_column(var="genes") %>% 
  pivot_longer(!genes, names_to = "celltype", values_to = "Z_score") %>% 
  dplyr::rename(Gene=genes) %>% 
  filter(Gene %in% ref_table$V4) %>% 
  mutate(celltype=str_replace(celltype, "_transcript_per_million", ""))

gene_table_long <- gene_table_tpm %>% 
  pivot_longer(!genes, names_to = "celltype", values_to = "TPM") %>% 
  mutate(celltype=str_replace(celltype, "_transcript_per_million", "")) %>% 
  dplyr::rename(Gene=genes) %>% 
  filter(Gene %in% ref_table$V4) %>% 
  filter(TPM > 1)

cardiac_genes <- unique(gene_table_zscore$Gene)

tpm_zscore_df <- gene_table_long %>% 
  left_join(gene_table_zscore %>% select(Gene, celltype, Z_score), by = c("Gene", "celltype"))

output_df <- tpm_zscore_df %>% select(celltype, Gene, TPM, Z_score) 

# Only rank genes with Z score > 1 in each cell type
output_df_tmp <- output_df %>% 
  filter(Z_score > 1) %>% 
  group_by(celltype) %>% 
  mutate(TPMRank=rank(desc(TPM), ties.method="random")) %>% 
  ungroup()

output_df <- output_df %>% left_join(output_df_tmp, by = c("celltype", "Gene", "TPM", "Z_score"))
write.table(output_df, "celltype_gene_table.tsv", sep="\t", quote=F, row.names=F)

########################################################################################
# Disease Enrichment Loop
########################################################################################
for (TPM_thresh in c(101, 201, 301, 401, 501, 1001, Inf)) {
  
  # Use a list to store rows (much faster than row-by-row rbind in R)
  results_list <- list()
  
  for (cat in unique(causal_genes$category)) {
    for (disease in unique(causal_genes %>% filter(category==cat) %>% pull(Trait))) {
      print(disease)
      
      true_disease_genes <- intersect(causal_genes %>% filter(Trait==disease) %>% pull(Gene) %>% unique(), cardiac_genes)
      num_true_disease <- length(true_disease_genes)
      
      if (num_true_disease >= 10) {
        for (ct in unique(gene_table_zscore$celltype)) {
          
          celltype_tpm_zscore_df <- tpm_zscore_df %>% 
            filter(celltype==ct) %>% 
            filter(Z_score > 1) %>% 
            mutate(TPMRank=rank(desc(TPM), ties.method="random")) %>% 
            filter(TPMRank < TPM_thresh)
          
          celltype_genes <- unique(celltype_tpm_zscore_df$Gene)
          
          in_Celltype_true_disease <- length(intersect(celltype_genes, true_disease_genes))
          in_Celltype_false_disease <- length(setdiff(celltype_genes, true_disease_genes))
          not_in_Celltype_true_disease <- length(setdiff(true_disease_genes, celltype_genes))
          bg_genes <- setdiff(cardiac_genes, celltype_genes)
          not_in_Celltype_false_disease <- length(setdiff(bg_genes, true_disease_genes))
          
          dat <- data.frame(
            true_disease = c(in_Celltype_true_disease, not_in_Celltype_true_disease), 
            false_disease = c(in_Celltype_false_disease, not_in_Celltype_false_disease),
            row.names = c("in_Celltype", "not_in_Celltype")
          )
          
          if (!(any(is.na(dat)))) {
            # Calculate fisher test ONCE per iteration
            dat_test <- fisher.test(dat, conf.int=TRUE, alternative="greater")
            
            odds_ratio <- dat_test$estimate
            CI_lower <- dat_test$conf.int[1]
            CI_upper <- dat_test$conf.int[2]
            P_value <- dat_test$p.value
          } else {
            odds_ratio <- NA
            CI_lower <- NA
            CI_upper <- NA
            P_value <- NA
          }
          
          # Safely get min Z-score to prevent Inf warnings
          min_z <- ifelse(nrow(celltype_tpm_zscore_df) > 0, min(celltype_tpm_zscore_df$Z_score), NA)
          
          # Append to list
          results_list[[length(results_list) + 1]] <- data.frame(
            Category = cat,
            Disease = disease,
            Celltype = ct,
            in_Celltype_true_disease = in_Celltype_true_disease,
            in_Celltype_false_disease = in_Celltype_false_disease,
            not_in_Celltype_true_disease = not_in_Celltype_true_disease,
            not_in_Celltype_false_disease = not_in_Celltype_false_disease,
            odds_ratio = odds_ratio,
            CI_lower = CI_lower,
            CI_upper = CI_upper,
            P_value = P_value,
            NumberOfTrueDiseaseGenes = num_true_disease,
            recall = round(in_Celltype_true_disease / num_true_disease, digits=3),
            disease_genes_in_Celltype = paste(intersect(celltype_genes, true_disease_genes), collapse="|"),
            minZscore = min_z
          )
        }
      }
    }
  }
  
  # Bind all rows for this threshold into a single dataframe
  enrichment_recall_table <- bind_rows(results_list)
  
  # Filter, calculate FDR/Bonferroni, and rename cell types
  enrichment_recall_table <- enrichment_recall_table %>% 
    filter(!is.na(odds_ratio)) %>% 
    group_by(Category) %>% 
    mutate(
      P_fdr = p.adjust(P_value, method="fdr"),
      P_bonferroni = round(p.adjust(P_value, method="bonferroni"), digits=3)
    ) %>% 
    ungroup() %>% 
    left_join(rename_df, by="Celltype") %>% 
    select(-Celltype) %>% 
    rename(Celltype=New_Name) %>% 
    relocate(Celltype, .after=Disease)
  
  # ---------------------------------------------------------
  # Write Files
  # ---------------------------------------------------------
  # Write threshold-specific files
  write.table(enrichment_recall_table, file.path(paste0("disease_enrichment_recall_thresh", as.character(TPM_thresh), ".tsv")), sep="\t", quote=F, row.names=F)
  
  write.table(
    enrichment_recall_table %>% filter(P_bonferroni < 0.05, CI_lower > 1), 
    file.path(paste0("bonferroni_filtered_disease_enrichment_recall_thresh", as.character(TPM_thresh), ".tsv")), 
    sep="\t", quote=F, row.names=F
  )
  
  # If this is the "Inf" (unrestricted) threshold, export the master files
  if (is.infinite(TPM_thresh)) {
    write.table(enrichment_recall_table, file.path("disease_enrichment_recall.tsv"), sep="\t", quote=F, row.names=F)
    
    write.table(
      enrichment_recall_table %>% filter(P_fdr < 0.05, CI_lower > 1), 
      file.path("fdr_filtered_disease_enrichment_recall.tsv"), 
      sep="\t", quote=F, row.names=F
    )
  }
}