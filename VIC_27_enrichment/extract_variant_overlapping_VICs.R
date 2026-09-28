library(tidyverse)
library(data.table)

valve_traits <- read.table("config/valve_traits.tsv", sep = "\t", stringsAsFactors = FALSE, header = TRUE) %>% pull(traits)

causal_genes <- read.table("config/combined_disease_causal_genes.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>% group_by(Gene) %>% summarize(Traits = paste(sort(unique(Trait)), collapse = "|")) %>% select(Gene, Traits)

v2g <- fread("config/combined_Credibleset_gene_variant_info_E2G_with_info_allV2G.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>% 
  select(Trait, CredibleSet, GEP, TargetGene, Celltype) %>% 
  distinct() %>% 
  separate_longer_delim(GEP, "|") %>% 
  filter(str_detect(Celltype, "VIC")) %>% 
  filter(Trait %in% valve_traits)

gep_top <- read.table("config/top_gene_meta_long.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

bg_gene_list <- read.table("config/expressed_cardiac_genes.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>% 
  filter(str_detect(celltype, "VIC")) %>% 
  pull(Gene) %>% 
  unique()

gep_list <- unique(gep_top$GEP)

V2G_all <- v2g %>% pull(TargetGene) %>% unique()
V2G_all <- intersect(V2G_all, bg_gene_list)

results_list <- list()

for (gep in gep_list) {
  GEP_genes <- gep_top %>% 
    filter(GEP == gep) %>% 
    pull(Gene) %>% 
    unique()
  GEP_genes <- intersect(GEP_genes, bg_gene_list)
  
  V2G_G2P <- v2g %>% 
    filter(GEP == gep) %>% 
    pull(TargetGene) %>% 
    unique()
  V2G_G2P <- intersect(V2G_G2P, bg_gene_list)
  
  in_V2G_in_GEP    <- length(V2G_G2P)
  not_V2G_in_GEP   <- length(setdiff(GEP_genes, V2G_G2P))
  in_V2G_not_GEP   <- length(setdiff(V2G_all, V2G_G2P))
  not_V2G_not_GEP  <- length(setdiff(bg_gene_list, union(V2G_all, GEP_genes)))
  
  dat <- matrix(
    c(in_V2G_in_GEP, not_V2G_in_GEP,
      in_V2G_not_GEP, not_V2G_not_GEP),
    nrow = 2, ncol = 2,
    dimnames = list(
      V2G = c("in_V2G", "not_in_V2G"),
      GEP = c("in_GEP", "not_in_GEP")
    )
  )
  
  dat_test <- fisher.test(dat, alternative = "greater")
  
  results_list[[gep]] <- data.frame(
    GEP             = gep,
    n_overlap       = in_V2G_in_GEP,
    overlap_genes   = paste(sort(V2G_G2P), collapse = ", "),
    in_V2G_in_GEP   = in_V2G_in_GEP,
    not_V2G_in_GEP  = not_V2G_in_GEP,
    in_V2G_not_GEP  = in_V2G_not_GEP,
    not_V2G_not_GEP = not_V2G_not_GEP,
    odds_ratio      = unname(dat_test$estimate),
    ci_lower        = dat_test$conf.int[1],
    ci_upper        = dat_test$conf.int[2],
    p_value         = dat_test$p.value,
    stringsAsFactors = FALSE
  )
}

results_df <- bind_rows(results_list) %>% 
  mutate(FDR = p.adjust(p_value, method = "BH")) %>% 
  arrange(p_value)

write.table(results_df, "valve_trait_genes_enrichment_in_GEP.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
