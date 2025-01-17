suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))

option.list <- list(
    make_option("--sumstats_table", type="character"),
    make_option("--traits", type="character"),
    make_option("--results_dir", type="character"), 
    make_option("--celltypes", type="character"),
    make_option("--outdir", type="character"),
    #make_option("--lead_variant_count", type="character", default=NA),
    make_option("--infoSheet", type="character", default=NA),
    make_option("--annotation", type="character", default=NA),
    make_option("--leadSNPCount", type="character", default=NA)
)
opt <- parse_args(OptionParser(option_list=option.list))
celltypes=unlist(strsplit(opt$celltypes, ","))
traits=unlist(strsplit(opt$traits, ","))
df_list <- vector(mode = "list", length=length(traits))
df_full_list <-vector(mode = "list", length=length(traits))
df_full_enrichment_v1_list  <- vector(mode = "list", length=length(traits))
df_full_enrichment_v2_list  <- vector(mode = "list", length=length(traits))
df_full_enrichment_v3_list  <- vector(mode = "list", length=length(traits))
df_full_zscore_list  <- vector(mode = "list", length=length(traits))
df_full_logP_list  <- vector(mode = "list", length=length(traits))
df_full_P_list  <- vector(mode = "list", length=length(traits))

sumstats_table <- read.table(opt$sumstats_table, header=T, sep="\t", stringsAsFactors=F) %>% mutate(DZ=str_replace(basename(StanfordSummaryStats), ".sumstats.gz", "")) %>% select(Disease, DZ)

lead_variant_count <- read.table(opt$leadSNPCount, header=T, sep="\t", stringsAsFactors=F)
sufficient_variants <- lead_variant_count %>% filter(n_LeadSNPs>=10) %>% pull(trait)

for (i in 1:length(traits)){
    trait=traits[[i]]
    result_path=file.path(opt$results_dir, paste0(trait, ".results"))
    result_df=read.table(result_path, header=T, sep="\t", stringsAsFactors=F)
    result_df$Enrichment_p_fdr <- round(p.adjust(result_df$Enrichment_p, "bonferroni"),3)
    results_df_full_table=result_df %>% select(Category, Enrichment, Enrichment_std_error, Enrichment_p, Enrichment_p_fdr, Coefficient_z.score) %>% mutate(Category=str_replace(Category, "L2_0", "")) %>% filter(Category %in% celltypes) %>% mutate(DZ=trait) %>% inner_join(sumstats_table, by=c("DZ")) %>% rename(finemapped_trait=Disease)	%>% left_join(lead_variant_count, by=c("finemapped_trait"="trait"))
    df_full_list[[i]]=results_df_full_table
    
    results_df_subset=results_df_full_table %>% filter(Enrichment>1, Enrichment_p < 0.05) %>% inner_join(sumstats_table, by=c("DZ"))  %>% select(-DZ)
    
    df_list[[i]]=results_df_subset

    results_enrichment_subset_v1=result_df %>% select(Category, Enrichment, Enrichment_p) %>% mutate(Category=str_replace(Category, "L2_0", "")) %>% filter(Category %in% celltypes) %>% mutate(DZ=trait) %>% inner_join(sumstats_table, by=c("DZ")) %>% rename(finemapped_trait=Disease) %>% select(-DZ) %>% mutate(Enrichment=ifelse(Enrichment<1,1, Enrichment)) %>% select(-c(Enrichment_p))

    results_enrichment_subset_v2=result_df %>% select(Category, Enrichment, Enrichment_p_fdr) %>% mutate(Category=str_replace(Category, "L2_0", "")) %>% filter(Category %in% celltypes) %>% mutate(DZ=trait) %>% inner_join(sumstats_table, by=c("DZ")) %>% rename(finemapped_trait=Disease) %>% select(-DZ) %>% mutate(Enrichment=ifelse(Enrichment<1 |Enrichment_p_fdr>=0.05,1, Enrichment)) %>% select(-c(Enrichment_p_fdr))
    
    results_enrichment_subset_v3=result_df %>% select(Category, Enrichment, Enrichment_p_fdr) %>% mutate(Category=str_replace(Category, "L2_0", "")) %>% filter(Category %in% celltypes) %>% mutate(DZ=trait) %>% inner_join(sumstats_table, by=c("DZ")) %>% rename(finemapped_trait=Disease) %>% select(-DZ) %>% select(-c(Enrichment_p_fdr))

    results_zscore_subset=result_df %>% select(Category, Coefficient_z.score) %>% mutate(Category=str_replace(Category, "L2_0", "")) %>% filter(Category %in% celltypes) %>% mutate(DZ=trait) %>% inner_join(sumstats_table, by=c("DZ")) %>% rename(finemapped_trait=Disease) %>% select(-DZ)

    results_logP_subset=result_df %>% select(Category, Enrichment_p) %>% mutate(Category=str_replace(Category, "L2_0", "")) %>% filter(Category %in% celltypes) %>% mutate(DZ=trait) %>% inner_join(sumstats_table, by=c("DZ")) %>% rename(finemapped_trait=Disease) %>% select(-DZ) %>% mutate(logP=-log10(Enrichment_p)) %>% select(-Enrichment_p)

    
    results_P_subset=result_df %>% select(Category, Enrichment_p_fdr) %>% mutate(Category=str_replace(Category, "L2_0", "")) %>% filter(Category %in% celltypes) %>% mutate(DZ=trait) %>% inner_join(sumstats_table, by=c("DZ")) %>% rename(finemapped_trait=Disease) %>% select(-DZ) 

    df_full_enrichment_v1_list[[i]]=results_enrichment_subset_v1
    df_full_enrichment_v2_list[[i]]=results_enrichment_subset_v2
    df_full_enrichment_v3_list[[i]]=results_enrichment_subset_v3
    df_full_zscore_list[[i]]=results_zscore_subset
    df_full_logP_list[[i]]=results_logP_subset
    df_full_P_list[[i]]=results_P_subset 
}

infoSheet <- read.table(opt$infoSheet, header=T, sep="\t", stringsAsFactors=F) %>% select(Trait, NiceName)
annotation <- read.table(opt$annotation, header=T, sep="\t", stringsAsFactors=F) %>% select(celltype, plot_label) %>% distinct()
combined_df <- do.call("rbind", df_list) %>% left_join(infoSheet, by=c("finemapped_trait"="Trait"))%>% left_join(annotation, by=c("Category"="celltype")) %>% rename(coarse_cluster=plot_label) 
combined_df_full <- do.call("rbind", df_full_list) %>% left_join(infoSheet, by=c("finemapped_trait"="Trait"))%>% left_join(annotation, by=c("Category"="celltype")) %>% rename(coarse_cluster=plot_label) 


if (is.na(opt$leadSNPCount)){
    robust_traits=combined_df %>% filter(Enrichment_p_fdr<0.05) %>% pull(finemapped_trait) %>% unique()
    print(robust_traits)
} else {
    robust_traits=combined_df %>% filter(Enrichment_p_fdr<0.05, finemapped_trait %in% sufficient_variants) %>% pull(finemapped_trait) %>% unique()
    print(robust_traits)
}


combined_enrichment_v1_mtx <- do.call("rbind", df_full_enrichment_v1_list) %>% filter(finemapped_trait %in% robust_traits ) %>% pivot_wider(names_from=finemapped_trait, values_from=Enrichment)
combined_enrichment_v2_mtx <- do.call("rbind", df_full_enrichment_v2_list) %>% filter(finemapped_trait %in% robust_traits ) %>% pivot_wider(names_from=finemapped_trait, values_from=Enrichment)
combined_enrichment_v3_mtx <- do.call("rbind", df_full_enrichment_v3_list) %>% filter(finemapped_trait %in% robust_traits ) %>% pivot_wider(names_from=finemapped_trait, values_from=Enrichment)
combined_zscore_mtx <- do.call("rbind", df_full_zscore_list) %>% filter(finemapped_trait %in% robust_traits ) %>% pivot_wider(names_from=finemapped_trait, values_from=Coefficient_z.score)
combined_logP_mtx <- do.call("rbind", df_full_logP_list) %>% filter(finemapped_trait %in% robust_traits ) %>% pivot_wider(names_from=finemapped_trait, values_from=logP)
combined_P_mtx <- do.call("rbind", df_full_P_list) %>% filter(finemapped_trait %in% robust_traits) %>% pivot_wider(names_from=finemapped_trait, values_from=Enrichment_p_fdr)

write.table(combined_df_full, file.path(opt$outdir, "combined_LDSC_enrichment_unfiltered.tsv"), sep="\t", row.names=F, quote=F)
write.table(combined_df, file.path(opt$outdir, "combined_LDSC_enrichment.tsv"), sep="\t", row.names=F, quote=F)
write.table(combined_df %>% filter(Enrichment_p_fdr<0.05), file.path(opt$outdir, "combined_LDSC_enrichment_fdr_filtered.tsv"), sep="\t", row.names=F, quote=F)
write.table(combined_enrichment_v1_mtx, file.path(opt$outdir, "combined_LDSC_enrichment_below1_set_to_1.mtx"), sep="\t", row.names=F, quote=F)
write.table(combined_enrichment_v2_mtx, file.path(opt$outdir, "combined_LDSC_enrichment_below1_or_insig_set_to_1.mtx"), sep="\t", row.names=F, quote=F)
write.table(combined_enrichment_v3_mtx, file.path(opt$outdir, "combined_LDSC_enrichment.mtx"), sep="\t", row.names=F, quote=F)
write.table(combined_zscore_mtx,file.path(opt$outdir, "combined_LDSC_zscore.mtx"), sep="\t", row.names=F, quote=F)
write.table(combined_logP_mtx, file.path(opt$outdir, "combined_LDSC_logP.mtx"), sep="\t", row.names=F, quote=F) 
write.table(combined_P_mtx, file.path(opt$outdir, "combined_LDSC_FDR_P.mtx"), sep="\t", row.names=F, quote=F) 
