suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))

option.list <- list(
	make_option("--E2GOverlap", type="character"),
	make_option("--GEPInCells", type="character"),
	make_option("--GEPTopGene", type="character"),
	make_option("--EnhancerUniqueness", type="character"),
    make_option("--LDSC_enrichment", type="character"),
    make_option("--scoreCol", type="character"),
    make_option("--outdir", type="character"),
    make_option("--complete", type="logical"),
    make_option("--KnownGenes", type="character"),
    make_option("--GEP_TF_regulators", type="character"),
    make_option("--E2GFilter", type="numeric"),
    make_option("--variantList", type="character")
    )

opt <- parse_args(OptionParser(option_list=option.list))

e2g <- read.table(opt$E2GOverlap, header=T, sep="\t", stringsAsFactors=F) 
varList <- read.table(opt$variantList, header=T, sep="\t", stringsAsFactors=F) %>% select(rsid) %>% distinct()

e2g <- e2g %>% select(chr, start, end, LeadVariant, distance, !!sym(opt$scoreCol), CellType, TargetGene) %>% rename(enh_chr=chr, enh_start=start, enh_end=end, Celltype=CellType) %>% distinct() %>% group_by(enh_chr, enh_start, enh_end,  LeadVariant, Celltype) %>% mutate(TargetGeneRank=dense_rank(-!!sym(opt$scoreCol))) %>% ungroup() %>% filter(TargetGeneRank<=opt$E2GFilter)
print(colnames(e2g)) 

gepCelltype <- read.table(opt$GEPInCells, header=T, sep="\t", stringsAsFactors=F, fill=TRUE) %>% select(c(fine_celltype, GEP, Proportion_in_GEP)) %>% distinct() %>% filter(Proportion_in_GEP > 0.1)

enh_df <- e2g %>% left_join(gepCelltype, by=c("Celltype"="fine_celltype"))
print(colnames(enh_df)) 
########################################################################
#enh_df_gep=enh_df %>% separate_longer_delim(GEP, delim="|") %>% mutate(GEP=ifelse(GEP=="", NA, GEP))
enh_df_gep=enh_df 
gep_list=na.omit(unique(enh_df_gep$GEP))
gene_list=unique(enh_df_gep$TargetGene)
########################################################################
gene_beta_metadata <- read.table(opt$GEPTopGene, header=T, sep="\t", stringsAsFactors=F) %>% select(c(Gene, GEP,Rank)) %>% distinct() %>% rename(GeneRankInGEP=Rank)

gep_tf <- read.table(opt$GEP_TF_regulators, header=T, sep="\t", stringsAsFactors=F) %>% select(GEP, Motif) %>% group_by(GEP) %>% mutate(TFs=paste(sort(unique(Motif)), collapse=",")) %>% select(GEP, TFs) %>% distinct()

enh_df_gep <- enh_df_gep %>% inner_join(gene_beta_metadata, by=c("TargetGene"="Gene", "GEP")) %>% left_join(gep_tf ) %>% select(LeadVariant, Celltype, TargetGene, GEP,GeneRankInGEP, Proportion_in_GEP, TFs)%>% group_by(LeadVariant, Celltype, TargetGene) %>% arrange(desc(Proportion_in_GEP)) %>% mutate(GEP=paste(unique(GEP), collapse="|")) %>% mutate(GeneRankInGEP=paste(unique(GeneRankInGEP),collapse="|")) %>% mutate(Proportion_in_GEP=paste(unique(Proportion_in_GEP), collapse="|")) %>% mutate(TFs=paste(unique(TFs), collapse="|")) %>% distinct()
print(head(enh_df_gep))
enh_df <- enh_df %>% select(-c(GEP, Proportion_in_GEP)) %>% left_join(enh_df_gep)
########################################################################
enh_df <- enh_df %>% arrange(LeadVariant, -!!sym(opt$scoreCol))
print(head(enh_df))
########################################################################
knownGenes <- read.table(opt$KnownGenes, header=T, sep="\t", stringsAsFactors=F) 
enh_df <- enh_df %>% left_join(knownGenes, by=c("TargetGene"="Gene"))
########################################################################
write.table(enh_df, file.path(opt$outdir, "Credibleset_gene_variant_info_E2G_with_info.tsv"), sep="\t", quote=F, row.names=F)
if (opt$complete){
    enh_df_complete=enh_df[complete.cases(enh_df),]
} else {
    enh_df_complete=enh_df %>%filter(!is.na(GEP))
}
write.table(enh_df_complete %>% distinct(), file.path(opt$outdir, "Credibleset_gene_variant_info_E2G_with_info_complete.tsv"), sep="\t", quote=F, row.names=F)


