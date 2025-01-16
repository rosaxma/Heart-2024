library(tidyverse)
library(data.table)
valve_traits <- read.table("config/valve_traits.tsv", sep="\t", stringsAsFactors=F, header=T) %>% pull(traits)

causal_genes <- read.table("config/combined_disease_causal_genes.tsv", header=T, sep="\t", stringsAsFactors=F) %>% group_by(Gene) %>% summarize(Traits=paste(sort(unique(Trait)), collapse="|")) %>% select(Gene, Traits)

###
message("number of target genes from the valve loci")
v2g <- fread("config/combined_Credibleset_gene_variant_info_E2G_with_info_allV2G.tsv", header=T, sep="\t", stringsAsFactors=F) %>% select(Trait, CredibleSet, GEP, TargetGene,Celltype)  %>% distinct() %>% separate_longer_delim(GEP, "|") %>% filter(str_detect(Celltype, "VIC")) %>% filter(Trait %in% valve_traits)

gep_list=c("VIC_GEP_27")
for (gep in gep_list){
V2G=v2g %>% pull(TargetGene) %>% unique()
V2G_G2P=v2g %>% filter(GEP==gep) %>% pull(TargetGene) %>% unique()
print(paste(V2G_G2P, collapse=","))
V2G_not_G2P=setdiff(V2G, V2G_G2P)
gep_top <- read.table("config/top_gene_meta_long.tsv", header=T, sep="\t", stringsAsFactors=F)
GEP_genes=gep_top %>% filter(GEP==gep) %>% pull(Gene)%>% unique()
not_V2G_G2P=setdiff(GEP_genes, V2G_G2P)
V2G_not_G2P=setdiff(V2G, V2G_G2P)
bg_gene_list <- read.table("config/expressed_cardiac_genes.tsv", sep="\t", header=T, stringsAsFactors=F) %>% filter(str_detect(celltype, "VIC")) %>% pull(Gene) %>% unique()
not_V2G_not_G2P=setdiff(bg_gene_list, unique(c(V2G, GEP_genes)))
dat=data.frame(c(length(V2G_G2P),length(not_V2G_G2P)), c(length(V2G_not_G2P), length(not_V2G_not_G2P)))
colnames(dat)=c("in_GEP", "not_in_GEP")
rownames(dat)=c("in_V2G", "not_in_V2G")
dat_test <- fisher.test(dat, conf.int, alternative="greater")
print(gep)
print(dat)
print(dat_test$estimate)
print(dat_test$conf.int[1])
print(dat_test$conf.int[2])
print(dat_test$p.value)
}
