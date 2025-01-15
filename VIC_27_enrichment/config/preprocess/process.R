library(tidyverse)
library(data.table)
causal_genes <- read.table("combined_disease_causal_genes.tsv", header=T, sep="\t", stringsAsFactors=F) %>% mutate(Trait=ifelse(Trait=="Obstructive_lesions", "Valve_defects", Trait)) %>% filter(category=="CHD") %>% select(Gene, Trait) 
write.table(causal_genes, "../combined_disease_causal_genes.tsv", sep="\t", quote=F, row.names=F)

#library(tidyverse)
#ref_file="CollapsedGeneBounds.hg38.TSS500bp.bed"
#ref_table <- read.table(ref_file, header=F, sep="\t", stringsAsFactors=F)
#gene_table <- read.table("GEX_TPM.tsv.gz", header=T, sep="\t", stringsAsFactors=F)
#gene_table_long <- gene_table %>% pivot_longer(!genes, names_to = "celltype", values_to = "TPM") %>% filter(TPM>1) %>% rename(Gene=genes) %>% mutate(celltype=str_replace(celltype, "_transcript_per_million", "")) %>% filter(Gene %in% ref_table$V4)
#cardiac_celltypes <- read.table("cardiac_cells.tsv", header=T, sep="\t", stringsAsFactors=F) %>% filter(CardiacCellType) %>% pull(Celltype)
#gene_table_long <- gene_table_long %>% filter(Gene %in% ref_table$V4) %>% filter(celltype %in% cardiac_celltypes)
#write.table(gene_table_long, "../expressed_cardiac_genes.tsv", sep="\t", quote=F, row.names=F)



