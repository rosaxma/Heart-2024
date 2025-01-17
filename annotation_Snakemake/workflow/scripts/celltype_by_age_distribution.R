suppressPackageStartupMessages(library(optparse))
option.list <- list(
    make_option("--outdir", type="character"),
    make_option("--seurat", type="character"),
    make_option("--SampleInfo", type="character"),
    make_option("--cluster_col", type="character")
)
set.seed(42)
opt <- parse_args(OptionParser(option_list=option.list))
library(Seurat)
library(tidyverse)
s <- readRDS(opt$seurat)
metadata <- s@meta.data
sampleInfo <- read.table(opt$SampleInfo, header=T, sep="\t", stringsAsFactors=F)
sample_table <- metadata %>% select(cbc,source) %>% left_join(sampleInfo, by="source")
###############################
## Distribution table
distribution <- function(metadata, reference, target, outdir, name){
	count_table <- metadata %>% select(!!sym(reference), !!sym(target)) %>% group_by(!!sym(reference), !!sym(target)) %>% summarize(n=n()) %>% ungroup() %>% distinct() %>% pivot_wider(names_from=!!sym(target), values_from=n) %>% column_to_rownames(var=opt$cluster_col) %>% as.data.frame()
    write.table(count_table, file.path(outdir, paste0(reference, "_vs_",target, name,".tsv")), sep="\t", quote=F, row.names=T, col.names=NA)
}
###############################
## without consolidation by age
summary_table <- sample_table %>% group_by(source,Rounded.PCW) %>% summarize(n=n()) %>% ungroup() %>% filter(n>8000) %>% column_to_rownames('source')
min_cell_count=min(summary_table$n)
subsampled_cbc=c()
for (sample in rownames(summary_table)){
    if (summary_table[sample,"n"] == min_cell_count){
        subsampled_cbc=c(subsampled_cbc,metadata %>% filter(source==sample) %>% pull(cbc))
    } else{
        original_cbc=metadata %>% filter(source==sample) %>% pull(cbc)
        cbc=sample(original_cbc,size=min_cell_count,replace=FALSE)
        subsampled_cbc=c(subsampled_cbc,cbc)
    }
}
s_subset=subset(s,cells=subsampled_cbc)
metadata_subset <- s_subset@meta.data %>% left_join(sampleInfo %>% select(source,Rounded.PCW)) %>% mutate(Rounded.PCW_source=paste0(Rounded.PCW, "_",source))%>% mutate(Rounded.PCW_Phase=paste0(Rounded.PCW, "_",Phase))

metadata_subset_summary= metadata_subset %>% group_by(Rounded.PCW_source) %>% summarize(n=n()) %>% ungroup()

write.table(metadata_subset_summary, file.path(opt$outdir, "subset_metadata_summary.tsv"), sep="\t", row.names=F, quote=F)
target_features=c("Rounded.PCW", "Rounded.PCW_source","Rounded.PCW_Phase")
for (target_feature in target_features){
	distribution(metadata_subset, opt$cluster_col, target_feature, opt$outdir, "_subsample_by_sample")
}
###############################
##consolidate and subsample by age
summary_table <- sample_table %>% group_by(Rounded.PCW) %>% summarize(n=n()) %>% ungroup()%>% column_to_rownames('Rounded.PCW')
min_cell_count=min(summary_table$n)
subsampled_cbc=c()
metadata=metadata %>% left_join(sampleInfo %>% select(source,Rounded.PCW))
for (sample in rownames(summary_table)){
    if (summary_table[sample,"n"] == min_cell_count){
        subsampled_cbc=c(subsampled_cbc,metadata %>% filter(Rounded.PCW==sample) %>% pull(cbc))
    } else{
        original_cbc=metadata %>% filter(Rounded.PCW==sample) %>% pull(cbc)
        cbc=sample(original_cbc,size=min_cell_count,replace=FALSE)
        subsampled_cbc=c(subsampled_cbc,cbc)
    }
}
print(dim(s))
print(length(subsampled_cbc))
print(length(intersect(subsampled_cbc, colnames(s))))
s_subset=subset(s,cells=subsampled_cbc)
print(dim(s_subset))
metadata_subset <- s_subset@meta.data %>% left_join(sampleInfo %>% select(source,Rounded.PCW)) %>% mutate(Rounded.PCW_Phase=paste0(Rounded.PCW, "_",Phase))
target_features=c("Rounded.PCW","Rounded.PCW_Phase")
for (target_feature in target_features){
	distribution(metadata_subset, opt$cluster_col, target_feature,opt$outdir,"_subsample_by_age_group")
}
