library(tidyverse)
celltype2include=read.table("assignment.tsv", header=T, sep="\t", fill=TRUE, comment.char = "") %>% filter(ConsiderForCellTypeEnrichment) %>% select(Celltype, New_Name)
write.table(celltype2include, "assignment.new.tsv", sep="\t", row.names=F, quote=F)
