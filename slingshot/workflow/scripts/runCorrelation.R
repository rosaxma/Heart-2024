suppressPackageStartupMessages(library(optparse))
library(tidyverse)
library(scales)
library(ggrastr)
set.seed(42)
option.list <- list(
    make_option("--LineageTable", type="character", default=""),
    make_option("--GEP_table", type="character", default=""),
    make_option("--outdir", type="character", default=""),
    make_option("--sampleInfo", type="character", default=NA)
)
opt <- parse_args(OptionParser(option_list=option.list))

lineageTable <- read.table(opt$LineageTable, header=T, sep="\t", row.names=1, stringsAsFactors=F)
cbc_table <- read.table(opt$sampleInfo, header=T, sep="\t", stringsAsFactors=F) %>% select(sample_orig, GEX_Barcode, ATAC_barcode, Rounded.PCW) %>% mutate(gex_cbc=paste0(sample_orig, "_", GEX_Barcode), atac_cbc=paste0(sample_orig, "_",ATAC_barcode)) %>% select(gex_cbc, atac_cbc, Rounded.PCW) 
GEP_table <- read.table(opt$GEP_table, header=T, sep="\t", stringsAsFactors=F, row.names=1) %>% rownames_to_column(var="gex_cbc") %>% left_join(cbc_table, by=c("gex_cbc")) %>% select(-gex_cbc) %>% column_to_rownames(var="atac_cbc")
GEP_column_names = colnames(GEP_table)


mat = matrix(ncol = 5, nrow = 0) 
cor_df=data.frame(mat)
colnames(cor_df) <- c("Lineage", "GEP", "Correlation", "P", "NumberOfCells")

for (lineage in colnames(lineageTable)){
  sub_lineageTable <- lineageTable %>% select(!!sym(lineage)) %>% merge(GEP_table,by="row.names", all=TRUE) %>% filter(!is.na(!!sym(lineage)))
  print(head(sub_lineageTable))
  for (gep in GEP_column_names){
      pseudotime=sub_lineageTable %>% filter(!is.na(!!sym(gep))) %>% pull(!!sym(lineage))
      gep_expression=sub_lineageTable %>% filter(!is.na(!!sym(gep))) %>% pull(!!sym(gep))
      if (length(pseudotime)!=0) {
      cor.result=cor.test(pseudotime, gep_expression, method=c("spearman"))
      pvalue <- cor.result$p.value
      estimate <- cor.result$estimate
      if (abs(estimate) > 0.1){
          cor_df[nrow(cor_df) + 1,] = c(lineage,gep, estimate, pvalue, length(pseudotime))
      }
    }
  }
}
if (nrow(cor_df) > 0 ){ 
	write.table(cor_df, file.path(opt$outdir, "pseudotime_correlation.tsv"), sep="\t", row.names=F, quote=F)
}

