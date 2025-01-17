suppressPackageStartupMessages(library(optparse))
option.list <- list(
    make_option("--summary_sheet",type="character"),
    make_option("--atacInfoSheet", type="character"),
    make_option("--outdir",type="character"),
    make_option("--sample", type="character")
    )

opt <- parse_args(OptionParser(option_list=option.list))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(tidyverse))
set.seed(42)

input_doublets <- read.table(opt$summary_sheet, header=TRUE, sep = '\t', stringsAsFactors=FALSE)
atac_InfoSheet <- read.table(opt$atacInfoSheet, header=TRUE, sep = '\t', stringsAsFactors=FALSE)

atac_InfoSheet <- atac_InfoSheet %>% select(-c(AmuletMultiplets, ArchRDoublets)) %>% rename(ATAC_barcode=barcode) %>% inner_join(input_doublets, by=c("Sample"="orig.ident", "GEX_Barcodes"="cbc")) 

write.table(atac_InfoSheet, file.path(opt$outdir, paste0(opt$sample, "_ATAC_cell_quality_sheet.tsv")), row.names=F, quote=F, sep="\t")