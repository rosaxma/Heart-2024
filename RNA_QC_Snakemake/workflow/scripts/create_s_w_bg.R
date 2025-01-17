suppressPackageStartupMessages(library(optparse))

option.list <- list(
    make_option("--STARfilteredMatrices", type="character", default=""), 
    make_option("--outdir", type="character", default=1), 
    make_option("--sample", type="character", default=1) 
)
opt <- parse_args(OptionParser(option_list=option.list))
library(SeuratObject)
library(Seurat)

STAR_filtered_matrix <- Read10X(opt$STARfilteredMatrices)
s <- CreateSeuratObject(STAR_filtered_matrix, project=opt$sample)
saveRDS(s, file.path(opt$outdir, paste0(opt$sample, "_SeuratObject.RDS")))
