suppressPackageStartupMessages(library(optparse))

option.list <- list(
    make_option("--outdir", type="character", default="", help="Output directory"), 
    make_option("--SeuratObject", type="character", default="", help="")
    )

opt <- parse_args(OptionParser(option_list=option.list))
library(SeuratObject)
library(Seurat)
library(tidyverse)
library(Matrix)
set.seed(42)
s <- readRDS(opt$SeuratObject)


umi_mtx <- GetAssayData(object = s, slot = "counts", assay="RNA")
mtx <- Matrix(umi_mtx, sparse=T)
writeMM(obj=mtx,file= file.path(opt$outdir, "matrix.mtx"))
write(x = rownames(mtx), file = file.path(opt$outdir, "features.tsv"))
write(x = colnames(mtx), file = file.path(opt$outdir, "barcodes.tsv"))
