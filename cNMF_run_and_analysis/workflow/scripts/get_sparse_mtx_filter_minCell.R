library(SeuratObject)
library(Seurat)
library(Matrix)
library(optparse)
library(dplyr)
library(tidyr)
option.list <- list(
    make_option("--min_cell", type="numeric"),
    make_option("--outdir", type="character"),
    make_option("--scratch", type="character"),
    make_option("--gene_df", type="character")
)
opt <- parse_args(OptionParser(option_list=option.list))

#########################################################
s.update.mtx=Read10X(file.path(opt$scratch),gene.column=1)
s.update=CreateSeuratObject(s.update.mtx,min.cells=opt$min_cell)
print(dim(s.update))
#########################################################
umi_mtx <- GetAssayData(object = s.update, slot = "counts", assay="RNA")
mtx <- umi_mtx
writeMM(obj=mtx,file= file.path(opt$outdir, "count_matrix/matrix.mtx"))
write(x = rownames(mtx), file = file.path(opt$outdir, "count_matrix/features.tsv"))
write(x = colnames(mtx), file = file.path(opt$outdir, "count_matrix/barcodes.tsv"))
df <- read.table(opt$gene_df, header=T, sep="\t", stringsAsFactors=F)
df <- df %>% filter(Gene %in% rownames(s.update))
saveRDS(df, file.path(opt$outdir, "barcode.dict.RDS"))


