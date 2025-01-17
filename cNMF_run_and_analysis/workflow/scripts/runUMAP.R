library(SeuratObject)
library(Seurat)
library(optparse)
set.seed(42)
option.list <- list(
    make_option("--SeuratObject", type="character"),
    make_option("--outdir", type="character"),
    make_option("--dims", type="numeric")
)
opt <- parse_args(OptionParser(option_list=option.list))

s <- readRDS(opt$SeuratObject)
if (is.null(s@reductions$umap)){
    s <- RunUMAP(s, dims=1:opt$dims)
}
saveRDS(s, file.path(opt$outdir,"seurat_umap.RDS"))

