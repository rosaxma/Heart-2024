library(Seurat)
library(sceasy)
library(reticulate)
library(optparse)
option.list <- list(
    make_option("--pythonPath", type="character", default=""),
    make_option("--SeuratObject", type="character", default=""),
    make_option("--out_h5ad", type="character", default=""),
    make_option("--out_high_var_genes", type="character", default=""),
    make_option("--geneList", type="character", default="")
)
opt <- parse_args(OptionParser(option_list=option.list))
use_python(opt$pythonPath)
s <- readRDS(opt$SeuratObject)
sceasy::convertFormat(s, from="seurat", to="anndata", outFile=opt$out_h5ad)
highly_variable_genes = s@assays$RNA@var.features
write.table(highly_variable_genes,opt$out_high_var_genes,quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)