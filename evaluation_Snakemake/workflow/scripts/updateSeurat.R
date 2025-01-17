suppressPackageStartupMessages(library(optparse))
library(SeuratObject)
library(Seurat)
option.list <- list(
    make_option("--inputSeuratObject", type="character", default=""),
    make_option("--transformed_mtx", type="character", default="")
)
opt <- parse_args(OptionParser(option_list=option.list))
DATADIR <- opt$transformed_mtx
s <- readRDS(opt$inputSeuratObject)
## Read 10X Matrix
transformed_gex_mtx <- Read10X(DATADIR,gene.column=1)
print("debug")
s <- SetAssayData(s, assay="RNA", slot="data", transformed_gex_mtx)
saveRDS(s, gsub(".RDS", ".pflogpf.RDS", opt$inputSeuratObject))
