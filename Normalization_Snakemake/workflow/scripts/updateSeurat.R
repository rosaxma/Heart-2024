suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(SeuratObject))
suppressPackageStartupMessages(library(Seurat))
option.list <- list(
    make_option("--inputSeuratObject", type="character"),
    make_option("--transformed_mtx", type="character")
)
opt <- parse_args(OptionParser(option_list=option.list))
transformed_mtx <- opt$transformed_mtx
s <- readRDS(opt$inputSeuratObject)
## Read 10X Matrix
transformed_gex_mtx <- Read10X(transformed_mtx,gene.column=1)

s <- SetAssayData(s, assay="RNA", slot="data", transformed_gex_mtx)
saveRDS(s, gsub(".RDS", ".pflogpf.RDS", opt$inputSeuratObject))
