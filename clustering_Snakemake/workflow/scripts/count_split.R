.libPaths("/oak/stanford/groups/engreitz/Users/rosaxma/software/R_scSHC")
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(countsplit))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(sctransform))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(tidyverse))

option.list <- list(
	make_option("--inputSeuratObject", type="character", default=""),
	make_option("--threads", type="numeric", default=""),
	make_option("--Sample", type="character", default=""),
	make_option("--OutDir", type="character", default="")
)
opt <- parse_args(OptionParser(option_list=option.list))
library(future)
options(future.globals.maxSize = 8000 * 1024^2)
plan("multicore", workers = opt$threads)
set.seed(42)
s <- readRDS(opt$inputSeuratObject)
# gene by cell
counts=s[["RNA"]]$counts
# need: cell by gene
overdisps.est <- sctransform::vst(counts,n_genes=NULL)$model_pars[,1]
counts=counts[names(overdisps.est),]
#############################################################################################
save_mtx <- function(counts, OutDir){
	mtx <- Matrix(counts, sparse=T)
	writeMM(obj=mtx,file= file.path(OutDir,"matrix.mtx"))
	write(x = rownames(mtx), file = file.path(OutDir,"features.tsv"))
	write(x = colnames(mtx), file = file.path(OutDir,"barcodes.tsv"))
	}
#############################################################################################
if (identical(names(overdisps.est), rownames(counts))){
	split <- countsplit(t(counts),overdisps=overdisps.est)
	Xtrain <- split[[1]]
	Xtest <- split[[2]]
	print(str(split))
	metadata=s@meta.data
	train_metadata=metadata[colnames(t(Xtrain)),]
	test_metadata=metadata[colnames(t(Xtest)),]
	options(Seurat.object.assay.version = "v3")
	train_seurat <- CreateSeuratObject(counts=t(Xtrain),meta.data=train_metadata)
	print(str(train_seurat))
	test_seurat <- CreateSeuratObject(counts=t(Xtest),meta.data=test_metadata)
	print(str(test_seurat))
	saveRDS(train_seurat,file.path(opt$OutDir, "train_matrix", "train.RDS"))
	saveRDS(test_seurat, file.path(opt$OutDir, "test_matrix", "test.RDS"))
	##################################################
	save_mtx(t(Xtrain),file.path(opt$OutDir, "train_matrix"))
	save_mtx(t(Xtest),file.path(opt$OutDir,"test_matrix"))	
} else {
	print("Gene counts don't match.")
}
