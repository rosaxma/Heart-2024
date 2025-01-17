options(future.globals.maxSize = 4000 * 1024^5)
suppressPackageStartupMessages(library(optparse))

option.list <- list(
    make_option("--outdir", type="character", default="", help="Output directory"), 
    make_option("--allSeuratObjects", type="character", default="", help="Path to the Seurat Objects separated by comma"),
    make_option("--effectiveSheet", type="character")
    )

opt <- parse_args(OptionParser(option_list=option.list))
library(SeuratObject)
library(Seurat)
library(Matrix)
suppressPackageStartupMessages(library(tidyverse))
set.seed(42)
SeuratObjects=opt$allSeuratObjects %>% strsplit(split=",") %>% unlist()
print(SeuratObjects)
cell_count_list <- vector(mode = "list", length = length(SeuratObjects))
if (length(cell_count_list) > 1){
	first_object=readRDS(SeuratObjects[1])
	other_objects=lapply(SeuratObjects[c(-1)],readRDS)
	first_object_id=as.character(unique(first_object$orig.ident)[1])
	print(first_object_id)
	other_objects_id <- vector()
	for (object in other_objects){
    		id <- as.character(unique(object$orig.ident)[1])
    		print(id)
    		other_objects_id <- append(other_objects_id, id)
	}
	test_barcode=colnames(first_object)[[1]]
	if (grepl("_", test_barcode, fixed=TRUE)){
		print(test_barcode)
		print("check")
		combined_object=merge(first_object, other_objects, project="combined_objects")
	} else {
		print(test_barcode)
		print("check2")	
		combined_object=merge(first_object, other_objects, add.cell.ids=c(first_object_id,other_objects_id), project="combined_objects")
	}
} else {
	combined_object=readRDS(SeuratObjects[1])
}
effective_sheet <- read.table(opt$effectiveSheet, header=T, sep="\t", stringsAsFactors=F) %>% mutate(cbc=paste0(orig.ident,"_",cbc))
avail_cbc=intersect(effective_sheet$cbc,colnames(combined_object))
combined_object=subset(combined_object,cells=avail_cbc)
saveRDS(combined_object, file.path(opt$outdir,"combined_SeuratObject_doublet_removed.RDS"))

umi_mtx <- GetAssayData(object = combined_object, slot = "counts", assay="RNA")
mtx <- Matrix(umi_mtx, sparse=T)
writeMM(obj=mtx,file= file.path(opt$outdir, "combined_matrix","matrix.mtx"))
write(x = rownames(mtx), file = file.path(opt$outdir, "combined_matrix","features.tsv"))
write(x = colnames(mtx), file = file.path(opt$outdir, "combined_matrix","barcodes.tsv"))


