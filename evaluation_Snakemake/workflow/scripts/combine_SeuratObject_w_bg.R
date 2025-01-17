suppressPackageStartupMessages(library(optparse))

option.list <- list(
    make_option("--outdir", type="character", default="", help="Output directory"), 
    make_option("--allSeuratObjects", type="character", default="", help="Path to the Seurat Objects separated by comma")
    )

opt <- parse_args(OptionParser(option_list=option.list))
library(SeuratObject)
library(Seurat)
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
	combined_object=merge(first_object, other_objects, add.cell.ids=c(first_object_id,other_objects_id), project="combined_objects")
} else {
	combined_object=readRDS(SeuratObjects[1])
}
saveRDS(combined_object, file.path(opt$outdir,"combined_SeuratObject_w_bg.RDS"))
