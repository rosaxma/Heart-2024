suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(Matrix))
option.list <- list(
	make_option("--inputSeuratObject", type="character", default=""),
	make_option("--SampleInfoTable", type="character", default=""),
	make_option("--Sample", type="character", default=""),
	make_option("--OutDir", type="character", default=""),
	make_option("--iPSC", type="logical", default=""),
	make_option("--SVD", type="logical", default=""),
	make_option("--effectiveSheet", type="character", default=""),
	make_option("--minCell", type="numeric", default="")
)
opt <- parse_args(OptionParser(option_list=option.list))
s <- readRDS(opt$inputSeuratObject)
table <- read.table(opt$SampleInfoTable, sep="\t", header=T, stringsAsFactors=F)
set.seed(42)

metadata=s@meta.data
print(head(metadata))
#FT030_group=c("FT030_1","FT030_2","FT030_3")
#FT074_group=c("FT074_1","FT074_2","FT074_3")

#metadata$experiment <- sub("_[^_]+$", "", metadata$orig.ident)
#metadata <- metadata %>% mutate(source=ifelse(source %in% FT030_group, "FT030", source))
#metadata <- metadata %>% mutate(source=ifelse(source %in% FT074_group, "FT074", source))
#metadata <- metadata %>% mutate(source=ifelse(source %in% c("a","b","c"), paste0(experiment,"_",source), source))

for (col in colnames(table)){
	if ( col != "source"){
		if (col %in% colnames(metadata)) {
			metadata <- metadata %>% select(-!!sym(col))
		}	
	}
}
if ("source" %in% colnames(metadata)){
	metadata <- metadata %>% left_join(table, by=c("source"))
	rownames(metadata) <- metadata$cbc
} else{
	effective_sheet <- read.table(opt$effectiveSheet, sep="\t", header=T, stringsAsFactors=F) %>% mutate(cbc=paste0(orig.ident, "_", cbc)) %>% select(cbc, source)
	metadata$cbc=rownames(metadata)
	metadata <- metadata %>% left_join(effective_sheet, by=c("cbc")) %>% left_join(table, by=c("source"))
	rownames(metadata) <- metadata$cbc
}
s@meta.data <- metadata
print(head(metadata))

Idents(s) <- "source"
if (opt$iPSC) {
	intersect_source=intersect(c("FT014", "FT021"), unique(s$source))
} else if (opt$SVD) {
	intersect_source=intersect(c("CS_iPSCs_1", "CS_iPSCs_2"), unique(s$source))
} else {
	intersect_source=intersect(c("FT014", "FT021","CS_iPSCs_1", "CS_iPSCs_2"), unique(s$source))
}
s_subset <- subset(x = s, idents=intersect_source,invert = TRUE)
counts=GetAssayData(object = s_subset, assay = "RNA", slot = "counts")
keep_genes <- Matrix::rowSums(counts)>=opt$minCell
filtered_counts <- counts[keep_genes,]

filtered_s <- CreateSeuratObject(counts=filtered_counts,meta.data=s_subset@meta.data)
saveRDS(filtered_s,file.path(opt$OutDir, paste0(opt$Sample, ".sample.info.RDS")))
save_mtx <- function(counts, outdir){
		print("here")
		print(file.path(outdir, "matrix.mtx"))
        mtx <- Matrix(counts, sparse=T)
        writeMM(obj=mtx,file= file.path(outdir, "matrix.mtx"))
        write(x = rownames(mtx), file = file.path(outdir, "features.tsv"))
        write(x = colnames(mtx), file = file.path(outdir, "barcodes.tsv"))
        }
save_mtx(filtered_counts,file.path(opt$OutDir, "matrix"))



