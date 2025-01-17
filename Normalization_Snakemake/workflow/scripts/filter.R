suppressPackageStartupMessages(library(optparse))

option.list <- list(
    make_option("--inputSeuratObject", type="character", default=""), 
    make_option("--min_nUMI",type="numeric", default=0), 
    make_option("--max_nUMI",type="numeric", default=99999999999),
    make_option("--min_nGene", type="numeric", default=0),
    make_option("--max_nGene", type="numeric", default=99999999999),
    make_option("--ratio_mt", type="numeric", default=1),
    make_option("--outdir", type="character"),
    make_option("--sample", type="character"),
    make_option("--atacBarcodes", type="character", default=NA)
)
opt <- parse_args(OptionParser(option_list=option.list))
library(SeuratObject)
library(Seurat)
library(Matrix)

s <- readRDS(opt$inputSeuratObject)
cellnumber_before_filtering=length(colnames(s))
if (!(is.na(opt$atacBarcodes))) {
    atac_barcode_table <- read.table(opt$atacBarcode, header=T, sep="\t", stringsAsFactors=F)
    s=subset(s, cells=atac_barcode_table$GEX_Barcodes)
    cellnumber_before_filtering_after_ATAC_filtering=length(colnames(s))
}

print(colnames(s@meta.data))
print(opt$ratio_mt)
s <- subset(s, subset = nUMI > opt$min_nUMI & nUMI < opt$max_nUMI & nGene > opt$min_nGene & nGene < opt$max_nGene & ratio.mt < opt$ratio_mt)
cellnumber_after_filtering=length(colnames(s))
### assigning cell cycle status
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
s <- CellCycleScoring(s, s.features = s.genes, g2m.features = g2m.genes)

saveRDS(s, file.path(dirname(opt$outdir),"QC", gsub(".RDS", ".filtered.RDS", basename(opt$inputSeuratObject))))

if (is.na(opt$atacBarcodes)){
    df <- data.frame(cellnumber_before_filtering,cellnumber_after_filtering)
    rownames(df)=opt$sample
    write.table(df, file.path(opt$outdir,"cell_count_before_vs_after_filtering.tsv"), sep="\t", row.names=T, quote=FALSE)
} else {
    df <- data.frame(cellnumber_before_filtering,cellnumber_before_filtering_after_ATAC_filtering,cellnumber_after_filtering)
    rownames(df)=opt$sample
    write.table(df, file.path(opt$outdir,"cell_count_before_vs_after_filtering.tsv"), sep="\t", row.names=T, quote=FALSE)
}


umi_mtx <- GetAssayData(object = s, slot = "counts", assay="RNA")
mtx <- Matrix(as.matrix(umi_mtx), sparse=T)
writeMM(obj=mtx,file= file.path(opt$outdir, "matrix.mtx"))
write(x = rownames(mtx), file = file.path(opt$outdir, "features.tsv"))
write(x = colnames(mtx), file = file.path(opt$outdir, "barcodes.tsv"))
