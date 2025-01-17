suppressPackageStartupMessages(library(optparse))
option.list <- list(
    make_option("--inputSeuratObject", type="character"),
    make_option("--doublet_to_exclude",type="character"),
    make_option("--outdir",type="character"),
    make_option("--sample", type="character")
    )

opt <- parse_args(OptionParser(option_list=option.list))
suppressPackageStartupMessages(library(SeuratObject))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
input_doublets <- read.table(opt$doublet_to_exclude, header=TRUE, sep = ',', stringsAsFactors=FALSE)$Barcodes
s <- readRDS(opt$inputSeuratObject)
##############################################################################
barcodes <- as.data.frame(colnames(s))
write.table(barcodes,paste0(opt$outdir, "/", opt$sample,"_cbc_before_doublet_removal.tsv"), quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
##############################################################################
avail_cells=colnames(s)
doublets=intersect(input_doublets,avail_cells)
s <- subset(s, cells=doublets, invert=TRUE)
metadata <- s@meta.data
metadata$orig.ident=opt$sample
s@meta.data=metadata
print(colnames(s@meta.data))
saveRDS(s, paste0(opt$outdir,"/", opt$sample, "_SeuratObject_doublet_removed.RDS"))
##############################################################################
metadata <- s@meta.data
df <- metadata %>% summarise(mean_UMI = mean(nUMI), median_UMI=median(nUMI),mean_gene_counts=mean(nGene), median_gene_counts=median(nGene), mean_ratio_mt=mean(ratio.mt), mean_nUMI_non_MT=mean(nUMI_non_MT), median_nUMI_non_MT=median(nUMI_non_MT), mean_ratio_mt=mean(ratio.mt), gex_cell_number=n()) %>% as.data.frame() %>% ungroup()
write.table(df,paste0(opt$outdir, "/", opt$sample,"_gex_data_QC_post_doublet_removal.tsv"), quote=FALSE, sep="\t", row.names=FALSE)