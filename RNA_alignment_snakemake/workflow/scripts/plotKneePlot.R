suppressPackageStartupMessages(library(optparse))
option.list <- list(
    make_option("--outdir", type="character"), 
    make_option("--matrixdir",type="character"),
    make_option("--sample", type="character")
)
opt <- parse_args(OptionParser(option_list=option.list))
suppressPackageStartupMessages(library(SeuratObject))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggridges))
suppressPackageStartupMessages(library(ggrastr))

cell_matrix <- Read10X(opt$matrixdir)
s <- CreateSeuratObject(cell_matrix)
df <- rownames(s@meta.data) %>% as.data.frame()
colnames(df) <- "CellBarcodes"
df$UMI <- s@meta.data$nCount_RNA
df <- df %>% mutate(rank=rank(-UMI))
pdf(paste0(opt$outdir,"/",opt$sample,"Kneeplot.pdf"))
p <- ggplot(df, aes(x=rank, y=UMI)) + rasterize(geom_point(), dpi=300) + scale_x_log10(limits = c(1e0, 1e6),n.breaks=14) + scale_y_log10(limits = c(1e0, 1e6),n.breaks=7) + ggtitle(opt$sample)
print(p)
dev.off()
