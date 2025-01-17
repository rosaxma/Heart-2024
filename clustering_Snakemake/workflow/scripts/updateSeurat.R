suppressPackageStartupMessages(library(optparse))
library(SeuratObject)
library(Seurat)
library(tidyverse)
option.list <- list(
    make_option("--inputSeuratObject", type="character", default=""),
    make_option("--transformed_mtx", type="character", default=""),
    make_option("--outdir", type="character", default=""),
    make_option("--variants_to_regress", type="character", default=NA),
    make_option("--sex_genes", type="character", default=NA)
)
opt <- parse_args(OptionParser(option_list=option.list))
DATADIR <- opt$transformed_mtx
s <- readRDS(opt$inputSeuratObject)
## Read 10X Matrix
transformed_gex_mtx <- Read10X(DATADIR,gene.column=1)
print("debug")
print(str(transformed_gex_mtx))
print(str(s))
s <- SetAssayData(s, assay="RNA", slot="data", transformed_gex_mtx)

if (!is.na(opt$variants_to_regress)) {
    regress_var = opt$variants_to_regress%>% strsplit(split=",") %>% unlist()
    s$CC.difference = s$S.Score - s$G2M.Score
    regress_var = c(regress_var, "CC.difference")
    print(regress_var)
    s <- FindVariableFeatures(object = s)
    s <- ScaleData(s, vars.to.regress=regress_var)
} else {
    s <- FindVariableFeatures(object = s, nfeatures=2000)
    HVG=VariableFeatures(s)
    mito_genes=HVG[grep("^MT-.*", HVG)]
    ribo_genes=HVG[grep("^RPS.*", HVG)]
    ribo_genes=c(ribo_genes, HVG[grep("^RPL.*", HVG)])
    cell_cycle_genes=c(cc.genes$s.genes,cc.genes$g2m.gene)
        if (!is.na(opt$sex_genes)){
        sex_genes_df=read.table(opt$sex_genes, header=T, sep="\t", stringsAsFactors=F)
        print(dim(sex_genes_df))
        sex_genes=intersect(sex_genes_df$x,HVG)
        genes2remove=c(mito_genes, ribo_genes, cell_cycle_genes,sex_genes)
    } else {
        genes2remove=c(mito_genes, ribo_genes, cell_cycle_genes)
    }
    write.table(genes2remove, file.path(opt$outdir, "genes_removed_from_HVG.tsv"), quote=F, row.names=F)
    print(length(HVG))
    print(length(genes2remove))
    HVG_keep=setdiff(HVG, genes2remove)
    print(length(HVG_keep))
    VariableFeatures(s)=HVG_keep
    s <- ScaleData(s)
}
s <- RunPCA(s,features = VariableFeatures(object = s), verbose = FALSE)
pdf(file.path(opt$outdir,"elbow_plot.pdf"))
p <- ElbowPlot(s)
print(p)
dev.off()
s <- FindNeighbors(s, dims = 1:15)
s=DietSeurat(s,scale.data = FALSE, dimreducs="pca",graphs=c("RNA_nn","RNA_snn"))
s <- RunUMAP(s, dims=1:15)
filename=basename(opt$inputSeuratObject)
new_filename=gsub(".RDS", ".pflogpf.RDS", filename)
saveRDS(s,file.path(opt$outdir, new_filename))

