suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(SeuratObject))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(tidyverse))
option.list <- list(
    make_option("--inputSeuratObject", type="character"),
    make_option("--transformed_mtx", type="character"),
    make_option("--outdir", type="character"),
    make_option("--res_range", type="character", default="0,1"),
    make_option("--increment", type="numeric", default="0.1"),
    make_option("--variants_to_regress", type="character", default=NA),
    make_option("--sex_genes", type="character", default=NA)
)
opt <- parse_args(OptionParser(option_list=option.list))
transformed_mtx <- opt$transformed_mtx
s <- readRDS(opt$inputSeuratObject)
## Read 10X Matrix
transformed_gex_mtx <- Read10X(transformed_mtx,gene.column=1)
s <- SetAssayData(s, assay="RNA", slot="data", new.data =transformed_gex_mtx)
saveRDS(s, gsub(".RDS", ".pflogpf.RDS", opt$inputSeuratObject))
########################################################################
meta <- s@meta.data
meta <- meta %>% select(-starts_with("RNA_snn")) %>% select(-starts_with("res_"))
s@meta.data <- meta
########################################################################
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
        print(opt$sex_genes)
        sex_genes_df=read.table(opt$sex_genes, header=T, sep="\t", stringsAsFactors=F)
        print(head(sex_genes_df))
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
########################################################################
start=as.numeric(unlist(strsplit(opt$res_range, split=","))[1])
end=as.numeric(unlist(strsplit(opt$res_range, split=","))[2])
res_list <- seq(start, end, by=as.numeric(opt$increment))
cluster_count <- c()
for (res in res_list){
    s <- FindClusters(s, resolution=as.numeric(res))
    metadata <- s@meta.data
    rownames <- rownames(metadata)
    colname <- paste0("res_", as.character(res))
    count <- length(unique(metadata$seurat_clusters))
    metadata <- metadata %>% rename(!!colname :="seurat_clusters")
    rownames(metadata) <- rownames
    s@meta.data <- metadata
    cluster_count <- c(cluster_count,count)
}
s <- RunUMAP(s, dims=1:15)
filename=basename(opt$inputSeuratObject)
new_filename=gsub(".RDS", ".cluster.pick.res.RDS", filename)
saveRDS(s,file.path(opt$outdir, new_filename))
#############################################################################
res_count=length(res_list)
pdf(file.path(opt$outdir, "clustree.pdf"), height=20, width=30)
library(clustree)
print(str(s))
if (res_count>1) {
	clustree(s, prefix = "res_", layout = "sugiyama")
}
dev.off()
df <- data.frame(res_list,cluster_count)
write.table(df, file.path(opt$outdir,"cluster_count_res.tsv"), quote=FALSE, sep = "\t", row.names = FALSE)
#############################################################################
metadata_clustering <- s@meta.data %>% select(starts_with("res_"))
write.table(metadata_clustering, file.path(opt$outdir,"subclustering_results.tsv"), quote=FALSE, sep = "\t", row.names = TRUE)
