library(SeuratObject)
library(Seurat)
library(Matrix)
library(optparse)
library(dplyr)
library(tidyr)
option.list <- list(
    make_option("--SeuratObject", type="character"),
    make_option("--min_gene", type="numeric"),
    make_option("--min_umi", type="numeric"),
    make_option("--outdir", type="character"),
    make_option("--scratch", type="character")
)
opt <- parse_args(OptionParser(option_list=option.list))

#########################################################
greplany <- function(patterns, v) {
  match <- rep(FALSE, length(v))
  for (pattern in patterns) {
    match <- match | grepl(pattern, v)
  }
  return(match)
}
#########################################################
s <- readRDS(opt$SeuratObject)
tokeep <- which(!(greplany(c("^LINC","-AS","^MIR","RNU","^LOC"), rownames(s))))
df <- as.data.frame(rownames(s)[tokeep])
colnames(df) <- "Gene"
df$Gene <- as.character(df$Gene)
#########################################################
# keep genes with unique ensemblIDs and entrezIDs
library("AnnotationDbi")
library("org.Hs.eg.db")
df$ensemblID <- mapIds(org.Hs.eg.db, keys=df$Gene, column="ENSEMBL", keytype="SYMBOL", multiVals="first")
df <- df %>% filter(!is.na(ensemblID))%>% distinct(ensemblID, .keep_all=TRUE)
df$entrezID <- mapIds(org.Hs.eg.db, keys=df$ensemblID, column="ENTREZID", keytype="ENSEMBL", multiVals="first")
df <- df %>% filter(!is.na(entrezID))%>% distinct(entrezID , .keep_all=TRUE)
#########################################################
s.subset <- s[rownames(s) %in% df$Gene,]
rm(s)
s.subset <- subset(s.subset, subset= nCount_RNA > opt$min_gene & nFeature_RNA > opt$min_umi)
umi_mtx <- GetAssayData(object = s.subset, slot = "counts", assay="RNA")
#mtx <- Matrix(as.matrix(umi_mtx), sparse=T)
mtx <- umi_mtx
print(mtx[1:10,1:10])
print(colnames(mtx)[1:10])
print(rownames(mtx)[1:10])
writeMM(obj=mtx,file= file.path(opt$scratch, "matrix.mtx"))
write(x = rownames(mtx), file = file.path(opt$scratch,  "features.tsv"))
write(x = colnames(mtx), file = file.path(opt$scratch,  "barcodes.tsv"))
write.table(df, file.path(opt$scratch,  "gene_df.tsv"), sep="\t", quote=F, row.names=F)
