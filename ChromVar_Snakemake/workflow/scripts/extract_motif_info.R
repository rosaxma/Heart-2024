suppressPackageStartupMessages(library(optparse))
option.list <- list(
    make_option("--GenomeAssembly", type="character"),
    make_option("--ArchR", type="character"),
    make_option("--Outdir", type="character"),
    make_option("--SampleName", type="character"),
    make_option("--Threads", type="numeric")
)

opt <- parse_args(OptionParser(option_list=option.list))
suppressPackageStartupMessages(library(ArchR))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(Seurat))
set.seed(42)

addArchRGenome(opt$GenomeAssembly)
addArchRThreads(threads=opt$Threads)

#######################################################################
proj <- loadArchRProject(opt$ArchR) 
#####################################################################
motif_mtx <- getMatrixFromProject(ArchRProj=proj, useMatrix="MotifMatrix")
motif_z_mtx <- assays(motif_mtx)$z
motif_list <- sub("_[^_]+$", "",rownames(motif_z_mtx))
cbc_list <- colnames(motif_z_mtx)
print(motif_list[1:10])
print(cbc_list[1:10])
print(motif_z_mtx[1:5, 1:5])
write.table(motif_list, file.path(opt$Outdir, "motif_list.tsv"), row.names=F, col.names=F, sep="\t", quote=F)
write.table(cbc_list, file.path(opt$Outdir, "atac_cbc_list.tsv"), row.names=F, col.names=F, sep="\t", quote=F)
write.table(motif_z_mtx, file.path(opt$Outdir, "motif_deviation_z_score.tsv"),row.names=T, col.names=NA, sep="\t",quote=F)