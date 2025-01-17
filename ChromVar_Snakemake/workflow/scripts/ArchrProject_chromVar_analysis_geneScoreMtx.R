suppressPackageStartupMessages(library(optparse))
option.list <- list(
    make_option("--GenomeAssembly", type="character"),
    make_option("--GEXArchR", type="character"),
    make_option("--Outdir", type="character"),
    make_option("--SampleName", type="character"),
    make_option("--Threads", type="numeric")
)

opt <- parse_args(OptionParser(option_list=option.list))
suppressPackageStartupMessages(library(ArchR))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(cowplot))
set.seed(42)

addArchRGenome(opt$GenomeAssembly)
addArchRThreads(threads=opt$Threads)

#######################################################################
proj <- loadArchRProject(opt$GEXArchR) 
#######################################################################
seGroupMotif <- getGroupSE(ArchRProj = proj, useMatrix = "MotifMatrix", groupBy = "Celltype")
seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]
rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
    rowMaxs(assay(seZ) - assay(seZ)[,x])
    }) %>% Reduce("cbind", .) %>% rowMaxs
#######################################################################
corGEM_MM <- correlateMatrices(
    ArchRProj = proj,
    useMatrix1 = "GeneScoreMatrix",
    useMatrix2 = "MotifMatrix",
    reducedDims = "IterativeLSI"
)
corGEM_MM$maxDelta <- rowData(seZ)[match(corGEM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]

corGEM_MM <- corGEM_MM[order(abs(corGEM_MM$cor), decreasing = TRUE), ]
corGEM_MM <- corGEM_MM[which(!duplicated(gsub("\\-.*","",corGEM_MM[,"MotifMatrix_name"]))), ]
corGEM_MM$TFRegulator <- "NO"
corGEM_MM$TFRegulator[which(corGEM_MM$cor > 0.5 & corGEM_MM$padj < 0.05 & corGEM_MM$maxDelta > quantile(corGEM_MM$maxDelta, 0.5))] <- "YES"
write.table(corGEM_MM, file.path(opt$Outdir, paste0(opt$SampleName,"_motif_TF_expression_correlation.tsv")), sep="\t", row.names=F, quote=F)
motifs=sort(corGEM_MM[corGEM_MM$TFRegulator=="YES",1])
print(motifs)
#######################################################################
ncol=3
if ((length(motifs) %% ncol) != 0){
        nrow=length(motifs)/ncol + 1
} else {
        nrow=length(motifs)/ncol
}
#######################################################################
markerMotifs <- getFeatures(proj, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs <- grep("z:", markerMotifs, value = TRUE)

p <- plotGroups(ArchRProj = proj, 
  groupBy = "Celltype", 
  colorBy = "MotifMatrix", 
  name = markerMotifs,
  imputeWeights = getImputeWeights(proj)
)

p2 <- lapply(seq_along(p), function(x){
  if(x != 1){
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6) +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
    theme(
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
    ) + ylab("")
  }else{
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6) +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
    theme(
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
    ) + ylab("")
  }
})

pdf(file.path(opt$Outdir, paste0(opt$SampleName,"_Marker_TF_expression_in_cluster.pdf")), width=3*length(motifs), height=length(unique(proj$Celltype))*1.5)
print(do.call(cowplot::plot_grid, c(list(nrow = 1, rel_widths = c(3, rep(2, length(p2) - 1))),p2)))
dev.off()


p <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "MotifMatrix", 
    name = sort(markerMotifs), 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(proj)
)
p2 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
    )
})

pdf(file.path(opt$Outdir, paste0(opt$SampleName,"_Motif_zscore_on_UMAP.pdf")), width=ncol*2, height=nrow*2)
print(do.call(cowplot::plot_grid, c(list(ncol = ncol),p2)))
dev.off()

