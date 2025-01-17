suppressPackageStartupMessages(library(optparse))
option.list <- list(
    make_option("--GenomeAssembly", type="character", default="hg38"),
    make_option("--WorkDir", type="character"),
    make_option("--ArrowFile", type="character"),
    make_option("--Threads", type="numeric"), 
    make_option("--OutDir", type="character"),
    make_option("--SampleName", type="character"), 
    make_option("--clustering_info", type="character"),
    make_option("--clustering_col", type="character")
)

opt <- parse_args(OptionParser(option_list=option.list))
setwd(opt$WorkDir)
suppressPackageStartupMessages(library(ArchR))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggpubr))
set.seed(42)

addArchRGenome(opt$GenomeAssembly)
addArchRThreads(threads=opt$Threads)
#######################################################################
arrowfile=opt$ArrowFile
proj <- ArchRProject(
    ArrowFiles = arrowfile,
    outputDirectory = file.path(opt$WorkDir, "ArchRProject"),
    copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)
#######################################################################
clustering_table <- read.table(opt$clustering_info, header=T, sep="\t", stringsAsFactors=F) %>% select(sample_orig, ATAC_barcode, !!sym(opt$clustering_col), source) %>% mutate(ATAC_cbc=paste0(sample_orig, "_", ATAC_barcode)) %>% select(ATAC_cbc, !!sym(opt$clustering_col), source)
print(head(clustering_table))

archr_df <- getCellColData(proj, select = c("Sample","nFrags", "TSSEnrichment"))
rownames(archr_df) <- proj$cellNames
archr_df <- archr_df %>% as.data.frame() %>% rownames_to_column(var="cbc") %>% separate_wider_delim(cbc, delim="#", names=c("cell_id", "ATAC_cbc")) %>% select(-cell_id) %>% inner_join(clustering_table, by=c("ATAC_cbc")) %>% as.data.frame()

proj$Celltype=archr_df %>% pull(!!sym(opt$clustering_col))

proj$Source=archr_df$source
#######################################################################
qc_table <- archr_df %>% select(-source) %>% as.data.frame() %>% group_by(!!sym(opt$clustering_col)) %>% summarise(mean_ATAC_frags=mean(nFrags), median_ATAC_frags=median(nFrags), mean_ATAC_TSS_enrichment=mean(TSSEnrichment), median_ATAC_TSS_enrichment=median(TSSEnrichment), atac_cell_number=n()) %>% as.data.frame() %>% ungroup()

write.table(qc_table,file.path(opt$OutDir, "QC_table.tsv"), quote=FALSE, sep="\t", row.names=FALSE)
#######################################################################
proj  <- addIterativeLSI(
    ArchRProj = proj ,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30
)
proj <- addClusters(
    input = proj,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.1
)

proj <- addUMAP(
    ArchRProj =proj, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
)
pdf(file.path(opt$OutDir, paste0(opt$SampleName, "_UMAP_Celltype_source.pdf")), width=5, height=5)
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Celltype", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Source", embedding = "UMAP")
print(p1)
print(p2)
dev.off()
#######################################################################
proj <- saveArchRProject(ArchRProj = proj,load = TRUE)
#######################################################################
proj <- addGroupCoverages(ArchRProj=proj, groupBy="Celltype", useLabels=FALSE)
pathToMacs2 <- findMacs2()
print(pathToMacs2)
proj <- addReproduciblePeakSet(
    ArchRProj = proj, 
    groupBy = "Celltype", 
    pathToMacs2 = pathToMacs2
)
proj <- addPeakMatrix(proj)
proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif")
proj <- addBgdPeaks(proj)
proj <- addDeviationsMatrix(
  ArchRProj = proj,
  peakAnnotation = "Motif",
  force = TRUE
)
saveArchRProject(ArchRProj = proj,load = FALSE)
#######################################################################
plotVarDev <- getVarDeviations(proj, name = "MotifMatrix", plot = FALSE)
#dev_df=annoGR2DF(plotVarDev)
dev_df=plotVarDev
write.table(dev_df, file.path(opt$OutDir, "chromvar_output.tsv"), quote=FALSE, sep="\t", row.names=FALSE)
plotVarDev <- getVarDeviations(proj, name = "MotifMatrix", plot = TRUE)
plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

