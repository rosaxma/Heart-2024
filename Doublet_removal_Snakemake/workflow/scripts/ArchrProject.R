suppressPackageStartupMessages(library(optparse))
option.list <- list(
    make_option("--GenomeAssembly", type="character", default="hg38"),
    make_option("--barcode_sheet", type="character"),
    make_option("--WorkDir", type="character"),
    make_option("--ArrowFile", type="character"),
    make_option("--Threads", type="numeric"), 
    make_option("--OutDir", type="character"),
    make_option("--SampleName", type="character")
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
input_barcodes <- read.table(opt$barcode_sheet, header=TRUE, sep = '\t', stringsAsFactors=FALSE, comment.char="") 
cells_to_keep= input_barcodes%>% pull(ATAC_Barcodes_achr)
arrowfile=opt$ArrowFile
proj <- ArchRProject(
    ArrowFiles = arrowfile,
    outputDirectory = file.path(opt$WorkDir, "ArchRProject_prefilter"),
    copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)

proj_filtered <- subsetArchRProject(
  ArchRProj = proj,
  cells = cells_to_keep,
  outputDirectory = file.path(opt$WorkDir, "ArchRProject_multiome"),
  dropCells = TRUE,
  logFile = NULL,
  force=TRUE
)

#plot UMAP
proj_filtered  <- addIterativeLSI(
    ArchRProj = proj_filtered ,
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

proj_filtered <- addClusters(
    input = proj_filtered,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.8
)

proj_filtered <- addUMAP(
    ArchRProj =proj_filtered, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
)
pdf(file.path(opt$OutDir, paste0(opt$SampleName, "_UMAP.pdf")))
p1 <- plotEmbedding(ArchRProj = proj_filtered, colorBy = "cellColData", name = "Clusters", embedding = "UMAP") + ggtitle(opt$SampleName)
print(p1)
dev.off()
saveArchRProject(ArchRProj=proj_filtered, outputDirectory=file.path(opt$WorkDir, "ArchRProject_multiome"), load=FALSE)
