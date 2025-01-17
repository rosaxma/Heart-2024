suppressPackageStartupMessages(library(optparse))
option.list <- list(
    make_option("--GenomeAssembly", type="character", default="hg38"),
    make_option("--WorkDir", type="character"),
    make_option("--ArrowFile", type="character"),
    make_option("--Threads", type="numeric"), 
    make_option("--OutDir", type="character"),
    make_option("--SampleName", type="character"), 
    make_option("--markdupFile", type="character")
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
doubScores <- addDoubletScores(
    input = arrowfile, 
    k=10, 
    knnMethod = "UMAP", 
    LSIMethod = 1
)
proj <- ArchRProject(
    ArrowFiles = arrowfile,
    outputDirectory = file.path(opt$WorkDir, "ArchRProject_basic"),
    copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)

#######################################################################
atac_bc_prefiltering <- proj$cellNames
proj_filtered <- filterDoublets(proj)
atac_bc_postfiltering <- proj_filtered$cellNames

isDoublets <- unlist(lapply(atac_bc_prefiltering, function(x) !(x %in% atac_bc_postfiltering)))
names(isDoublets) = proj$cellNames
proj$isDoublets <- isDoublets
#######################################################################
df <- getCellColData(proj, select = c("Sample","nFrags", "TSSEnrichment"))
rownames(df) <- proj$cellNames
qc_table <- df %>% as.data.frame() %>% group_by(Sample) %>% summarise(mean_ATAC_frags=mean(nFrags), median_ATAC_frags=median(nFrags), mean_ATAC_TSS_enrichment=mean(TSSEnrichment), median_ATAC_TSS_enrichment=median(TSSEnrichment), atac_cell_number=n()) %>% as.data.frame() %>% ungroup()

if (!(opt$markdupFile == "" )){
    temp <- readLines(opt$markdupFile)
    start=grep('ESTIMATED_LIBRARY_SIZE', temp)
    end=start+1
    markdup_file=read.table(text=paste0(trimws(temp[start:end]), collapse="\n"),sep="\t", header=T)
    qc_table <- qc_table %>% mutate(percent_duplication=as.numeric(markdup_file[1,"PERCENT_DUPLICATION"])) %>% mutate(estimated_library_size=as.numeric(markdup_file[1,"ESTIMATED_LIBRARY_SIZE"])) %>% mutate(read_pairs_examined=as.numeric(markdup_file[1,"READ_PAIRS_EXAMINED"]))
}

write.table(qc_table,file.path(opt$OutDir, "ATAC_QC_table.tsv"), quote=FALSE, sep="\t", row.names=FALSE)


#######################################################################
#save cell nfrag and TSS info
df_cell_info <- df %>% as.data.frame() 
df_cell_info$barcode <- word(proj$cellNames,2,sep = "#")
write.table(df_cell_info,file.path(opt$OutDir, "ATAC_cell_info.csv"), quote=FALSE, sep=",", row.names=FALSE)
#######################################################################
df_doublets <- getCellColData(proj, select = c("isDoublets"))
df_doublets<- df_doublets %>% as.data.frame()
df_doublets$barcodes <- word(proj$cellNames,2,sep = "#")
df_doublets <- df_doublets %>% filter(isDoublets)
write.table(df_doublets [,c("barcodes")],file.path(opt$OutDir,"ArchR_doublet_barcodes_in_ATAC.csv"), quote=FALSE,sep=",", col.names= FALSE, row.names=FALSE)

#######################################################################
atac_cbc<-as.data.frame(proj$cellNames)
atac_cbc$barcode <- word(proj$cellNames,2,sep = "#")
atac_cbc$is__cell_barcode <- 1

write.table(atac_cbc[,c("barcode","is__cell_barcode")],file.path(opt$OutDir,"filtered_barcodes_in_ATAC.csv"), quote=FALSE,sep=",", row.names=FALSE)
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
    resolution = 0.8
)

proj <- addUMAP(
    ArchRProj =proj, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
)
pdf(file.path(opt$OutDir, paste0(opt$SampleName, "_rawUMAP.pdf")), width=5, height=5)
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
print(p1)
dev.off()
#######################################################################
#plot
p1 <- plotGroups(
    ArchRProj = proj, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "ridges"
   )

p2 <- plotGroups(
    ArchRProj = proj, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
)

p3 <- plotGroups(
    ArchRProj = proj, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "log10(nFrags)",
    plotAs = "ridges"
   )

p4 <- plotGroups(
    ArchRProj = proj, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "log10(nFrags)",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )

p5 <- plotFragmentSizes(ArchRProj = proj)

p6 <- plotTSSEnrichment(ArchRProj = proj)

#######################################################################
p7 <- ggplot(df%>% as.data.frame(), aes(log10(nFrags))) + stat_ecdf(geom = "point") + ggtitle("Cumulative distribution of nFrags") + theme_classic() + theme(axis.text = element_text(size = 10)) 
p8 <- ggplot(df%>% as.data.frame(), aes(TSSEnrichment)) + stat_ecdf(geom = "point") + ggtitle("Cumulative distribution of TSSEnrichment") + theme_classic() + theme(axis.text = element_text(size = 10)) 
#######################################################################
pdf(file.path(opt$OutDir, paste0(opt$SampleName, "_ATAC_QC.pdf")), width=10, height=10)
figure <- ggarrange(ggarrange(p1,p2, ncol=2, labels=rep("TSS enrichment", 2),widths = c(1, 1)), ggarrange(p3, p4, ncol = 2, labels=rep("Unique nuclear fragment", 2), widths=c(1,1)), ggarrange(p5, p6, ncol=2, labels=c("Fragment Size Distribution", "TSS enrichment profiles"), widths=c(1,1)), ggarrange(p7, p8, ncol=2, widths=c(1,1)),nrow = 4)
print(figure)
dev.off()

