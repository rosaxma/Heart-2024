suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(SeuratObject))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(hdf5r))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggridges))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggrastr))
mytheme <- theme(axis.text = element_text(size = 16), axis.title = element_text(size = 16), plot.title = element_text(hjust = 0.5, face="bold"))
#######################################################################
option.list <- list(
    make_option("--bg_removed_h5", type="character", default=""), 
    make_option("--OutDir",type="character", default=""),
    make_option("--STARfilteredMatrices", type="character"),
    make_option("--min.cells", type="integer", default=10), 
    make_option("--sample", type="character")
)
opt <- parse_args(OptionParser(option_list=option.list))
OUTDIR <- opt$OutDir
bg_removed_mtx <- Read10X_h5(opt$bg_removed_h5, use.names = TRUE)
print(dim(bg_removed_mtx))
bg_removed_mtx_filtered <- bg_removed_mtx[rowSums(bg_removed_mtx> 0) > opt$min.cells, ]
print(dim(bg_removed_mtx_filtered))

STAR_filtered_matrix <- Read10X(opt$STARfilteredMatrices)
STAR_filtered_s <- CreateSeuratObject(STAR_filtered_matrix)
STAR_filtered_cbc <- colnames(STAR_filtered_s)
rm(STAR_filtered_matrix, STAR_filtered_s)
s <- CreateSeuratObject(bg_removed_mtx)
s <- subset(s, cells = STAR_filtered_cbc)
s_filtered <- CreateSeuratObject(bg_removed_mtx_filtered)
s_filtered <- subset(s_filtered,cells = STAR_filtered_cbc)
print(length(STAR_filtered_cbc))
print(length(rownames(s)))
print(length(rownames(s_filtered)))
#######################################################################
modifyMetaData <- function(s, sampleid){
    s[["ratio.mt"]] <- PercentageFeatureSet(s, pattern = "^MT-")/100
    s[["ratio.ribo"]] <- PercentageFeatureSet(s, pattern = "^RPS|^RPL")/100
    metadata <- s@meta.data
    metadata$cells <- rownames(metadata)
    metadata <- metadata %>% dplyr::rename(nUMI=nCount_RNA,nGene=nFeature_RNA) %>% mutate(nUMI_non_MT=nUMI*(1-ratio.mt))
    metadata$sample_id <- sampleid
    s@meta.data <- metadata
    s$log10GenesPerUMI <- log10(s$nGene)/log10(s$nUMI)
    return(s)
}
s_filtered <- modifyMetaData(s_filtered, opt$sample)
#######################################################################
options(scipen=999)
plotQC <- function(s, min_cell){
    pdf(paste0(OUTDIR,"/", opt$sample, "_QC_min_cell_",min_cell,".pdf"), width=20, height=40)
    # number of cells in each cluster
    metadata <- s@meta.data
    cell_count <- table(metadata$sample_id) %>% as.data.frame()
    colnames(cell_count) <- (c("classification","cell.count"))

    p1 <- ggplot(cell_count, aes(x=classification, y=cell.count, fill=classification)) + geom_bar(stat="identity") + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position="none") + mytheme + ggtitle("nCells") + geom_text(aes(label=cell.count), position=position_dodge(width=0.9), vjust=-0.25, size=5)

    cell_count$gene.count <- length(rownames(s))
    p_gene_counts <- ggplot(cell_count, aes(x=classification, y=gene.count, fill=classification)) + geom_bar(stat="identity") + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),legend.position="none") + mytheme + ggtitle("nGenes") + geom_text(aes(label=gene.count), position=position_dodge(width=0.9), vjust=-0.25,size=5)

    # distribution of UMI
    p2 <- ggplot(metadata, aes(x=nUMI)) +
    geom_density(size=1)+scale_x_log10(limits = c(100, 70000), n.breaks=20) + theme_classic() +ylab("Cell density") + ggtitle("nUMI density") + mytheme + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

    # distribution of non_MT_UMI
    p3 <- ggplot(metadata, aes(x=nUMI_non_MT)) + geom_density(size=1) + scale_x_log10(limits = c(100, 70000), n.breaks=20) + theme_classic() +
    ylab("Cell density") + ggtitle("non MT nUMI density") + mytheme + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

    # distribution of gene number per cell
    p4 <-  ggplot(metadata, aes(x=nGene)) +
    geom_density(size=1) + theme_classic() +
    scale_x_log10(limits = c(100, 25000),n.breaks=20) + ggtitle("nGene distribution") + mytheme + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    
    # distribution of gene number per cell log transformed
    p5 <- ggplot(metadata, aes(x=sample_id, y=log10(nGene), fill=sample_id)) + geom_boxplot() +theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +theme(plot.title = element_text(hjust=0.5, face="bold")) + ggtitle("nGene distribution (Log10)") + mytheme

    # Gene counts vs. UMI counts, colored by mitochondria ratio
    p6 <- metadata %>% arrange(ratio.mt) %>%
            ggplot(aes(x=nUMI, y=nGene, color=ratio.mt)) +
            rasterize(geom_point(), dpi=300) +
            scale_colour_gradient(low = "gray90", high = "black") +
            stat_smooth(method=lm) +
            scale_x_log10(limits = c(100, 70000), n.breaks=20) +
            scale_y_log10(limits = c(100, 25000), n.breaks=20) +
            theme_classic() + mytheme  + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +theme(plot.title = element_text(hjust=0.5, face="bold"))


    # Gene counts vs. UMI counts, colored by ribosome ratio
    p7 <- metadata %>% arrange(ratio.ribo) %>%
            ggplot(aes(x=nUMI, y=nGene, color=ratio.ribo)) +
            rasterize(geom_point(), dpi=300) +
            scale_colour_gradient(low = "gray90", high = "black") +
            stat_smooth(method=lm) +
            scale_x_log10(limits = c(10, 70000), n.breaks=20) +
            scale_y_log10(limits = c(100, 25000),n.breaks=20) +
            theme_classic() + mytheme  + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +theme(plot.title = element_text(hjust=0.5, face="bold"))

    # Distribution of mitochondria ratio
    p8 <- metadata %>%
            ggplot(aes(x=ratio.mt)) + geom_density(size=1) +
            scale_x_log10(n.breaks=20) +
            theme_classic() + mytheme + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

    # Distribution of ribosome ratio
    p9 <- metadata %>%
            ggplot(aes(x=ratio.ribo)) + geom_density(size=1) +
            scale_x_log10(n.breaks=20) +
            theme_classic() + mytheme + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    # Distribution of gene per UMI (log10)
    p10 <- metadata %>%
            ggplot(aes(x=log10GenesPerUMI)) +
            geom_density(size=1) +
            theme_classic() +
            geom_vline(xintercept = 0.8, color="forestgreen", linetype="longdash")+ mytheme

    figure <- ggarrange(ggarrange(p1, p_gene_counts, p2, p3, ncol=4, labels=c("", "","",""),widths = c(1, 1, 4,4)), ggarrange(p4, p5, ncol = 2, labels = c("", ""), widths=c(2,1)), ggarrange(p6, p7, ncol=2, labels=rep("nGeneVsnUMI",2)), ggarrange(p8, p9, p10, ncol=3,labels=c("mt.ratio","ribo.ratio"), widths=c(1,1,1)),nrow = 4)
    print(figure)
    dev.off()

    plot_list=list()
    ID_set <- metadata %>% pull(sample_id) %>% unique() %>% sort()
    counter=0

    for (ID in ID_set){
            counter=counter+1
            print(ID)
            p <- ggplot(subset(metadata %>% arrange(ratio.ribo),sample_id %in% ID),aes(x=nUMI, y=ratio.mt, color=ratio.ribo)) + scale_x_log10(limits = c(100, 70000), n.breaks=20) + rasterize(geom_point(size=2), dpi=300) + theme_classic() + theme(axis.text = element_text(size = 15), axis.title = element_text(size = 0), plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + ggtitle(ID) + scale_colour_gradient(low = "gray90", high = "red") + scale_y_continuous(breaks=seq(0,max(metadata$ratio.mt,na.rm = TRUE),by=0.05))
            plot_list[[counter]]=p
    }
    pdf(paste0(OUTDIR,"/", opt$sample,"_QC_UMI_vs_mito_",min_cell,".pdf"), width=6, height=6)
    if (length(plot_list) %% 4 ==0 ){
        nrow=length(plot_list)/4
    } else if (length(plot_list) <= 4){
        nrow=1
    } else {
        nrow=length(plot_list)%/%4 + 1
    }

    if (length(plot_list) <=4 ){
        ncol=length(plot_list)
    } else {
        ncol=4
    }
    figure <- ggarrange(plotlist=plot_list, nrow = nrow, ncol = ncol)
    annotated_figure <- annotate_figure(figure,
                    bottom = text_grob("nUMI", color = "black", face = "bold", size = 15),
                    left = text_grob("ratio.mt", color = "black", face="bold",size=15, rot = 90))
    print(annotated_figure)
    dev.off()

if (min_cell > 0){
    saveRDS(s,paste0(OUTDIR,"/",opt$sample, "_BgRemoved_SeuratObject_min_cell_", min_cell, ".RDS"))
}
}

#plotQC(s,min_cell=0)
plotQC(s_filtered,opt$min.cells)
