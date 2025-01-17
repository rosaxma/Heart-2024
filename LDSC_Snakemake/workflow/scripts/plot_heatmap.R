suppressPackageStartupMessages(library(optparse))
option.list <- list(
    make_option("--outdir", type="character"), 
    make_option("--mtx", type="character"),
    make_option("--annotation", type="character"),
    make_option("--lead_variant_count", type="character", default=NA),
    make_option("--infoSheet", type="character", default=NA),
    make_option("--p_value", type="character", default=NA)
)

opt <- parse_args(OptionParser(option_list=option.list))

suppressPackageStartupMessages(library(utils))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Matrix))
library(ComplexHeatmap)
library(circlize)
set.seed(42)
mtx_list=unlist(strsplit(opt$mtx, ","))
annotation <- read.table(opt$annotation, header=T, sep="\t", stringsAsFactors=F) %>% select(celltype, plot_label) %>% distinct()
infoSheet <- read.table(opt$infoSheet, header=T, sep="\t", stringsAsFactors=F) %>% select(Trait, NiceName, DataGroup, Property, EnhancerRelevance, Source)
mat_fdr <- read.table(opt$p_value, header=T, sep="\t", stringsAsFactors=F) %>% column_to_rownames(var="Category") %>% t()
print(head(mat_fdr))
for (mtx in mtx_list){
        all_df <- read.table(mtx, header=T, sep="\t", stringsAsFactors=F) 
        all_df <- all_df %>% column_to_rownames(var="Category")
        trait_annot <- data.frame(Trait=colnames(all_df)) %>% left_join(infoSheet)
        colnames(all_df)=trait_annot$NiceName
        mat <- all_df %>% drop_na() %>% as.matrix()
        mat <- mat[!rowSums(is.na(mat)),]
        mat <- t(mat)
        print(length(colnames(mat)))
        Category=colnames(mat)
        celltypes=data.frame(Category) %>% left_join(annotation, by=c("Category"="celltype")) %>% rename(coarse_cluster=plot_label) %>% select(coarse_cluster)

        filename=basename(mtx)
        filename=gsub("combined_LDSC_", "", filename)
        filename=gsub(".mtx", "", filename)
        print("check1")
        ha = HeatmapAnnotation(coarse_cluster=celltypes$coarse_cluster, col=list(coarse_cluster=c("ACM"="skyblue2", "CoreConductionCells"="orchid1", "Endocardial"="green4", "Endothelial"="#FF7F00", "Epicardial"="gold1", "FB"="yellow4", "LymphoidCells"="#FB9A99", "MyeloidCells"="steelblue4", "NC"="deeppink1", "MuralCells"="darkturquoise","SchwannCells"="blue1", "SympatheticNeuron"="dodgerblue2", "TzConductionCells"="mediumslateblue", "VCM"="khaki2")))
        print("check2")
        print(head(trait_annot %>% select(DataGroup, Property, EnhancerRelevance,Source)))
        row_ha=rowAnnotation(df=trait_annot %>% select(DataGroup, Property, EnhancerRelevance, Source), col=list(DataGroup=c("Aorta"="#FF7F00","AortaPulmonaryArtery"="orchid1","AorticValve"="green4","CHD"="gold1","DescendingAorta"="yellow4","HeartRhythm"="#FB9A99","LeftHeart"="steelblue4","Other"="darkturquoise","Other"="blue1","PulmonaryArtery"="dodgerblue2","RightHeart"="mediumslateblue","RightLeftHeart"="palegreen2","Valve"="deeppink1","Vasculature"="#E31A1C"), Property=c("Disease"="skyblue2","Structure_Function"="green4"), EnhancerRelevance=c("Adult"="#FF7F00", "Fetal"="darkturquoise", "Fetal_Adult"="mediumslateblue"),Property=c("Disease"="skyblue2","Structure_Function"="green4"), Source=c("Cordell"="#FF7F00","JPP_v1"="orchid1","JPP_v2"="green4","Other"="gold1","Priest"="yellow4","UKB"="#FB9A99", "FINNGEN"="dodgerblue2")), annotation_width = 0.2)

        quantile_table=quantile(mat, c(0.1, 0.95), na.rm=TRUE)
        quantile_low=as.numeric(quantile_table["10%"]) 
        quantile_high=as.numeric(quantile_table["95%"]) 
        if (filename=="zscore"){
               col_fun=colorRamp2(c(quantile_low,0, quantile_high), c("blue", "white", "#A50026")) 
        } else {
                col_fun=colorRamp2(c(0, quantile_high), c("white", "#A50026"))  
        }
        print(head(mat))
        for (i in (rownames(mat_fdr))) {
              for (j in (colnames(mat_fdr))) {
                                print(i)
                                print(j)
                                 if(mat_fdr[i, j] < 0.001 ) {
                                        print("*")
                                 }
              }
        }
        p <- Heatmap(mat, name=filename,
                cluster_columns = TRUE,
                show_column_dend = TRUE,
                show_row_dend = TRUE,
                cluster_rows=TRUE,
                row_names_gp = gpar(fontsize = 10),
                column_title_gp = gpar(fontsize = 2),
                show_column_names = TRUE,
                use_raster = TRUE,
                top_annotation=ha,
                right_annotation = row_ha,
                col=col_fun, 
                column_title_side = "top",
                raster_quality = 4, 
                row_dend_width=unit(30, "mm"),
                column_dend_height=unit(30, "mm"),
                cell_fun = function(j, i, x, y, w, h, fill) {
                        if(mat_fdr[i, j] < 0.001 ) {
                                grid.text("***", x, y)
                        } else if(mat_fdr[i, j] < 0.01) {
                grid.text("**", x, y)
                } else if (mat_fdr[i, j] < 0.05) {
                grid.text("*", x, y)
                }
        })
        print(dim(mat))
        pdf(file.path(opt$outdir,paste0(filename, ".heatmap.pdf")), height=nrow(mat)*0.2+4, width=ncol(mat)*0.2+4)
        print(p)      
        dev.off()
}
