suppressPackageStartupMessages(library(optparse))
library(SeuratObject)
library(Seurat)
library(Matrix)
library(optparse)
library(dplyr)
library(tidyr)
option.list <- list(
    make_option("--SeuratObject", type="character"),
    make_option("--usageScoreTable", type="character"),
    make_option("--features", type="character"),
    make_option("--numK", type="character"),
    make_option("--outdir", type="character"),
    make_option("--additionalCellInfo", type="character", default=NA),
    make_option("--cbc_colname", type="character", default=NA),
    make_option("--additional_metrics", type="character", default=NA)
)
opt <- parse_args(OptionParser(option_list=option.list))
s <- readRDS(opt$SeuratObject)

feature_list <- opt$features%>% strsplit(split=",") %>% unlist()
print(s@meta.data%>%head(n=2))
og_metadata <- s@meta.data %>% select(-cbc) 
usage_table <- read.table(opt$usageScoreTable,header=TRUE, sep="\t",stringsAsFactors=FALSE, row.names=1) %>% tibble :: rownames_to_column("cbc")
print(head(usage_table))
print(feature_list)
for (feature in feature_list){
    cluster_assignment=feature
    if (cluster_assignment %in% colnames(og_metadata)){
        metadata <- og_metadata %>% tibble :: rownames_to_column("cbc") %>% select(cbc, orig.ident, !!sym(cluster_assignment)) %>% mutate(!!sym(cluster_assignment) := paste0("cluster_", !!sym(cluster_assignment)))
        table <- metadata %>% left_join(usage_table, by="cbc")
        print(head(table))
        if (!is.na(opt$additionalCellInfo)){
            metrics <- opt$additional_metrics %>% strsplit(split=",") %>% unlist()
            #at least contain "cbc","orig.ident"
            cell_info <- read.table(opt$additionalCellInfo, header=T, sep="\t", stringsAsFactors=F) %>% mutate(cbc=paste0(orig.ident, "_", !!sym(opt$cbc_colname)))
            print(head(cell_info))
            table <- table %>% left_join(cell_info,by=c("cbc","orig.ident"))
        }

print(dim(table))
print(head(table))
##########################################################
        nk=as.numeric(opt$numK)
        usage_columns=grep("GEP_", colnames(table),value=TRUE)
        if (!is.na(opt$additionalCellInfo)){
            nrow=length(unique(metadata$orig.ident))+length(unique(metadata[[cluster_assignment]])) 
            rownames=c(unique(metadata$orig.ident), unique(metadata[[cluster_assignment]]))

            nrow_2=length(unique(metadata$orig.ident))
            rownames_2=unique(metadata$orig.ident)
            for (metric in metrics){
		        print(metric)
                nrow=nrow+length(unique(table %>% pull(!!sym(metric))))
                rownames=c(rownames, paste0(metric, "_", unique(table %>% pull(!!sym(metric)))))

                nrow_2=nrow_2+length(unique(table %>% pull(!!sym(metric))))
                rownames_2=c(rownames_2, paste0(metric, "_", unique(table %>% pull(!!sym(metric)))))
            }

        } else {
            nrow=length(unique(metadata$orig.ident))+length(unique(metadata[[cluster_assignment]]))
            rownames=c(unique(metadata$orig.ident), unique(metadata[[cluster_assignment]]))

            nrow_2=length(unique(metadata$orig.ident))
            rownames_2=unique(metadata$orig.ident)
        }
        print(rownames)
        topic_assignment_df=setNames(data.frame(matrix(nrow = nrow, ncol= nk)), usage_columns)
        rownames(topic_assignment_df) <- rownames
        topic_expression_df=setNames(data.frame(matrix(nrow = nrow, ncol= nk)), usage_columns)
        rownames(topic_expression_df) <- rownames
        topic_content_df=setNames(data.frame(matrix(nrow = nrow_2, ncol= nk)), usage_columns)
        rownames(topic_content_df) <- rownames_2


        
        for (k in usage_columns){
            # filter for usage larger than 10%
            gep_cells = table %>% filter(as.numeric(!!sym(k)) > 10) %>% pull(cbc)
            for (celltype in unique(table %>% pull(!!sym(cluster_assignment)))) {
                celltype=as.character(celltype)
                print(table %>% filter(!!sym(cluster_assignment)==celltype))
                cells = table %>% filter(!!sym(cluster_assignment)==celltype) %>% pull(cbc)
                shared_cells = intersect(cells, gep_cells)
                percentage = length(shared_cells)/length(cells)*100
                avg_expression = sum(table %>% filter(cbc %in% shared_cells) %>% pull(!!sym(k)))/length(shared_cells)
                topic_assignment_df[celltype, k] = percentage
                topic_expression_df[celltype, k] = avg_expression
                gep_content=length(shared_cells)/length(gep_cells)*100
                topic_content_df[celltype, k] =gep_content 
            }
            for (ident in unique(table$orig.ident)){
                cells = table %>% filter(orig.ident==ident) %>% pull(cbc)
                shared_cells= intersect(cells, gep_cells)
                percentage = length(shared_cells)/length(cells)*100
                avg_expression = sum(table %>% filter(cbc %in% shared_cells) %>% pull(!!sym(k)))/length(shared_cells)
                topic_assignment_df[ident, k] = percentage
                topic_expression_df[ident, k] = avg_expression

                gep_content=length(shared_cells)/length(gep_cells)*100
                topic_content_df[ident, k] =gep_content
            }
            if (!is.na(opt$additionalCellInfo)) {
            	print("additional")
            	metrics <- opt$additional_metrics %>% strsplit(split=",") %>% unlist()
            	for (metric in metrics){
                	print(metric)
                	print(unique(table[[metric]]))
                	for (item in unique(table[[metric]])){
                	print(item)
                	cells = table %>% filter(!!sym(metric)==item) %>% pull(cbc)
                	shared_cells = intersect(cells, gep_cells)
                	percentage = length(shared_cells)/length(cells)*100
                	avg_expression = sum(table %>% filter(cbc %in% shared_cells) %>% pull(!!sym(k)))/length(shared_cells)
                	gep_content=length(shared_cells)/length(gep_cells)*100

                	topic_assignment_df[paste0(metric, "_", item), k] = percentage
                	topic_expression_df[paste0(metric, "_", item), k] = avg_expression
                	topic_content_df[paste0(metric, "_", item), k] =gep_content
                	}
            	}
            	}
        }
        write.table(topic_assignment_df, file.path(opt$outdir, paste0("k_", as.character(opt$numK), "_topic_in_cell_feature_", as.character(feature),".tsv")), row.names=TRUE, quote=F, sep="\t", col.names=NA)
        print("save")
        write.table(topic_expression_df, file.path(opt$outdir, paste0("k_", as.character(opt$numK), "_topic_exp_in_cell_feature_", as.character(feature), ".tsv")), row.names=TRUE, quote=F, sep="\t", col.names=NA)
        write.table(topic_content_df, file.path(opt$outdir, paste0("k_", as.character(opt$numK), "_topic_content_by_feature_", as.character(feature), ".tsv")), row.names=TRUE, quote=F, sep="\t", col.names=NA)
        write.table(table, file.path(opt$outdir, paste0("k_", as.character(opt$numK), "_metatable.tsv")), row.names=TRUE, quote=F, sep="\t", col.names=NA)
    }
}
write.table("job_complete", file.path(opt$outdir,"job_complete.tsv"))
