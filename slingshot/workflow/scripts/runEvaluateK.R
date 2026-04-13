suppressPackageStartupMessages(library(optparse))
library(tidyverse)
library(slingshot)
library(Seurat)
library(tradeSeq)
library(cowplot)
library(BiocParallel)
set.seed(42)
option.list <- list(
    make_option("--sce_object", type="character", default=""),
    make_option("--outdir", type="character", default=""),
    make_option("--celltype", type="character", default=NA)
)
opt <- parse_args(OptionParser(option_list=option.list))
sce <- readRDS(opt$sce_object)
rv <- rowVars(assay(sce, "counts"))
top_indices <- order(rv, decreasing=TRUE, na.last=TRUE)[1:2000]
ic_genes <- rownames(sce)[top_indices]
pst <- slingPseudotime(sce)
weights <- slingCurveWeights(sce) 
for (lin in colnames(pst)){ 
	lineage_cells <- weights[, lin] > 0
	if(sum(lineage_cells) > 100) {
		lineage_counts <- assays(sce)$counts[ic_genes, lineage_cells]
		prop_expressed <- rowMeans(lineage_counts > 0)
		active_genes <- names(prop_expressed[prop_expressed > 0.05])
		if(length(active_genes) > 500) {
            		lineage_vars <- rowVars(as.matrix(lineage_counts[active_genes, ]))
            		names(lineage_vars) <- active_genes
            		active_genes <- names(sort(lineage_vars, decreasing = TRUE))[1:500]
        		}	
		if(length(active_genes) < 100) {
            		message("Too few active genes for stable evaluation. Skipping.")
            		next
        		}	
			
		sub_counts <- as.matrix(lineage_counts[active_genes, ])
        	sub_pst <- pst[lineage_cells, lin]
        	sub_weights <- weights[lineage_cells, lin]
		message(paste("Evaluating Lineage:", lin, "with", ncol(sub_counts), "cells"))			
		eval_results <- evaluateK(counts = sub_counts,
                          pseudotime = sub_pst, 
                          cellWeights = sub_weights,  
                          k = 3:10, 
                          nGenes = length(active_genes)) 
	saveRDS(eval_results, file = file.path(opt$outdir, paste0("tradeSeq_evaluateK_results_", opt$celltype, "_", lin, ".rds")))
	eval_table <- as.data.frame(eval_results)
	eval_table$gene <- rownames(eval_table)
	write.table(eval_table, file.path(opt$outdir, paste0("tradeSeq_knot_metrics_", opt$celltype, "_", lin, ".tsv")), row.names = FALSE, sep="\t", quote=F)
	}
}
