suppressPackageStartupMessages(library(optparse))
library(tidyverse)
library(slingshot)
library(Seurat)
library(tradeSeq)
library(cowplot)
library(BiocParallel)
set.seed(42)
option.list <- list(
    make_option("--inputSeuratObject", type="character", default=""),
    make_option("--sce_object", type="character", default=""),
    make_option("--outdir", type="character", default=""),
    make_option("--sample", type="character", default=NA),
	make_option("--geneList", type="character", default=NA),
	make_option("--knots_table", type="character", default=NA)
)
opt <- parse_args(OptionParser(option_list=option.list))
gene_list <- read.table(opt$geneList, header=T, sep="\t", stringsAsFactors=F) %>% pull(Gene) %>% unique()
sce <- readRDS(opt$sce_object)
n_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1))
BPPARAM <- MulticoreParam(workers = n_cores)
pst <- slingPseudotime(sce)
weights <- slingCurveWeights(sce)
knots_table <- read.table(opt$knots_table, header=T, sep="\t", stringsAsFactors=F)
lin_list=knots_table %>% filter(celltype==opt$sample) %>% pull(lineage) %>% unique()
object <- readRDS(opt$inputSeuratObject)
counts <- as.matrix(object@assays$RNA@counts)
genes_in_data <- intersect(gene_list, rownames(sce))
full_counts <- counts(sce)

plot_differential_expression <- function(feature_id, sce_obj, original_sce, lin_idx, cell_idx) {
    pdf(file.path(opt$outdir, paste0(opt$sample, "_", feature_id, "_lin", lin_idx, "_association.pdf")), width = 12, height = 5)
    sce_sub <- original_sce[, cell_idx]
    sds_sub <- SlingshotDataSet(sce_sub)
    p1 <- plotGeneCount(sds_sub, 
                        counts = assays(sce_sub)$counts, 
                        gene = feature_id) + 
          ggplot2::ggtitle(paste("Lineage", lin_idx, ":", feature_id)) +
          ggplot2::theme_minimal() +
          ggplot2::theme(legend.position = "right")
    p2 <- plotSmoothers(sce_obj, 
                        counts = assays(sce_sub)$counts, 
                        gene = feature_id) +
          ggplot2::ggtitle("Fitted GAM Progression") +
          ggplot2::theme_minimal()
    print(cowplot::plot_grid(p1, p2, ncol = 2, rel_widths = c(1.2, 1)))
    dev.off()
}

library_sizes <- colSums(assays(sce)$counts)
out_list <- list()
for (lin in lin_list){
	cell_idx <- which(weights[, lin] > 0)
	sub_counts <- as.matrix(full_counts[genes_in_data, cell_idx])
	sub_pst <- pst[cell_idx, lin]
    	sub_weights <- weights[cell_idx, lin]

	prop_expressed <- rowMeans(sub_counts > 0)
	genes_to_test <- names(prop_expressed[prop_expressed > 0.05])
	sub_offset  <- log(library_sizes[cell_idx])

	nknots_val <- knots_table %>% 
        filter(celltype == opt$sample, lineage == lin) %>% 
        pull(knot) %>% as.numeric()

	# returns a list of gamObject for each gene for plotting the heatmap
	sce_list <- fitGAM(counts = sub_counts[genes_to_test, , drop=FALSE],
				pseudotime = sub_pst,
				cellWeights = sub_weights,
				offset = sub_offset,
				nknots = nknots_val,
				parallel = FALSE, sce=FALSE)
	saveRDS(sce_list, file.path(opt$outdir, paste0(opt$sample, "_fit_GAM_", lin, ".RDS")))

	# returns a singleCellExperiment object 
	sce_lin <- fitGAM(counts = sub_counts[genes_to_test, , drop=FALSE],
                      pseudotime = sub_pst, 
                      cellWeights = sub_weights,
					  offset = sub_offset,
                      nknots = nknots_val,
                      parallel = TRUE, BPPARAM = BPPARAM)

	res <- associationTest(sce_lin, lineages=FALSE) %>% rownames_to_column(var="gene")
	res$fdr <- p.adjust(res$pvalue, method = "fdr")
	sig_genes <- res %>% filter(fdr < 0.05) %>% pull(gene)

	sce_subset <- sce[, cell_idx]
	if(!all(colnames(sce_lin) == colnames(sce_subset))) {
    		stop("Cell order mismatch detected between fitGAM output and original subset!")
	} else {
		message("Cell orders aligned")
	}
	metadata(sce_lin)$slingshot <- metadata(sce_subset)$slingshot
	lin_idx <- as.numeric(gsub("\\D", "", lin))
	out_df <- res %>% mutate(lineage=lin)
	out_list[[lin]] = out_df
	for (g in sig_genes) {
			plot_differential_expression(g, sce_lin, sce, lin_idx, cell_idx)

	}
	}
	combined_df <- bind_rows(out_list) 
	final_output_path <- file.path(opt$outdir, paste0(opt$sample, "_tradeSeq_results.tsv"))
	if (nrow(combined_df) > 0) {
		dev_values <- sapply(sce_list, function(m) {
			return(summary(m)$dev.expl)
        })
        dev_df <- data.frame(
                    gene = names(sce_list),
                    deviance_explained = dev_values, stringsAsFactors=FALSE
                    )
		combined_df <- combined_df %>% left_join(dev_df, by="gene") %>% filter(deviance_explained > 0.05) %>% arrange(fdr)
		write.table(combined_df, final_output_path, sep="\t", quote=F, row.names=F)
	} else {
		message("No significant genes found across any lineages.")
	}
