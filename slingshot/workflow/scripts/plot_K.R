library(ggplot2)
library(optparse)
library(tools)

option.list <- list(
    make_option("--outdir", type="character", default=""),
    make_option("--celltype", type="character", default=NA)
)
opt <- parse_args(OptionParser(option_list=option.list))
celltype_var <- opt$celltype
pattern_str <- paste0("^tradeSeq_knot_metrics_", celltype_var, "_Lineage\\d+\\.tsv$")

files <- list.files(path = opt$outdir, 
                    pattern = pattern_str, 
                    full.names = TRUE)

for (f in files) {
    current_data <- read.table(f, header = TRUE, sep = "\t")
    name_label <- file_path_sans_ext(basename(f))
    eval_table <- read.table(f, header = TRUE, sep = "\t", check.names = FALSE)
    k_cols <- setdiff(colnames(eval_table), "gene")
    mean_AIC <- colMeans(eval_table[, k_cols], na.rm = TRUE)
    k_values <- as.numeric(gsub("k: ", "", names(mean_AIC)))
    plot_data <- data.frame(k = k_values, AIC = mean_AIC)
    p <- ggplot(plot_data, aes(x = factor(k), y = AIC, group = 1)) +
         geom_line(color = "steelblue", size = 1) +
         geom_point(size = 3) +
         theme_minimal() +
         labs(title = paste(name_label),
         x = "Number of Knots (K)",
         y = "Mean AIC") 
    pdf(file.path(opt$outdir, paste0(name_label, ".pdf")))
    print(p)
    dev.off()
}
