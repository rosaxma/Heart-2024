suppressPackageStartupMessages(library(optparse))
option.list <- list(
    make_option("--betaScoreTable", type="character"),
    make_option("--knees", type="character"), 
    make_option("--numK", type="character"),
    make_option("--outdir", type="character")
)
opt <- parse_args(OptionParser(option_list=option.list))

library(Matrix)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(ggpubr)
library(ggrepel)
library(stringr)
suppressPackageStartupMessages(library(ggrastr))
set.seed(42)

knee_table <- read.table(opt$knees,header=TRUE, sep="\t",stringsAsFactors=FALSE)
beta_table <- read.table(opt$betaScoreTable,header=TRUE, sep="\t",stringsAsFactors=FALSE, row.names=1) %>% rownames_to_column("Gene")
colnames(knee_table) <- c("GEP", "Knee")
k=opt$numK
ncol=6
if ((as.numeric(k) %% ncol) != 0){
        nrow=as.numeric(k)/ncol + 1
} else {
        nrow=as.numeric(k)/ncol
}
plot_list=vector(mode="list", length=as.numeric(k))
beta_colnames=grep("GEP", colnames(beta_table),value=TRUE)

max_length=max(knee_table$Knee)
top_beta_score <- data.frame(matrix(ncol = length(beta_colnames)*2, nrow = max_length))
colnames(top_beta_score) <- sort(c(paste0(beta_colnames, "_genes"), beta_colnames))

top_100_beta_score <- data.frame(matrix(ncol = length(beta_colnames)*2, nrow = 100))
top_200_beta_score <- data.frame(matrix(ncol = length(beta_colnames)*2, nrow = 200))
top_300_beta_score <- data.frame(matrix(ncol = length(beta_colnames)*2, nrow = 300))
colnames(top_100_beta_score) <- sort(c(paste0(beta_colnames, "_genes"), beta_colnames))
colnames(top_200_beta_score) <- sort(c(paste0(beta_colnames, "_genes"), beta_colnames))
colnames(top_300_beta_score) <- sort(c(paste0(beta_colnames, "_genes"), beta_colnames))

for (n in 1:length(beta_colnames)){
        topic <- beta_colnames[[n]]
        topic_df <- beta_table %>% select(!!sym(topic), Gene) 
        colnames(topic_df) <- c("Beta", "Gene")
        topic_df <- topic_df %>% mutate(rank=rank(-Beta, na.last=TRUE,ties.method="first")) %>% mutate(rank_rev=rank(Beta, na.last=TRUE, ties.method="first")) %>% filter(!is.na(Gene))
	knee=knee_table[n, "Knee"]
        p <- ggplot(topic_df, aes(x=rank, y=Beta, label=Gene)) + rasterize(geom_point(), dpi=300) + ggtitle(topic) + geom_text_repel(data=subset(topic_df, rank<20 | rank_rev <20),max.overlaps=15,nudge_x=0.0003, nudge_y=0.0001 ) + geom_vline(xintercept = knee, color="blue", size=1)
        plot_list[[n]] <- p
        ####################################
        top_beta_df <- topic_df %>% filter(rank < knee) %>% arrange(desc(Beta)) 
        top_beta_score[[paste0(topic, "_genes")]] <- c(top_beta_df$Gene, rep(NA, max_length-length(top_beta_df$Gene)))
        top_beta_score[[topic]] <- c(top_beta_df$Beta, rep(NA, max_length-length(top_beta_df$Beta)))
        ####################################
        top_beta_df_100 <- topic_df %>% filter(rank <= 100) %>% arrange(desc(Beta))
        top_100_beta_score[[paste0(topic, "_genes")]] <- top_beta_df_100$Gene
        top_100_beta_score[[topic]] <- top_beta_df_100$Beta

        top_beta_df_200 <- topic_df %>% filter(rank <= 200) %>% arrange(desc(Beta))
        top_200_beta_score[[paste0(topic, "_genes")]] <- top_beta_df_200$Gene
        top_200_beta_score[[topic]] <- top_beta_df_200$beta

        top_beta_df_300 <- topic_df %>% filter(rank <= 300) %>% arrange(desc(Beta))
        top_300_beta_score[[paste0(topic, "_genes")]] <- top_beta_df_300$Gene
        top_300_beta_score[[topic]] <- top_beta_df_300$Beta
    }

pdf(file.path(opt$outdir, "figures", paste0("K_", k, "_Beta_distribution.pdf")), height=nrow*6, width=ncol*6)
p <- ggarrange(plotlist=plot_list, ncol=ncol, nrow=nrow)
print(p)
dev.off()

write.table(top_beta_score, file.path(opt$outdir,"tables", paste0("K_", k, "_topGenes_in_topics.tsv")), quote=F, sep="\t", row.names=F)
write.table(top_100_beta_score, file.path(opt$outdir,"tables", paste0("K_", k, "_top100Genes_in_topics.tsv")), quote=F, sep="\t", row.names=F)
write.table(top_200_beta_score, file.path(opt$outdir,"tables", paste0("K_", k, "_top200Genes_in_topics.tsv")), quote=F, sep="\t", row.names=F)
write.table(top_300_beta_score, file.path(opt$outdir,"tables", paste0("K_", k, "_top300Genes_in_topics.tsv")), quote=F, sep="\t", row.names=F)

