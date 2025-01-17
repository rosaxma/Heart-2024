library(optparse)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggpubr)
option.list = list(
    make_option("--pathwayCount", type="character"),
    make_option("--outdir", type="character"),
    make_option("--sample", type="character")
)
opt = parse_args(OptionParser(option_list=option.list))

pathway_count_df= opt$pathwayCount %>% strsplit(split=",") %>% unlist()
HM_K_count = as.data.frame(matrix(nrow=length(pathway_count_df),ncol=5))
colnames(HM_K_count) = c("K", "nHypergeometric", "nGSEA","nHypergeometric_unique","nGSEA_unique")
counter=0
for (df in pathway_count_df){
    counter = counter +1
    countTable <- read.table(df, header=TRUE, sep="\t",stringsAsFactors=FALSE, row.names=1)
    HM_K_count[counter, "K"] = as.numeric(unique(countTable$K))
    HM_K_count[counter, "nHypergeometric"] = as.numeric(countTable["Hypergeometric","HallmarkGeneSets"])
    HM_K_count[counter, "nGSEA"] =  as.numeric(countTable["GSEA","HallmarkGeneSets"])
    HM_K_count[counter, "nHypergeometric_unique"] =  as.numeric(countTable["Hypergeometric","Hallmark_AverageUniquePathwaysPerGEP"])
    HM_K_count[counter, "nGSEA_unique"] =  as.numeric(countTable["GSEA","Hallmark_AverageUniquePathwaysPerGEP"])
}


HM_plot_list=vector(mode="list", length=as.numeric(4))
HM_plot_list[[1]] <- ggplot(data=HM_K_count %>% arrange(as.numeric(K)), aes(x=as.numeric(K), y=nHypergeometric)) + geom_line() +
geom_point() + theme_bw() + xlab("K") + ylab("Number of unique pathways (Hypergeometric)")
HM_plot_list[[2]] <- ggplot(data=HM_K_count%>% arrange(as.numeric(K)), aes(x=as.numeric(K), y=nGSEA)) + geom_line() +
geom_point() + theme_bw() + xlab("K") + ylab("Number of unique pathways (GSEA)")
HM_plot_list[[3]] <- ggplot(data=HM_K_count%>% arrange(as.numeric(K)), aes(x=as.numeric(K), y=nHypergeometric_unique)) + geom_line() +
geom_point() + theme_bw() + xlab("K") + ylab("Number of average unique pathways per program (Hypergeometric)")
HM_plot_list[[4]] <- ggplot(data=HM_K_count%>% arrange(as.numeric(K)), aes(x=as.numeric(K), y=nGSEA_unique)) + geom_line() +
geom_point() + theme_bw() + xlab("K") + ylab("Number of average unique pathways per program(GSEA)")

GO_K_count = as.data.frame(matrix(nrow=length(pathway_count_df),ncol=5))
colnames(GO_K_count) = c("K", "nHypergeometric", "nGSEA","nHypergeometric_unique","nGSEA_unique")
counter=0
for (df in pathway_count_df){
    counter = counter +1
    countTable <- read.table(df, header=TRUE, sep="\t",stringsAsFactors=FALSE, row.names=1)
    GO_K_count[counter, "K"] = as.numeric(unique(countTable$K))
    GO_K_count[counter, "nHypergeometric"] = as.numeric(countTable["Hypergeometric","GO"])
    GO_K_count[counter, "nGSEA"] =  as.numeric(countTable["GSEA","GO"])
    GO_K_count[counter, "nHypergeometric_unique"] =  as.numeric(countTable["Hypergeometric","GO_AverageUniqueGOPerGEP"])
    GO_K_count[counter, "nGSEA_unique"] =  as.numeric(countTable["GSEA","GO_AverageUniqueGOPerGEP"])
}


GO_plot_list=vector(mode="list", length=as.numeric(4))
GO_plot_list[[1]] <- ggplot(data=GO_K_count %>% arrange(as.numeric(K)), aes(x=as.numeric(K), y=nHypergeometric,label=K)) + geom_line() + geom_point() + theme_bw() + xlab("K") + ylab(paste0("Number of unique pathways", "\n" , "(Hypergeometric)"))+theme(axis.text=element_text(size=14),axis.title=element_text(size=12,face="bold")) + geom_text(vjust = -0.5)+ scale_y_continuous(breaks = seq(800, 3000, by=200), limits=c(800,3000))

GO_plot_list[[2]] <- ggplot(data=GO_K_count %>% arrange(as.numeric(K)), aes(x=as.numeric(K), y=nGSEA,label=K)) + geom_line() +
geom_point() + theme_bw() + xlab("K") + ylab(paste0("Number of unique pathways", "\n", "(GSEA)"))+theme(axis.text=element_text(size=14),axis.title=element_text(size=12,face="bold"))+ geom_text(vjust = -0.5)+ scale_y_continuous(breaks = seq(800, 3000, by=200), limits=c(800,3000))


GO_plot_list[[3]] <- ggplot(data=GO_K_count%>% arrange(as.numeric(K)), aes(x=as.numeric(K), y=nHypergeometric_unique,label=K)) + geom_line() +geom_point() + theme_bw() + xlab("K") + ylab(paste0("Number of average unique pathways", "\n", "per program (Hypergeometric)"))+theme(axis.text=element_text(size=14),axis.title=element_text(size=12,face="bold"))+ geom_text(vjust = -0.5)+ scale_y_continuous(breaks = seq(0, 300, by=50), limits=c(0,300))


GO_plot_list[[4]] <- ggplot(data=GO_K_count%>% arrange(as.numeric(K)), aes(x=as.numeric(K), y=nGSEA_unique,label=K)) + geom_line() +
geom_point() + theme_bw() + xlab("K") +theme(axis.text=element_text(size=14),axis.title=element_text(size=2,face="bold"))+ ylab(paste0("Number of average unique pathways", "\n", "per program (GSEA)"))+ geom_text(vjust = -0.5)+ scale_y_continuous(breaks = seq(0, 300, by=50), limits=c(0,300))


pdf(file.path(opt$outdir, paste0(opt$sample, "_pathway_enrichment_hallmarkGeneSets.pdf")), width=length(pathway_count_df)*1, height=8)
p <- ggarrange(plotlist=HM_plot_list, ncol=2, nrow=2)
annotate_figure(p, top = text_grob(opt$sample, face = "bold", size = 18))
print(p)
dev.off()

pdf(file.path(opt$outdir, paste0(opt$sample, "_pathway_enrichment_GO.pdf")), width=length(pathway_count_df)*1, height=8)
p <- ggarrange(plotlist=GO_plot_list, ncol=2, nrow=2)
annotate_figure(p, top = text_grob(opt$sample, face = "bold", size = 18))
print(p)
dev.off()

#####################################################
GO_K_count_sorted=GO_K_count %>% arrange(as.numeric(K))
GO_K_count_delta <- as.data.frame(matrix(NA, nrow=length(GO_K_count_sorted$K)-1, ncol=0))
GO_K_count_delta$curr_K=GO_K_count_sorted$K[-length(GO_K_count_sorted$K)]
GO_K_count_delta$next_K=GO_K_count_sorted$K[-1]
GO_K_count_delta$delta_nHypergeometric=diff(GO_K_count_sorted$nHypergeometric)
GO_K_count_delta$delta_nGSEA=diff(GO_K_count_sorted$nGSEA)
GO_K_count_delta$delta_nHypergeometric_unique=diff(GO_K_count_sorted$nHypergeometric_unique)
GO_K_count_delta$delta_nGSEA_unique=diff(GO_K_count_sorted$nGSEA_unique)
GO_K_count_delta= GO_K_count_delta %>% mutate(delta=paste0(next_K, "-", curr_K)) %>% relocate(delta) %>% select(-c(curr_K, next_K)) %>% mutate(delta_nHypergeometric_unique=round(delta_nHypergeometric_unique, digits=2), delta_nGSEA_unique=round(delta_nGSEA_unique, digits=2))
write.table(GO_K_count_sorted, file.path(opt$outdir,paste0(opt$sample, "_pathway_enrichment_GO.tsv")), row.names=F, quote=F, sep="\t")
write.table(GO_K_count_delta, file.path(opt$outdir,paste0(opt$sample, "_pathway_enrichment_GO_delta.tsv")), row.names=F, quote=F, sep="\t")
#####################################################
GO_plot_list=vector(mode="list", length=as.numeric(4))
GO_plot_list[[1]] <- ggplot(data=GO_K_count_delta, aes(x=delta, y=delta_nHypergeometric,label=delta_nHypergeometric)) + geom_path(group=1) + geom_point() + theme_bw() + xlab("Delta") + ylab(paste0("Delta_Number of unique pathways", "\n" , "(Hypergeometric)"))+theme(axis.text=element_text(size=14),axis.title=element_text(size=12,face="bold")) + geom_text(vjust = 1, size=6)+ scale_y_continuous(breaks = seq(-150, 500, by=50), limits=c(-150,500))+ geom_hline(yintercept=0,linetype="dashed", color = "red")

GO_plot_list[[2]] <- ggplot(data=GO_K_count_delta, aes(x=delta, y=delta_nGSEA,label=delta_nGSEA)) + geom_path(group=1) +
geom_point() + theme_bw() + xlab("Delta") + ylab(paste0("Delta_Number of unique pathways", "\n", "(GSEA)"))+theme(axis.text=element_text(size=14),axis.title=element_text(size=12,face="bold"))+ geom_text(vjust = 1, size=6)+ scale_y_continuous(breaks = seq(-200, 500, by=50), limits=c(-200,450))+ geom_hline(yintercept=0,linetype="dashed", color = "red")


GO_plot_list[[3]] <- ggplot(data=GO_K_count_delta, aes(x=delta, y=delta_nHypergeometric_unique,label=delta_nHypergeometric_unique)) + geom_path(group=1) +geom_point() + theme_bw() + xlab("Delta") + ylab(paste0("Delta_Number of average unique pathways", "\n", "per program (Hypergeometric)"))+theme(axis.text=element_text(size=14),axis.title=element_text(size=12,face="bold"))+ geom_text(vjust = 1, size=6) + geom_hline(yintercept=0,linetype="dashed", color = "red")+ scale_y_continuous(breaks = seq(-45, 25, by=10), limits=c(-45,25))+ geom_hline(yintercept=0,linetype="dashed", color = "red")

GO_plot_list[[4]] <- ggplot(data=GO_K_count_delta, aes(x=delta, y=delta_nGSEA_unique,label=delta_nGSEA_unique)) + geom_path(group=1)+
geom_point() + theme_bw() + xlab("Delta") +theme(axis.text=element_text(size=14),axis.title=element_text(size=12,face="bold"))+ ylab(paste0("Delta_Number of average unique pathways", "\n", "per program (GSEA)"))+ geom_text(vjust = 1, size=6)+ geom_hline(yintercept=0,linetype="dashed", color = "red")+ scale_y_continuous(breaks = seq(-150, 500, by=50), limits=c(-150,450))+ geom_hline(yintercept=0,linetype="dashed", color = "red")+ scale_y_continuous(breaks = seq(-40, 20, by=20), limits=c(-40,20))+ geom_hline(yintercept=0,linetype="dashed", color = "red")
pdf(file.path(opt$outdir, paste0(opt$sample, "_pathway_enrichment_GO_delta.pdf")), width=length(pathway_count_df)*2, height=8)
p <- ggarrange(plotlist=GO_plot_list, ncol=2, nrow=2)
annotate_figure(p, top = text_grob(opt$sample, face = "bold", size = 18))
print(p)
dev.off()