library(optparse)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggpubr)
option.list = list(
    make_option("--varianceExplained", type="character"),
    make_option("--outputFigure", type="character"),
    make_option("--outputFigure_delta", type="character"),
    make_option("--outputTable", type="character"),
    make_option("--outputTable_delta", type="character"),
    make_option("--sample", type="character")
)
opt = parse_args(OptionParser(option_list=option.list))

varianceExplained_df = opt$varianceExplained %>% strsplit(split=",") %>% unlist()

percentage_df = as.data.frame(matrix(nrow=length(varianceExplained_df),ncol=2))
colnames(percentage_df) = c("VarianceExplained", "K")
counter=0
for (df in varianceExplained_df){
    counter = counter +1
    varianceTable <- read.table(df, header=TRUE, sep="\t",stringsAsFactors=FALSE)
    percentage_df[counter, "VarianceExplained"] = round(as.numeric(varianceTable[1,"VarianceExplained"]), digits=4)
    percentage_df[counter, "K"] = as.numeric(varianceTable[1,"K"])
}

pdf(opt$outputFigure, width=5, height=3)
p <- ggplot(data=percentage_df %>% arrange(as.numeric(K)), aes(x=as.numeric(K), y=VarianceExplained, label=K)) + geom_line() +
geom_point() + theme_bw() + xlab("K") + ylab("Variance Explained")+ geom_text(vjust = -0.5, size=6) + ggtitle(opt$sample)+ scale_y_continuous(breaks = seq(0, 0.5, by=0.1), limits=c(0.1,0.4))
print(p)
dev.off()

write.table(percentage_df, opt$outputTable, sep="\t", quote=FALSE, row.names=FALSE)

percentage_df_sorted=percentage_df %>% arrange(as.numeric(K))
delta_percentage_df <- as.data.frame(matrix(NA, nrow=length(percentage_df$K)-1, ncol=0))
delta_percentage_df$curr_K=percentage_df_sorted$K[-length(percentage_df_sorted$K)]
delta_percentage_df$next_K=percentage_df_sorted$K[-1]
delta_percentage_df$delta_variance=diff(percentage_df_sorted$VarianceExplained)
print(delta_percentage_df)
delta_percentage_df <- delta_percentage_df %>% mutate(delta=paste0(next_K, "-", curr_K))%>% relocate(delta) %>% select(-c(curr_K, next_K)) %>% mutate(delta_variance=round(delta_variance, digits=5))

pdf(opt$outputFigure_delta, width=8, height=5)
p <- ggplot(data=delta_percentage_df, aes(x=delta, y=delta_variance, label=delta_variance)) + geom_path(group=1) +
geom_point() + theme_bw() + xlab("Delta") + ylab("Variance Explained Delta")+ geom_text(vjust = -0.5, size=3) + ggtitle(opt$sample)+ geom_hline(yintercept=0,linetype="dashed", color = "red")+ scale_y_continuous(breaks = seq(-0.05, 0.05, by=0.05), limits=c(-0.05,0.05))
print(p)
dev.off()

write.table(percentage_df, opt$outputTable_delta, sep="\t", quote=FALSE, row.names=FALSE)
