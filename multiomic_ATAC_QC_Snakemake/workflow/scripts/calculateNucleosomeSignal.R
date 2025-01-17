suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(optparse))
option.list <- list(
    make_option("--FragmentFile", type="character"),
    make_option("--nuc_signal_threshold", type="numeric", default=2),
    make_option("--nuc_signal_output_table", type="character"),
    make_option("--nuc_signal_distribution", type="character"),
    make_option("--sample", type="character")
)
# nucleosome_signal: The nucleosome signal refers to the ratio of mono-nucleosomal to nucloesome-free fragments and can also be interpreted as a signal-to-noise ratio in each cell (more details below).
# The obtained scores in this data set range from 0 to 3. As a rule of thumb, previous analysis projects chose a threshold between 2 and 4 to label low quality cells. 
opt <- parse_args(OptionParser(option_list=option.list))
frag <- fread(opt$FragmentFile, header=F, sep="\t", stringsAsFactors=F)
frag_mod <- frag %>% mutate(length=V3-V2) %>% filter(length<294) %>% mutate(state=ifelse(length<147, "NucleosomeFree", "MonoNuclesome")) %>% group_by(V4, state) %>% summarize(n=n()) %>% ungroup() %>% as.data.frame()
frag_mod_wide <- frag_mod %>% pivot_wider(names_from = state, values_from = n, values_fill = 0) %>% as.data.frame()
frag_mod_wide <- frag_mod_wide %>% mutate(ratio=ifelse(NucleosomeFree!=0, MonoNuclesome/NucleosomeFree, NA)) %>% mutate(keep=ifelse(is.na(ratio)|ratio >= opt$nuc_signal_threshold, FALSE, TRUE))

output_df <- frag_mod_wide %>% select(V4, ratio, keep)
colnames(output_df) <- c("cbc", "nuc_signal", "keep")

p <- ggplot(output_df, aes(x=nuc_signal))+ geom_density() + geom_vline(aes(xintercept=opt$nuc_signal_threshold), color="blue", linetype="dashed", size=1) + ggtitle(opt$sample)

write.table(output_df, opt$nuc_signal_output_table, sep="\t", quote=F, row.names=F)
pdf(opt$nuc_signal_distribution)
print(p)
dev.off()