suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))

option.list <- list(
  make_option("--combined_E2G_info", type="character"),
  make_option("--combined_enh", type="character"),
  make_option("--outfile", type="character")
)
opt <- parse_args(OptionParser(option_list=option.list))

df_enh <- read.table(opt$combined_enh, header=F, sep="\t", stringsAsFactors=F)
colnames(df_enh) <- c("merged_chr", "merged_enh_start", "merged_enh_end", "enhancers")

df_enh <- df_enh %>% separate_longer_delim(enhancers, delim=";") %>% separate_wider_delim(enhancers, delim="|", names=c("enh_chr", "enh_start", "enh_end")) %>% mutate(enh_start=as.numeric(enh_start)) %>% mutate(enh_end=as.numeric(enh_end))

table <- read.table(opt$combined_E2G_info, header=T, sep="\t", stringsAsFactors=F)

table <- table %>% left_join(df_enh) %>% relocate(merged_enh_end) %>% relocate(merged_enh_start) %>% relocate(merged_chr)

write.table(table, opt$outfile, row.names=F, sep="\t", quote=F)
