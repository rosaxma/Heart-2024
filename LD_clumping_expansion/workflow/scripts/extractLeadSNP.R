library(tidyverse)
library(optparse)
option.list <- list(
    make_option("--LD_clump", type="character"),
    make_option("--outdir", type="character"),
    make_option("--scratchdir", type="character")
)
opt <- parse_args(OptionParser(option_list=option.list))
ld_table <- read.table(opt$LD_clump, header=T) %>% select(CHR, SNP, BP, P)
write.table(ld_table, file.path(opt$outdir, "LeadSNPs.tsv"), quote=F, sep="\t", row.names=F)
for (i in seq(1, 22)) {
	chr_table <- ld_table %>% filter(CHR==i) %>% select(SNP)
	write.table(chr_table, file.path(opt$scratchdir, paste0("chr",i, ".LeadSNPs.tsv")), col.names=F, row.names=F, quote=F)
}

