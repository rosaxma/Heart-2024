suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))

option.list <- list(
    make_option("--sumstats", type="character"),
    make_option("--LD_expand", type="character"),
    make_option("--outdir", type="character"),
    make_option("--trait", type="character")
)
opt <- parse_args(OptionParser(option_list=option.list))
##############################################################################
project=opt$trait
##############################################################################
## Load GWAS data and credible sets
LD_expand <- read.table(opt$LD_expand, header=T, stringsAsFactors=F)
LD_expand  <- LD_expand %>%
  mutate(
    LocusID=LEADSNP,
    Disease=project,
    Variant=SNP_B,
    Chr=paste0("chr",CHR),
    Position=as.numeric(BP_B),
    ) %>% distinct() %>% select(LocusID, Disease, Variant, Chr, Position)
##############################################################################
sumstats_table <- read.table(opt$sumstats, header=T, sep="\t", stringsAsFactors=F) %>% select(SNP, P)
if (nrow(LD_expand) != 0){
  LD_expand <- LD_expand %>% left_join(sumstats_table, by=c("Variant"="SNP")) %>% rename(variant_P_value=P) %>% left_join(sumstats_table, by=c("LocusID"="SNP")) %>% rename(LeadSNP_P_value=P)
}
##############################################################################
write.table(LD_expand, file=file.path(opt$outdir, "expanded.variant.txt"),sep="\t", quote=F, row.names=F)
