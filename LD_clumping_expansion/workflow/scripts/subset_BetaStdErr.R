suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
option.list <- list(
  make_option("--beta", type="character"),
  make_option("--variant_list", type="character"),
  make_option("--trait", type="character"),
  make_option("--outfile", type="character")
)
opt <- parse_args(OptionParser(option_list=option.list))

var_df <- read.table(opt$variant_list, header=T, sep="\t", stringsAsFactors=F) %>% select(LocusID, variant) %>% rename(LeadVariant=LocusID, rsid=variant) %>% distinct()
beta_df <- read.table(opt$beta, header=T, sep="\t", stringsAsFactors=F) %>% filter(rsid %in% unique(c(var_df$LeadVariant, var_df$rsid))) %>% mutate(CI_lower=Beta-1.96*StdErr, CI_upper=Beta+1.96*StdErr)
beta_df_lead <- beta_df %>% filter(rsid %in% unique(c(var_df$LeadVariant))) %>% rename(LeadVariant=rsid, Lead_Beta=Beta, Lead_CI_lower=CI_lower, Lead_CI_upper=CI_upper, Lead_A1=A1, Lead_A2=A2)
beta_df_var <- beta_df %>% filter(rsid %in% unique(c(var_df$rsid))) %>% rename(Var_Beta=Beta, Var_CI_lower=CI_lower, Var_CI_upper=CI_upper, Var_A1=A1, Var_A2=A2)
var_df <- var_df %>% left_join(beta_df_lead)
var_df <- var_df %>% left_join(beta_df_var) %>% mutate(Trait=opt$trait)
write.table(var_df, gzfile(opt$outfile), sep="\t", quote=F, row.names=F)

