suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
option.list <- list(
  make_option("--odds_ratio", type="character"),
  make_option("--variant_list", type="character"),
  make_option("--trait", type="character"),
  make_option("--outfile", type="character")
)
opt <- parse_args(OptionParser(option_list=option.list))

var_df <- read.table(opt$variant_list, header=T, sep="\t", stringsAsFactors=F) %>% select(LocusID, variant) %>% rename(LeadVariant=LocusID, rsid=variant) %>% distinct()
odds_ratio_df <- read.table(opt$odds_ratio, header=T, sep="\t", stringsAsFactors=F) %>% filter(rsid %in% unique(c(var_df$LeadVariant, var_df$rsid)))
odds_ratio_df_lead <- odds_ratio_df %>% filter(rsid %in% unique(c(var_df$LeadVariant))) %>% rename(LeadVariant=rsid, Lead_OR=OR, Lead_CI_lower=CI_lower, Lead_CI_upper=CI_upper, Lead_A1=A1, Lead_A2=A2)
print(head(odds_ratio_df_lead))
odds_ratio_df_var <- odds_ratio_df %>% filter(rsid %in% unique(c(var_df$rsid))) %>% rename(Var_OR=OR, Var_CI_lower=CI_lower, Var_CI_upper=CI_upper, Var_A1=A1, Var_A2=A2)
print(head(odds_ratio_df_var))
print(head(var_df))
var_df <- var_df %>% left_join(odds_ratio_df_lead)
print(head(var_df))
var_df <- var_df %>% left_join(odds_ratio_df_var) %>% mutate(Trait=opt$trait)
write.table(var_df, gzfile(opt$outfile), sep="\t", quote=F, row.names=F)

