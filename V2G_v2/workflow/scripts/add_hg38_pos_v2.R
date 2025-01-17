suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
option.list <- list(
    make_option("--commonbgSNPs", type="character"),
    make_option("--finemappedSNPs", type="character"),
    make_option("--cs", type="character"),
    make_option("--StandardSNPColName", type="character"),
    make_option("--StandardPosColName", type="character"),
    make_option("--StandardChrColName", type="character"),
    make_option("--finemappedSNPColName", type="character"),
    make_option("--finemappedChrColName", type="character"), 
    make_option("--csSNPColName", type="character"),
    make_option("--csChrColName", type="character"), 
    make_option("--outputVariant", type="character"),
    make_option("--outputCS", type="character"),
    make_option("--traits", type="character")
)

opt <- parse_args(OptionParser(option_list=option.list))
##################################################################
standardSNP_df <- fread(opt$commonbgSNPs, sep="\t", header=TRUE, stringsAsFactors=FALSE) %>%  rename(position := !!sym(opt$StandardPosColName)) %>% rename(chr := !!sym(opt$StandardChrColName)) %>% rename(rsid := !!sym(opt$StandardSNPColName)) %>% select(chr, position, rsid) %>% mutate(chr=ifelse(str_detect(chr, "^chr"), chr, paste0("chr",chr)))
print(head(standardSNP_df, n=2))

finemappedSNPs_df <- read.table(opt$finemappedSNPs, sep="\t", header=TRUE, stringsAsFactors=FALSE) %>% select(-c(!!sym(opt$finemappedChrColName))) %>% select(-any_of(c("chr","start", "end", "position"))) %>% rename(rsid := !!sym(opt$finemappedSNPColName))

cs_df <- read.table(opt$cs, sep="\t", header=TRUE, stringsAsFactors=FALSE) %>% select(-c(!!sym(opt$csChrColName))) %>% select(-any_of(c("chr","start", "end", "position"))) %>% rename(LeadVariant:= !!sym(opt$csSNPColName))

finemappedSNPs <- finemappedSNPs_df %>% inner_join(standardSNP_df, by=c("rsid")) %>% mutate(start=position-1) %>% rename(end=position) %>% relocate(chr) %>% relocate(start, .after=chr) %>% relocate(end, .after=start)

cs_df <- cs_df %>% inner_join(standardSNP_df, by=c("LeadVariant"="rsid")) %>% mutate(start=position-1) %>% rename(end=position) %>% relocate(chr, .after=CredibleSet) %>% relocate(start, .after=chr) %>% relocate(end, .after=start)
##################################################################
traits <- read.table(opt$traits, header=T, sep="\t", stringsAsFactors=F) %>% filter(!REMOVE) %>% pull(Trait)
write.table(finemappedSNPs %>% filter(Trait %in% traits), opt$outputVariant,sep="\t", quote=F, row.names=FALSE)
write.table(cs_df %>% filter(Trait %in% traits), opt$outputCS,sep="\t", quote=F, row.names=FALSE)
