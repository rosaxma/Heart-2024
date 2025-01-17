suppressPackageStartupMessages(library(optparse))

option.list <- list(
    make_option("--effectiveSheet", type="character", default="", help="Output directory"), 
    make_option("--assignmentSheet", type="character", default=NA),
    make_option("--outdir", type="character"),
    make_option("--sample", type="character")
)

opt <- parse_args(OptionParser(option_list=option.list))
suppressPackageStartupMessages(library(tidyverse))

effectiveSheet <- read.table(opt$effectiveSheet, header=T, sep="\t", stringsAsFactors=F) %>% filter(EffectiveCells=="True") %>% mutate_if(is.numeric, as.character)
assignmentSheet <- read.table(opt$assignmentSheet, header=T, sep="\t", stringsAsFactors=F) %>% mutate_if(is.numeric, as.character)
if ("assignment" %in% colnames(effectiveSheet)){
    effectiveSheet <- effectiveSheet %>% left_join(assignmentSheet, by="assignment") %>% select(-assignment) %>% mutate(orig.ident=opt$sample) %>% relocate(orig.ident, .before=source)
} else{
    effectiveSheet <- effectiveSheet %>% mutate(assignment="0") %>% left_join(assignmentSheet, by="assignment") %>% select(-assignment) %>% mutate(orig.ident=opt$sample) %>% relocate(orig.ident, .before=source)
}

write.table(effectiveSheet, file.path(opt$outdir, paste0(opt$sample, "_doublet_summary_sheet_with_sample_origin.tsv")), sep="\t", quote=F, row.names=F)
