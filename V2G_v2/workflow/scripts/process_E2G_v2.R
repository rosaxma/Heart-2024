suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))
option.list <- list(
	make_option("--E2GOverlap", type="character"),
    make_option("--outdir", type="character"),
    make_option("--removeNonCoding", type="logical", default=TRUE),
    make_option("--UBQ_genes", type="character"),
    make_option("--TPM", type="character"),
    make_option("--LeadCol", type="character"),
    make_option("--GeneReference", type="character"),
    make_option("--CelltypeReference", type="character")
    )

opt <- parse_args(OptionParser(option_list=option.list))

#######################################################
greplany <- function(patterns, v) {
  match <- rep(FALSE, length(v))
  for (pattern in patterns) {
    match <- match | grepl(pattern, v)
  }
  return(match)
}
#######################################################
E2G <- read.table(opt$E2GOverlap, header=TRUE, sep = '\t', stringsAsFactors=FALSE)
geneRef <- read.table(opt$GeneReference, header=FALSE, sep="\t", stringsAsFactors=F)
celltypeRef <- read.table(opt$CelltypeReference, header=T, sep="\t", fill=TRUE, comment.char = "") %>% filter(ConsiderForVariantOverlapping) %>% pull(Celltype)
E2G <- E2G %>% filter(TargetGene %in% geneRef$V4) %>% filter(CellType %in% celltypeRef)

UBQ.genes <- read.table(opt$UBQ_genes, header=FALSE,stringsAsFactors=FALSE)
TPM_df <- read.table(opt$TPM, header=T, sep="\t", stringsAsFactors=F) %>% mutate(Celltype=str_replace(Celltype, "_transcript_per_million", ""))

if (opt$removeNonCoding){
	E2G<- E2G %>% filter(!greplany(c("^LINC","-AS","^MIR","RNU","^LOC"),TargetGene))
}
print(colnames(E2G))

colnames(UBQ.genes) <- "UbiquitouslyExpressedGenes"
UBQ.genes$IsThisGeneUbiquitouslyExpressed <- "TRUE"

E2G <- left_join(E2G,UBQ.genes, by=c("TargetGene" = "UbiquitouslyExpressedGenes"))
E2G <- E2G %>% mutate(IsThisGeneUbiquitouslyExpressed=ifelse(is.na(IsThisGeneUbiquitouslyExpressed),"FALSE","TRUE"))
E2G <- E2G %>% select(-c(TargetGeneEnsembl_ID, normalizedATAC_prom, numCandidateEnhGene,numTSSEnhGene,ubiqExpressed,E2G.Score.qnorm.ignoreTPM)) %>% rename(varChr=chr.1, varStart=start.1, varEnd=end.1, LeadVariant=!!sym(opt$LeadCol)) %>% inner_join(TPM_df %>% filter(TPM>=1), by=c("CellType"="Celltype", "TargetGene"="Gene"))

write.table(E2G %>% distinct(), file.path(opt$outdir, "E2GOverlap_filtered.tsv"), sep="\t", quote=F, row.names=F)

if ("rsid" %in% colnames(E2G)){
    var_info_bed=E2G %>% distinct() %>% select(chr, varEnd, rsid) %>% mutate(start=varEnd-1) %>% select(chr, start, varEnd, rsid) %>% distinct()
} else {
    var_info_bed=E2G %>% distinct() %>% select(chr, varEnd, LeadVariant) %>% mutate(start=varEnd-1) %>% select(chr, start, varEnd, LeadVariant) %>% distinct()
}
write.table(var_info_bed, file.path(opt$outdir, "Credibleset_gene_variant_info_E2G_only.bed"), sep="\t", quote=F, row.names=F, col.names=F)
