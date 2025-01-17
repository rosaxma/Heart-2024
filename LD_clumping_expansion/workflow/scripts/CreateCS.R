suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))

option.list <- list(
  make_option("--helperFunctions", type="character"),
  make_option("--fineMappedVariants", type="character"),
  make_option("--genes", type="character"),
  make_option("--promoters", type="character"),
  make_option("--trait", type="character"),
  make_option("--outdir", type="character"),
  make_option("--leadVariantCol", type="character"),
  make_option("--variantCol", type="character"),
  make_option("--chrCol", type="character"),
  make_option("--positionCol", type="character"),
  make_option("--leadSNPPCol", type="character"),
  make_option("--variantPCol", type="character"),
  make_option("--LeadPIPCol", type="character"),
  make_option("--PIPCol", type="character"),
  make_option("--Source", type="character", help="Source of the data"),
  make_option("--zeroIndexed", type="logical", help="Whether the variant position is zero indexed. i.e. is the variant 1 bp off from the record in dbSNP"), 
  make_option("--excludeVariants", type="character", help="Variant to exclude")
)
opt <- parse_args(OptionParser(option_list=option.list))
source(opt$helperFunctions)
##############################################################################
variants=read.table(opt$fineMappedVariants, header=T, stringsAsFactors=F, sep="\t")
##############################################################################
## Load common data
genes <- readBed(opt$genes)
genes$symbol <- unlist(lapply(strsplit(as.character(as.matrix(genes$name)), ";"), "[", 1))
promoters <- readBed(opt$promoters)
promoters$symbol <- unlist(lapply(strsplit(as.character(as.matrix(promoters$name)), ";"), "[", 1))
promoters.gr <- GRangesFromBed(promoters)
############################################
names(variants)[names(variants) == opt$variantCol] <-  "variant"
names(variants)[names(variants) == opt$chrCol] <- "chr"
names(variants)[names(variants) == opt$positionCol] <- "position"

variants=variants %>% subset(., select=which(!duplicated(names(.))))

if (opt$leadSNPPCol != "None"){
    names(variants)[names(variants) == opt$leadSNPPCol] <- "LeadVariant_P_value"
} 
if (opt$variantPCol != "None"){
    names(variants)[names(variants) == opt$variantPCol] <- "variant_P_value"
}
if (opt$leadVariantCol!= "None"){
    names(variants)[names(variants) == opt$leadVariantCol] <- "LeadVariant"
}
print(head(variants,n=2))
if (opt$leadSNPPCol == "None"){
    if (opt$variantPCol != "None" & opt$leadVariantCol!= "None") {
      leadSNP_p=variants %>% filter(variant==LeadVariant) %>% select(LeadVariant, variant_P_value) %>% distinct() %>% mutate(LeadVariant_P_value=variant_P_value) %>% select(-variant_P_value)
      variants <- variants %>% left_join(leadSNP_p, by=c("LeadVariant"))
    } else if (opt$variantPCol != "None" & opt$PIPCol == "None" ){
        leadSNP_p =variants %>% select(CredibleSet, variant, variant_P_value) %>% group_by(CredibleSet) %>% slice_min(variant_P_value, n=1) %>% ungroup() %>% as.data.frame() %>% distinct() %>% mutate(LeadVariant=variant)%>% mutate(LeadVariant_P_value=variant_P_value) %>% select(CredibleSet, LeadVariant,LeadVariant_P_value)
        variants <- variants %>% left_join(leadSNP_p, by=c("CredibleSet"))
    } else {
      variants <- variants %>% mutate(variant_P_value=NA)%>% mutate(LeadVariant_P_value=NA)
    }
}

print(head(variants,n=2))
if (opt$LeadPIPCol != "None"){
	names(variants)[names(variants) == opt$LeadPIPCol] <- "LeadSNP_PosteriorProb"
} 
if (opt$PIPCol != "None"){
	names(variants)[names(variants) == opt$PIPCol] <- "PosteriorProb"
}
print(head(variants,n=2))
if (opt$LeadPIPCol == "None") {
    if (opt$PIPCol != "None" & opt$leadVariantCol!= "None"){
      leadSNP_pip=variants %>% filter(variant==LeadVariant) %>% select(LeadVariant, PosteriorProb) %>% distinct() %>% mutate(LeadSNP_PosteriorProb=PosteriorProb) %>% select(-PosteriorProb)
      variants <- variants %>% left_join(leadSNP_pip, by=c("LeadVariant"))
      } else if (opt$leadVariantCol== "None"){
      leadSNP_pip=variants %>% select(CredibleSet, variant, PosteriorProb) %>% group_by(CredibleSet) %>% distinct() %>% slice_max(PosteriorProb, n=1) %>% ungroup() %>% as.data.frame() %>% distinct() %>% mutate(LeadVariant=variant) %>% mutate(LeadSNP_PosteriorProb=PosteriorProb) %>% select(CredibleSet, LeadVariant, LeadSNP_PosteriorProb)
      variants <- variants %>% left_join(leadSNP_pip, by=c("CredibleSet"))
      } else {
	  variants <- variants %>% mutate(PosteriorProb=NA)%>% mutate(LeadSNP_PosteriorProb=NA)
      }
}

variants=variants %>% subset(., select=which(!duplicated(names(.))))

if (opt$excludeVariants!="None"){
    excludeVariants=opt$excludeVariants %>% strsplit(split=",") %>% unlist()    
    variants <- variants %>% filter(!(LeadVariant %in% excludeVariants))
}

print(head(variants,n=2))

if (opt$zeroIndexed){
    variants <- variants %>% mutate(
    CredibleSet=paste0(opt$trait, "-", LeadVariant),
    Trait=opt$trait,
    chr=ifelse(str_starts(chr, "Chr|chr"), chr, paste0("chr", chr)),
    position=as.numeric(position)+1) %>% group_by(LeadVariant) %>% mutate(LeadVariantPos_tmp=ifelse(variant==LeadVariant, position, NA)) %>% mutate(LeadVariantPos=max(as.numeric(LeadVariantPos_tmp), na.rm = TRUE)) %>% ungroup() %>% select(-LeadVariantPos_tmp)%>% filter(!is.infinite(LeadVariantPos))
} else {
    variants <- variants %>% mutate(
    CredibleSet=paste0(opt$trait,"-", LeadVariant),
    Trait=opt$trait,
    chr=ifelse(str_starts(chr, "Chr|chr"), chr, paste0("chr", chr)),
    position=as.numeric(position)) %>% group_by(LeadVariant) %>% mutate(LeadVariantPos_tmp=ifelse(variant==LeadVariant, position, NA)) %>% mutate(LeadVariantPos=max(as.numeric(LeadVariantPos_tmp), na.rm = TRUE)) %>% ungroup() %>% select(-LeadVariantPos_tmp) %>% filter(!is.infinite(LeadVariantPos))
}


variant.gr <- with(variants, GRangesFromBed(data.frame(chr=chr, start=position, end=position)))
variants2 <- annotateVariants(variants, variant.gr, promoters.gr, include.names=T)

variants <- variants2 %>% mutate(start=position, end=position, LocusID=LeadVariant) %>%
  dplyr:::select(
    chr,
    start, 
    end,
    position,
    LocusID, 
    LeadVariant_P_value,
    variant,
    variant_P_value,
    LeadSNP_PosteriorProb,
    PosteriorProb,
    CredibleSet,
    Trait,
    Coding,
    CodingVariantGene,
    SpliceSite,
    SpliceSiteVariantGene,
    Promoter,
    PromoterVariantGene
    ) 
##############################################################################
write.table(variants, file=file.path(opt$outdir,"variant.list.txt"), row.names=F, col.names=T, sep='\t', quote=F)
print("debug2")
variants_pip <- variants2 %>% mutate(start=position, end=position, LocusID=LeadVariant) %>%
    dplyr:::select(
      chr,
      start,
      end,
      position,
      LocusID, 
      LeadVariant_P_value, 
      variant,
      variant_P_value,
      LeadSNP_PosteriorProb,
      PosteriorProb,
      CredibleSet,
      Trait,
      Coding,
      CodingVariantGene,
      SpliceSite,
      SpliceSiteVariantGene,
      Promoter,
      PromoterVariantGene
      ) 
if (opt$PIPCol == "None"){
    variants_pip=variants_pip %>% mutate(LeadSNP_PosteriorProb=0, PosteriorProb=0)
} 
print(head(variants_pip %>% as.data.frame()))
all.cs <- makeCredibleSets(variants_pip, genes, Source=opt$Source,include.names=TRUE)
write.table(all.cs %>% distinct(), file=file.path(opt$outdir, "all.cs.txt"),  row.names=F, col.names=T, sep='\t', quote=F)
##############################################################################
tmp <- writeVariantBed(variants, file.path(opt$outdir, "Variants.bed"))
if (opt$PIPCol != "None"){
  variants %>% dplyr:::select(chr, start, end, PosteriorProb) %>% writeBed(file=file.path(opt$outdir,"Variants.bedgraph"))
}
tmp <- writeCredibleSetBed(all.cs, variants, file=file.path(opt$outdir, "CredibleSets.bed"))
