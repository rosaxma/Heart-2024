suppressPackageStartupMessages(library(optparse))
option.list <- list(
    make_option("--GenomeAssembly", type="character", default="hg38"),
    make_option("--SampleName", type="character"),
    make_option("--FragmentFile", type="character"),
    make_option("--minTSS", type="numeric", default=7),
    make_option("--minFrags", type="numeric"),
    make_option("--WorkDir", type="character"),
    make_option("--Threads", type="numeric"),
    make_option("--ATAC_barcodes", type="character"),
    make_option("--nuc_signal_table", type="character")
)

opt <- parse_args(OptionParser(option_list=option.list))
setwd(opt$WorkDir)

suppressPackageStartupMessages(library(ArchR))
suppressPackageStartupMessages(library(tidyverse))
set.seed(42)


NucSigFile = read.table(opt$nuc_signal_table, header=T, sep="\t", stringsAsFactors=F)
validcbc=NucSigFile %>% filter(keep) %>% pull(cbc)
invalidcbc=NucSigFile %>% filter(!keep) %>% pull(cbc)
#######################################################################
addArchRGenome(opt$GenomeAssembly)
addArchRThreads(threads=opt$Threads)
print(opt$minFrags)
print(opt$minTSS)
cbc <- read.table(opt$ATAC_barcodes, sep="\t", header=T, stringsAsFactors=F)
arrowfile <- createArrowFiles(
   inputFiles = opt$FragmentFile,
   sampleNames = opt$SampleName,
   minTSS = as.numeric(opt$minTSS), 
   minFrags = as.numeric(opt$minFrags),
   addTileMat = T,
   addGeneScoreMat = T,
   force=T,
   validBarcodes=intersect(cbc$ATAC_cbc,validcbc)
)

if (length(intersect(cbc$ATAC_cbc,invalidcbc)) != 0){
   dir.create(paste0(opt$WorkDir, "_invalid"), showWarnings = FALSE)
   setwd(paste0(opt$WorkDir, "_invalid"))
   rrowfile <- createArrowFiles(
   inputFiles = opt$FragmentFile,
   sampleNames = paste0(opt$SampleName, "_invalid"),
   validBarcodes=intersect(cbc$ATAC_cbc,invalidcbc),
   minTSS = as.numeric(opt$minTSS), 
   minFrags = as.numeric(opt$minFrags),
   addTileMat = T,
   addGeneScoreMat = T,
   force = T
)
}



