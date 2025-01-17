suppressPackageStartupMessages(library(optparse))
option.list <- list(
    make_option("--GenomeAssembly", type="character", default="hg38"),
    make_option("--SampleName", type="character"),
    make_option("--FragmentFile", type="character"),
    make_option("--minTSS", type="numeric", default=7),
    make_option("--minFrags", type="numeric"),
    make_option("--WorkDir", type="character"),
    make_option("--Threads", type="numeric")
)

opt <- parse_args(OptionParser(option_list=option.list))
setwd(opt$WorkDir)

suppressPackageStartupMessages(library(ArchR))
suppressPackageStartupMessages(library(tidyverse))
set.seed(42)

addArchRGenome(opt$GenomeAssembly)
addArchRThreads(threads=opt$Threads)
print(opt$minFrags)
print(opt$minTSS)

arrowfile <- createArrowFiles(
   inputFiles = opt$FragmentFile,
   sampleNames = opt$SampleName,
   minTSS = as.numeric(opt$minTSS), 
   minFrags = as.numeric(opt$minFrags),
   addTileMat = T,
   addGeneScoreMat = T,
   force=T
)



