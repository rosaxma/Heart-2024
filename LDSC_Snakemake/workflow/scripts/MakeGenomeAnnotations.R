##################################################################
suppressPackageStartupMessages(library("optparse"))

## parse options
option.list <- list(
  make_option("--input", type="character"),
  make_option("--output", type="character"),
  make_option("--codeDir", type="character"),
  make_option("--sizes", type="character"),
  make_option("--baseline", type="character")
)
opt <- parse_args(OptionParser(option_list=option.list))

suppressPackageStartupMessages(source(paste0(opt$codeDir, "/src/libs/JuicerUtilities.R")))
suppressPackageStartupMessages(library(GenomicRanges))


cat("Reading input...\n")
input <- unfactor(read.delim(opt$input))

cat("Reading sizes...\n")
sizes <- read.delim(opt$sizes, header=F); 
colnames(sizes) <- c("chr","length")
seqinfo <- Seqinfo(as.character(as.matrix(sizes$chr)), sizes$length)


save.image(paste0(opt$output,"/MakeGenomeAnnotations.params.RData"))

## Write annotation files for heritability partitioning using original, potentially overlapping categories
for (chr in 1:22) {
  cat(paste0("Starting chromosome ", chr, " ...\n"))
  annot <- read.delim(gzfile(paste0(opt$baseline, ".",chr,".annot.gz")), sep='')
  #print(head(annot))
  for (i in 1:nrow(input)) {
    tryCatch({
      bed <- subset(readBed(input$file[i]), chr %in% seqnames(seqinfo))
      bed$start <- bed$start+1

      curr<-GRanges(bed[c('chr','start','end')])
      #curr <- GRangesFromBed(bed, seqinfo=seqinfo)
      ixn <- findOverlaps(GRangesFromBed(with(annot, data.frame(chr=paste0("chr",chr), start=BP-1, end=BP))), curr)
      annot[,input$name[i]] <- 0
      annot[unique(queryHits(ixn)),input$name[i]] <- 1
    }, error = function(e) {
      print(e)
      cat(paste0("Skipping ", input$name[i], "\n"))
    })
  }
  gz <- gzfile(paste0(opt$output, "/overlapping.",chr,".annot.gz"), 'w')
  write.tab(annot, file=gz)
  close(gz)
}
