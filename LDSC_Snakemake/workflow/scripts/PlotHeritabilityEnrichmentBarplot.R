# Jesse Engreitz
# December 2019
##Rscript to run automatically after LDSC runs
suppressPackageStartupMessages(library("optparse"))

## To do -- make this more customizable
option.list <- list(
  make_option(c("-o", "--output"), default="/oak/stanford/groups/engreitz/Users/rosaxma/2109_LDSC_FT/FT05007010f3/AvgHic.ABC0.02.minus150/ldsc/dz/meta_asd.results.pdf", type="character", help="Output PDF file"),
  make_option(c("-i", "--input"), default="/oak/stanford/groups/engreitz/Users/rosaxma/2109_LDSC_FT/FT05007010f3/AvgHic.ABC0.02.minus150/ldsc/dz/meta_asd.results",type="character", help=".results.txt file output by LDSC"),
  make_option(c("-c", "--control"), type="character", help=".results.txt file for control set, e.g. for DNase-peaks only comparison"),
  make_option("--params", default="/oak/stanford/groups/engreitz/Users/rosaxma/2109_LDSC_FT/misc/CellParams.forABCGWAS.FT005007010f3.txt", type="character", help="Cell type params file"),
  make_option("--include-error-bars", default=TRUE, type="logical", help="skip error bars"),
  make_option("--ID", default=TRUE, type="character", help="disease ID")
)

opt <- parse_args(OptionParser(option_list=option.list))
library(ggplot2)
params <- read.delim(opt$params, stringsAsFactors=F)

readLdsc <- function(file) {
  he <- read.delim(file, stringsAsFactors=F)
  toplot <- subset(he, Category %in% paste0(params$`cell_type`,"L2_0"))
  toplot$Category <- gsub("L2_0","",toplot$Category)
  toplot <- toplot[order(toplot$Enrichment),]
  toplot$Category <- factor(toplot$Category, levels=toplot$Category)
  toplot$Significant <- toplot$Enrichment_p < 0.05
  return(toplot)
}

toplot <- readLdsc(opt$input)

if (!is.null(opt$control)) {
  print("Loading control")
  control <- readLdsc(opt$control)
  control$Enrichment_std_error <- 0
  control$Significant <- 'Control'
  toplot$Significant <- as.character(toplot$Significant)
  toplot <- rbind(toplot, control)
}


getEnrichmentPlot <- function(toplot, ylabel, ID) {
    #toplot$Significant <- ordered(as.character(as.matrix(toplot$Significant)), levels=c("gray","black","red"))
    #toplot <- toplot[toplot$Significant,]
    p <- ggplot(toplot, aes(x=Category)) + geom_point(aes(x=Category, y=Enrichment, colour=Significant), size=2)
    lowerbound = min(toplot$Enrichment-10)    
    
    if ("Enrichment_std_error" %in% colnames(toplot) & opt$`include-error-bars`) {
        p <- p + geom_errorbar(aes(ymin=Enrichment-Enrichment_std_error*1.96, ymax=Enrichment+Enrichment_std_error*1.96), width=.2, position=position_dodge(width=0.5))+geom_linerange(aes(ymin = Enrichment, ymax = Enrichment+Enrichment_std_error*1.96))
        }
    if ('TRUE' %in% toplot$Significant){
	      p <- p + geom_point(data=subset(toplot, Significant == 'TRUE'), aes(x=Category, y=Enrichment, colour=Significant, size=1.5), show.legend = FALSE)
    }
    p <- p + theme_classic() 
    p <- p + ylab(ylabel)
    p <- p + ylim(lowerbound,NA)
    p <- p + theme(
      axis.text.y = element_text(size=5), 
      axis.ticks.y = element_blank(), 
      axis.text.x = element_text(size=25), 
      axis.title=element_text(size=25)
      )
    p <- p + scale_color_manual(values=c('TRUE'='red','FALSE'='black'))
    p <- p + coord_flip()
    p <- p + geom_hline(yintercept=1, linetype='dashed', col='gray') + ggtitle(ID)
    return(p)
}
getCoefficientPlot <-function(toplot,ylabel,ID){
    print("here3.5")
    p <- ggplot(toplot, aes(x=Category)) + geom_point(aes(x=Category, y=Coefficient, colour= "black"), size=2,show.legend = FALSE) 
    lowerbound = min(toplot$Coefficient-4)
    if ("Coefficient_std_error" %in% colnames(toplot) & opt$`include-error-bars`) {
	    p <- p + geom_errorbar(aes(ymin=Coefficient-Coefficient_std_error*1.96, ymax=Coefficient+Coefficient_std_error*1.96), width=.1, position=position_dodge(width=0.5))+geom_linerange(aes(ymin=Coefficient, ymax = Coefficient+Coefficient_std_error*1.96))
      }
    p <- p + theme_classic()
    p <- p + ylab(ylabel)
    p <- p + theme(
      axis.text.y = element_text(size=5),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size=25),
      axis.title=element_text(size=25))
    p <- p + coord_flip()
    p <- p + geom_hline(yintercept=0, linetype='dashed', col='gray') + ggtitle(ID)
    return(p)
}
getCoefficientZScorePlot <-function(toplot,ylabel, ID){
    p <- ggplot(toplot, aes(x=Category)) + geom_point(aes(x=Category, y=Coefficient_z.score, colour= "black"), size=2,show.legend = FALSE)
    lowerbound = min(toplot$Enrichment-10)
    p <- p + theme_classic()
    p <- p + ylab(ylabel)
    p <- p + theme(
      axis.text.y = element_text(size=5),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size=25),
      axis.title=element_text(size=25))
    p <- p + coord_flip()
    p <- p + geom_hline(yintercept=0, linetype='dashed', col='gray') + ggtitle(ID)
    return(p)
}


pdf(file=opt$output, width=10, height=length(unique(toplot$Category))*0.2)
p <- getEnrichmentPlot(toplot, "Enrichment\n(% heritability / % SNPs)", opt$ID)
print(p)
p <- getCoefficientZScorePlot(toplot, "Coefficient Z score", opt$ID)
print(p)
p <-getCoefficientPlot(toplot,"Coefficient score", opt$ID)
print(p)
dev.off()


