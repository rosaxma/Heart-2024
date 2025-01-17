suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrastr))
option.list <- list(
    make_option("--outdir", type="character"),
    make_option("--FragmentFile",type="character"),
    make_option("--sample",type="character"),  
    make_option("--AmuletMultiplets", type="character"),
    make_option("--ArchRMultiplets", type="character")

)
opt <- parse_args(OptionParser(option_list=option.list))
frag_table <- read.table(opt$FragmentFile, sep="\t", header=FALSE, stringsAsFactor=FALSE) %>% distinct()
amulet_multiplets_table <- read.table(opt$AmuletMultiplets, sep="\t", header=FALSE, stringsAsFactor=FALSE) %>% distinct() %>% mutate(multiplets=TRUE)
colnames(amulet_multiplets_table) <- c("cbc", "amulet_multiplets")

colnames(frag_table)=c("chr", "start", "end", "cbc", "placeholder")

if (file.info(opt$ArchRMultiplets)$size != 0) {
archr_multiplets_table <- read.table(opt$ArchRMultiplets, sep="\t", header=FALSE, stringsAsFactor=FALSE) %>% distinct() %>% mutate(multiplets=TRUE)
colnames(archr_multiplets_table) <- c("cbc", "archr_multiplets")

count_table <- frag_table %>% distinct() %>% count(cbc) %>% rename(nfrag=n) %>% mutate(rank=rank(-nfrag)) %>% left_join(amulet_multiplets_table) %>% mutate(amulet_multiplets=ifelse(is.na(amulet_multiplets), "FALSE", "TRUE")) %>% left_join(archr_multiplets_table) %>% mutate(archr_multiplets=ifelse(is.na(archr_multiplets), "FALSE", "TRUE")) %>% rowwise () %>% mutate(multiplets=ifelse((amulet_multiplets=="TRUE" || archr_multiplets=="TRUE"), "TRUE", "FALSE"))

} else {
count_table <- frag_table %>% count(cbc) %>% rename(nfrag=n) %>% mutate(rank=rank(-nfrag)) %>% left_join(amulet_multiplets_table) %>% mutate(amulet_multiplets=ifelse(is.na(amulet_multiplets), "FALSE", "TRUE")) %>% mutate(archr_multiplets="FALSE") %>% rowwise () %>% mutate(multiplets=ifelse((amulet_multiplets=="TRUE" || archr_multiplets=="TRUE"), "TRUE", "FALSE"))
}


pdf(file.path(opt$outdir,paste0(opt$sample,"_Kneeplot_complete.pdf")))
boo_color=scale_color_manual(values = c("FALSE"= "#999999", "TRUE" = "#ff0000"))

p1 <- ggplot(count_table, aes(x=rank, y=nfrag, color=multiplets)) + rasterize(geom_point(),dpi=300) + scale_x_log10(limits = c(1e0, 1e6),n.breaks=14) + scale_y_log10(limits = c(1e0, 1e6),n.breaks=7) + boo_color + ggtitle(paste0(opt$sample, "_All multiplets counts: ", nrow(count_table)))

p2 <- ggplot(count_table %>% filter(amulet_multiplets=="TRUE"), aes(x=rank, y=nfrag, color=amulet_multiplets)) + rasterize(geom_point(),dpi=300) + scale_x_log10(limits = c(1e0, 1e6),n.breaks=14) + scale_y_log10(limits = c(1e0, 1e6),n.breaks=7) + boo_color + ggtitle(paste0("Amulet multiplets counts: ", nrow(count_table %>% filter(amulet_multiplets=="TRUE"))))

p3 <- ggplot(count_table %>% filter(archr_multiplets=="TRUE"), aes(x=rank, y=nfrag, color=archr_multiplets)) + rasterize(geom_point(),dpi=300) + scale_x_log10(limits = c(1e0, 1e6),n.breaks=14) + scale_y_log10(limits = c(1e0, 1e6),n.breaks=7) + boo_color+ ggtitle(paste0("Archr doublets  counts: ", nrow(count_table %>% filter(archr_multiplets=="TRUE"))))

p4 <- ggplot(count_table %>% filter((archr_multiplets=="TRUE" && amulet_multiplets=="TRUE")), aes(x=rank, y=nfrag, color=multiplets)) + rasterize(geom_point(),dpi=300) + scale_x_log10(limits = c(1e0, 1e6),n.breaks=14) + scale_y_log10(limits = c(1e0, 1e6),n.breaks=7) + boo_color+ ggtitle(paste0("Intersected multiplets counts: ", nrow(count_table %>% filter((archr_multiplets=="TRUE" && amulet_multiplets=="TRUE")))))

p5 <- ggplot(count_table %>% filter((archr_multiplets=="FALSE" && amulet_multiplets=="FALSE")), aes(x=rank, y=nfrag, color=multiplets)) + rasterize(geom_point(),dpi=300) + scale_x_log10(limits = c(1e0, 1e6),n.breaks=14) + scale_y_log10(limits = c(1e0, 1e6),n.breaks=7) + boo_color+ ggtitle(paste0("Non-multiplets counts: ", nrow(count_table %>% filter((archr_multiplets=="FALSE" && amulet_multiplets=="FALSE")))))

print(p1)
print(p2)
print(p3)
print(p4)
print(p5)
dev.off()

pdf(file.path(opt$outdir,paste0(opt$sample,"_Kneeplot.pdf")))
print(p1)
dev.off()