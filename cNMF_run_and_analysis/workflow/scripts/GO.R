library(clusterProfiler)
library(msigdbr)
library(optparse)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
option.list = list(
    make_option("--betaScoreTable", type="character"),
    make_option("--topGeneTable", type="character"),
    make_option("--GeneIDdict", type="character"),
    make_option("--outdir", type="character"), 
    make_option("--numK", type="character"),
    make_option("--topGenes2look", type="character")
)
opt = parse_args(OptionParser(option_list=option.list))


dict = readRDS(opt$GeneIDdict)
beta_table = read.table(opt$betaScoreTable,header=TRUE, sep="\t",stringsAsFactors=FALSE, row.names=1)
topGeneTable = read.table(opt$topGeneTable,header=TRUE, sep="\t",stringsAsFactors=FALSE)
# MSigDB: H; Hallmark gene sets: coherently expressed signatures derived by aggregating many MSigDB gene sets to represent well-defined biological states or processes
h_df = msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, entrez_gene)

dict_vec = dict$entrezID
names(dict_vec) = dict$Gene

# MSigDB: C8; curated from cluster markers identified in single cell sequencing studies of human tissue
c8_df = msigdbr(species = "Homo sapiens", category = "C8") %>% dplyr::select(gs_name, entrez_gene)
universe = as.character(dict$entrezID)
topGene_colnames=grep("genes", colnames(topGeneTable),value=TRUE)

# MSigDB: C5; curated from cluster markers identified in single cell sequencing studies of human tissue
c5_df = msigdbr(species = "Homo sapiens", category = "C5",subcategory = "BP") %>% dplyr::select(gs_name, entrez_gene)

results = as.data.frame(matrix(nrow=2,ncol=5))
colnames(results) =c("HallmarkGeneSets", "GO", "Hallmark_AverageUniquePathwaysPerGEP","GO_AverageUniqueGOPerGEP", "K")
rownames(results) = c("Hypergeometric", "GSEA")

# hypergeometric test
enricher_h_list =  vector(mode = "list", length = length(topGene_colnames))
enricher_c8_list =  vector(mode = "list", length = length(topGene_colnames))
enricher_c5_list =  vector(mode = "list", length = length(topGene_colnames))

is.not.null <- function(x) !is.null(x)

for (topic in topGene_colnames) {
    print(topic)
    gene_set = na.omit(topGeneTable[[topic]])
    gene_set = gene_set[1:as.numeric(opt$topGenes2look)]
    gene_set_entrezID = dict %>% filter(Gene %in% gene_set) %>% pull(entrezID)

    enricher_output = enricher(gene=gene_set_entrezID, universe=universe, TERM2GENE=h_df, pvalueCutoff = 1,  pAdjustMethod = "BH")
    if (is.not.null(enricher_output)){
        h_df_results = enricher_output@result %>% mutate(GEP=topic)
        enricher_h_list[[topic]] = h_df_results
    }

    enricher_output = enricher(gene=gene_set_entrezID, universe=universe, TERM2GENE=c8_df, pvalueCutoff = 1,  pAdjustMethod = "BH")
    if (is.not.null(enricher_output)){
        c8_df_results = enricher_output@result %>% mutate(GEP=topic)
        enricher_c8_list[[topic]] = c8_df_results  
    }

    enricher_output = enricher(gene=gene_set_entrezID, universe=universe, TERM2GENE=c5_df, pvalueCutoff = 1,  pAdjustMethod = "BH")
    if (is.not.null(enricher_output)){
        c5_df_results = enricher_output@result %>% mutate(GEP=topic)
        enricher_c5_list[[topic]] = c5_df_results
    }
}

enricher_h_df= do.call(what = rbind, args = enricher_h_list)
enricher_c8_df= do.call(what = rbind, args = enricher_c8_list)
enricher_c5_df= do.call(what = rbind, args = enricher_c5_list)

# GeneSetEnrichment
beta_colnames = grep("GEP", colnames(beta_table),value=TRUE)
GSEA_h_list =  vector(mode = "list", length = length(beta_colnames))
GSEA_c8_list =  vector(mode = "list", length = length(beta_colnames))
GSEA_c5_list =  vector(mode = "list", length = length(beta_colnames))

for (topic in beta_colnames) {
    print(topic)
    ranked_gene_df  = beta_table %>% select(!!sym(topic)) %>% rownames_to_column(var="Gene") %>% arrange(desc(!!sym(topic))) %>% left_join(dict, by="Gene") %>% select(entrezID,!!sym(topic))
    geneList = ranked_gene_df %>% pull(!!sym(topic))
    names(geneList) = ranked_gene_df %>% pull(entrezID)
    print(length(geneList))

    GSEA_output = GSEA(geneList=geneList, pvalueCutoff = 1, pAdjustMethod = "BH", TERM2GENE=h_df, seed=FALSE)
    if (is.not.null(GSEA_output)){
        h_df_results = GSEA_output@result %>% mutate(GEP=topic)
        GSEA_h_list[[topic]] = h_df_results
    }

    GSEA_output= GSEA(geneList=geneList, pvalueCutoff = 1, pAdjustMethod = "BH", TERM2GENE=c8_df, seed=FALSE)
    if (is.not.null(GSEA_output)){
        c8_df_results = GSEA_output@result %>% mutate(GEP=topic)
        GSEA_c8_list[[topic]] = c8_df_results
    }

    GSEA_output= GSEA(geneList=geneList, pvalueCutoff = 1, pAdjustMethod = "BH", TERM2GENE=c5_df, seed=FALSE)
    if (is.not.null(GSEA_output)){
        c5_df_results = GSEA_output@result %>% mutate(GEP=topic)
        GSEA_c5_list[[topic]] = c5_df_results
    }
}


enricher_h_df <- enricher_h_df %>% filter(p.adjust < 0.05) %>% mutate(geneID_tmp=str_split(geneID, "/")) %>% rowwise() %>% mutate(geneSymbol_tmp=list(names(dict_vec[dict_vec %in% geneID_tmp]))) %>% mutate(geneSymbol=paste(geneSymbol_tmp, collapse="/")) %>% select(-c(geneID_tmp, geneSymbol_tmp)) %>% relocate(geneSymbol, .after=geneID)

enricher_c8_df<- enricher_c8_df %>% filter(p.adjust < 0.05) %>% mutate(geneID_tmp=str_split(geneID, "/")) %>% rowwise() %>% mutate(geneSymbol_tmp=list(names(dict_vec[dict_vec %in% geneID_tmp]))) %>% mutate(geneSymbol=paste(geneSymbol_tmp, collapse="/")) %>% select(-c(geneID_tmp, geneSymbol_tmp)) %>% relocate(geneSymbol, .after=geneID)

enricher_c5_df<- enricher_c5_df %>% filter(p.adjust < 0.05) %>% mutate(geneID_tmp=str_split(geneID, "/")) %>% rowwise() %>% mutate(geneSymbol_tmp=list(names(dict_vec[dict_vec %in% geneID_tmp]))) %>% mutate(geneSymbol=paste(geneSymbol_tmp, collapse="/")) %>% select(-c(geneID_tmp, geneSymbol_tmp)) %>% relocate(geneSymbol, .after=geneID)


GSEA_h_df= do.call(what = rbind, args = GSEA_h_list) %>% filter(NES>0) %>% filter(p.adjust < 0.05) %>% mutate(geneID_tmp=str_split(core_enrichment, "/")) %>% rowwise() %>% mutate(geneSymbol_tmp=list(names(dict_vec[dict_vec %in% geneID_tmp]))) %>% mutate(geneSymbol=paste(geneSymbol_tmp, collapse="/")) %>% select(-c(geneID_tmp, geneSymbol_tmp)) %>% relocate(geneSymbol, .after=core_enrichment)

GSEA_c8_df= do.call(what = rbind, args = GSEA_c8_list) %>% filter(NES>0)%>% filter(p.adjust < 0.05) %>% mutate(geneID_tmp=str_split(core_enrichment, "/")) %>% rowwise() %>% mutate(geneSymbol_tmp=list(names(dict_vec[dict_vec %in% geneID_tmp]))) %>% mutate(geneSymbol=paste(geneSymbol_tmp, collapse="/")) %>% select(-c(geneID_tmp, geneSymbol_tmp)) %>% relocate(geneSymbol, .after=core_enrichment)


GSEA_c5_df= do.call(what = rbind, args = GSEA_c5_list) %>% filter(NES>0)%>% filter(p.adjust < 0.05) %>% mutate(geneID_tmp=str_split(core_enrichment, "/")) %>% rowwise() %>% mutate(geneSymbol_tmp=list(names(dict_vec[dict_vec %in% geneID_tmp]))) %>% mutate(geneSymbol=paste(geneSymbol_tmp, collapse="/")) %>% select(-c(geneID_tmp, geneSymbol_tmp)) %>% relocate(geneSymbol, .after=core_enrichment)


results["Hypergeometric", "HallmarkGeneSets"] = length(unique(enricher_h_df$ID))
results["Hypergeometric", "GO"] = length(unique(enricher_c5_df$ID))
results["Hypergeometric", "Hallmark_AverageUniquePathwaysPerGEP"] = enricher_h_df %>% select(ID, GEP) %>% group_by(GEP) %>% summarize(n_unique=n_distinct(ID)) %>% ungroup() %>% summarize(mean(n_unique))
results["Hypergeometric", "GO_AverageUniqueGOPerGEP"] = enricher_c5_df %>% select(ID, GEP) %>% group_by(GEP) %>% summarize(n_unique=n_distinct(ID)) %>% ungroup() %>% summarize(mean(n_unique))

results["GSEA", "HallmarkGeneSets"] = length(unique(GSEA_h_df$ID))
results["GSEA", "GO"] = length(unique(GSEA_c5_df$ID))
results["GSEA", "Hallmark_AverageUniquePathwaysPerGEP"] = GSEA_h_df %>% select(ID, GEP) %>% group_by(GEP) %>% summarize(n_unique=n_distinct(ID)) %>% ungroup() %>% summarize(mean(n_unique))
results["GSEA", "GO_AverageUniqueGOPerGEP"] = GSEA_c5_df %>% select(ID, GEP) %>% group_by(GEP) %>% summarize(n_unique=n_distinct(ID)) %>% ungroup() %>% summarize(mean(n_unique))
results$K = rep(opt$numK, 2)

write.table(enricher_h_df, file.path(opt$outdir, paste0("K_", opt$numK, "_Hypergeometric_HallmarkGeneSet.txt")), sep="\t", quote=F, row.names=F)
write.table(enricher_c8_df,file.path(opt$outdir, paste0("K_", opt$numK, "_Hypergeometric_ClusterMarkers.txt")), sep="\t", quote=F, row.names=F )
write.table(enricher_c5_df,file.path(opt$outdir, paste0("K_", opt$numK, "_Hypergeometric_GO.txt")), sep="\t", quote=F, row.names=F )
write.table(GSEA_h_df, file.path(opt$outdir, paste0("K_", opt$numK, "_GSEA_HallmarkGeneSet.txt")), sep="\t", quote=F, row.names=F)
write.table(GSEA_c8_df, file.path(opt$outdir, paste0("K_", opt$numK, "_GSEA_ClusterMarkers.txt")), sep="\t", quote=F, row.names=F)
write.table(GSEA_c5_df, file.path(opt$outdir, paste0("K_", opt$numK, "_GSEA_GO.txt")), sep="\t", quote=F, row.names=F)
write.table(results, file.path(opt$outdir, paste0("K_", opt$numK, "_uniquePathwayID_count.txt")), sep="\t", quote=F, row.names=T)



