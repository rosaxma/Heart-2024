proj <- loadArchRProject(path="/oak/stanford/groups/engreitz/Users/rosaxma/230828_heart_map/chromVar/240507/outdir/core_conduction/ArchR/ArchRProject")
getAvailableMatrices(ArchRProj=proj)
motif_mtx <- getMatrixFromProject(ArchRProj=proj, useMatrix="MotifMatrix")
motif_z_mtx <- assays(motif_mtx)$z
tf_regulators <- read.table("/oak/stanford/groups/engreitz/Users/rosaxma/230828_heart_map/chromVar/240507/outdir/core_conduction/analysis/TF/core_conduction_motif_TF_expression_correlation.tsv", header=T, sep="\t", stringsAsFactors=F)
tf <- tf_regulators %>% filter(TFRegulator=="YES") %>% pull(MotifMatrix_name)
motif_tf_z_mtx=motif_z_mtx[tf, ]
info <- read.table("/oak/stanford/groups/engreitz/Users/rosaxma/230828_heart_map/figures/figure1/tables/all_sample_info_table.tsv", header=T, sep="\t", stringsAsFactors=F)
info <- info %>% mutate(atac_cbc=paste0(sample_orig, "_", ATAC_barcode))
info <- info %>% mutate(atac_cbc=paste0("core_conduction", "#", atac_cbc))
mtx_cbc=colnames(motif_tf_z_mtx)
info <- info %>% filter(atac_cbc %in% mtx_cbc)
motif_tf_z_df <- as.data.frame(as.matrix(motif_tf_z_mtx)) %>% rownames_to_column(var="motif") %>% pivot_longer(cols=starts_with("core"), names_to="atac_cbc", values_to="chromvar_z_score")
motif_tf_z_df <- motif_tf_z_df %>% left_join(info %>% select(atac_cbc, celltype)) %>% select(-atac_cbc)
motif_tf_z_mtx_2 <- motif_tf_z_df %>% group_by(motif, celltype) %>% summarize(avg_chromvar_z_score=mean(chromvar_z_score))
motif_tf_z_mtx_2 <- motif_tf_z_mtx_2 %>% pivot_wider(names_from=celltype, values_from=avg_chromvar_z_score)
library(ComplexHeatmap)

