# Analysis code for *Molecular convergence of risk variants for congenital heart defects leveraging a regulatory map of the human fetal heart*

## envrionments
  ### environment for snakemake runs:
      envs/snakemake.yml 
  ### enivronemnt for running R scripts:
      envs/R.yml

## scRNA-seq alignment and background removal
  ### RNA_alignment_snakemake
      The repo
  ### CellBender_Snakemake
  
## scATAC-seq alignment
  ### scATAC_pipeline_lite 

## scATAC-seq quality control
  ### multiomic_ATAC_QC_Snakemake

## scRNA-seq quality control and normalization
  ### RNA_QC_Snakemake
  ### Normalization_Snakemake
  
## Doublet removal
  ### Doublet_removal_Snakemake

## Clustering and cell type annotation
  ### evaluation_Snakemake
  ### clustering_Snakemake
  ### annotation_Snakemake

## Defining gene programs with consensus non-negative matrix factorization
  ### cNMF_run_and_analysis/Snakefile_cNMF
  ### cNMF_run_and_analysis/Snakefile_cNMF_analysis

## Identifying transcription factors regulating gene programs 
  ### ChromVar_Snakemake
  ### get_motif_enrichment

## Testing disease genes for preferential expression in particular cell types
  ### Preferential expression of CHD genes

## Assessing cell type heritability enrichment with stratified LD score regression
  ### LDSC_Snakemake 

## Preparing variants for overlapping with scE2G-predicted enhancers
  ### LD_clumping_expansion

## Linking GWAS variants to enhancers, target genes, cell types, and gene programs
  ### V2G_v2

## Assessing enrichment of genes associated with valve traits in VIC program 27 
  ### Assessing enrichment of genes associated with valve traits in VIC program 27






