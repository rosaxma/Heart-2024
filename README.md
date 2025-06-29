# Analysis code for *Molecular convergence of risk variants for congenital heart defects leveraging a regulatory map of the human fetal heart*
The snakemake pipelines were designed to run on [Sherlock](https://www.sherlock.stanford.edu), an HPC platform at Stanford University. The config files are located in \[dir\]/config. The Methods section of the corresponding [manuscript](https://www.medrxiv.org/content/10.1101/2024.11.20.24317557v2.full) provides detailed explanations of each step. 
## Envrionments
  ### environment for snakemake runs:
  envs/snakemake.yml 
  ### enivronemnt for running R scripts:
  envs/R.yml

## scRNA-seq alignment and background removal and Sample demultiplexing
  ### RNA_alignment_snakemake
  The code for aligning the scRNA portion of single cell multiomic sequencing to a reference assembly using [Starsolo](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md) and demultiplexing samples using [Souporcell](https://github.com/wheaton5/souporcell). The input sample table requires information about the location of FASTQ files, the corresponding library name as found in the filenames, and the number of samples pooled in the library.
  ### CellBender_Snakemake 
  The code for removing ambient RNA (background) using CellBender. 
  
## scATAC-seq alignment
  ### scATAC_pipeline_lite
  The code for aligning the scATAC portion of single cell multiomic sequencing to a reference assembly. The codebase was forked from https://github.com/austintwang/scATAC_pipeline_lite (commit: 76a4bc545318b5685eb9843e94000dde0d3394e7) with modifications for file and result organization. 

## scATAC-seq quality control
  ### multiomic_ATAC_QC_Snakemake
  The code for performing quality control on the scATAC portion of the multiome libraries. As defined in the config file, the default filters are minTSS >= 6, and minFrags >= 1000, unless overriden in the sample sheet. Only cells deemed valid by [STARsolo](#rna_alignment_snakemake) were considered for this analysis. Both "addDoubletScores" from [ArchR](https://www.archrproject.com/bookdown/inferring-scatac-seq-doublets-with-archr.html) and [Amulet](https://github.com/UcarLab/AMULET) were used to identify doublets. 

## scRNA-seq quality control and normalization
  ### RNA_QC_Snakemake
  The code for generating the necessary QC metrics of scRNA-seq using the outputs from [ambient RNA removal](#cellbender_snakemake). 
  ### Normalization_Snakemake
  The code for filtering cells based on the thresholds determined by manual examination of the QC metrics generated by [RNA_QC_Snakemake](#rna_qc_snakemake). The filtered scRNA-seq count matrices are then normalized and [Scublet](https://github.com/swolock/scrublet) is used to identify doublets. 
  
## Doublet removal
  ### Doublet_removal_Snakemake
  The code for rerunning Scrublet if specified in the SampleSheet (when "Doublet_threshold" is not "Default"), removing the union of doublets identified by Scrublet, [Archr](#multiomic_atac_qc_Snakemake), [Amulet](#multiomic_atac_qc_Snakemake), and [Souporcell](#rna_alignment_snakemake). 
## Clustering and cell type annotation
  ### evaluation_Snakemake
  The code for merging filtered Seurat objects, normalizing the resulting count matrices, and performing cell clustering at various resolutions. The pipeline generates heatmaps reflecting marker gene expression in the clusters and tentatively labels clusters by mapping the cells to an annotated earlier version of the heart map. 
  ### clustering_Snakemake
  The code for performing clustering on each major cell type group. The code performs count spliting (splitting the count matrix into training and test sets), normalization, and clustering of the training set at various resolutions. Subsequently, the test set is labeled using the clustering results from the training set. Finally, the expression of marker genes in these clusters is highlighted in both the test dataset and the original dataset. 

  ### annotation_Snakemake
  The code for faciliating the annotation of fine-grained clusters. It generates heatmaps visualizaing the top 10 differentially expressed genes in each cluster at the resolution specified in the sample sheet. It also generates dot plots and feature plots of specified marker genes. It calculates the TPM of genes in each cluster. 

## Defining gene programs with consensus non-negative matrix factorization
  ### cNMF_run_and_analysis/Snakefile_cNMF
  The code for running [cNMF](https://github.com/dylkot/cNMF) on count matrices specified in the SampleSheet, using the array of K values specified in the config file. 
  ### cNMF_run_and_analysis/Snakefile_cNMF_analysis
  The code for performing further analysis on cNMF results, including: rerunning the consensus step with an updated local density filter; organizing output files; summarizing the distribution of gene programs across various features (e.g., cell annotation, sex); identifying Gene Ontology (GO) terms associated with gene programs; quantifying the variance in gene expression explained by the gene programs; and generating additional plots for K selection.

## Identifying transcription factors regulating gene programs 
  ### ChromVar_Snakemake
  The code for using [chromVar](https://github.com/GreenleafLab/chromVAR) to calculate per-cell motif accessibility and output the Z-scores of [bias corrected deviations](https://greenleaflab.github.io/chromVAR/articles/Articles/Deviations.html) of motif accessibility in each cell.
  ### get_motif_enrichment
  This directory contains code for calculating the correlation bewteen per-cell motif accessibility and the expression of gene programs. 

## Assessing cell type heritability enrichment with stratified LD score regression
  ### LDSC_Snakemake 
  This directory contains code for LD score regression.
## Preparing variants for overlapping with scE2G-predicted enhancers
  ### LD_clumping_expansion
  The code for calculating beta and odds ratios for variants in each dataset. It also performs LD clumping and expansion to identify the lead variants and variants in LD with the lead variants at each locus.
## Linking GWAS variants to enhancers, target genes, cell types, and gene programs
  ### V2G_v2
  The code for linking GWAS variants to cell-type specific enhancers, target genes and gene programs. 

## Testing disease genes for preferential expression in particular cell types
  ### CHD_enrichment
  
## Assessing enrichment of genes associated with valve traits in VIC program 27 
  ### VIC_27_enrichment





