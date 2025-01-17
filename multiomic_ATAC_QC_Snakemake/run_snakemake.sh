#!/bin/bash
project="/oak/stanford/groups/engreitz/Users/rosaxma/230828_heart_map/atac_alignment/pool_3"
snakemake="/oak/stanford/groups/engreitz/Users/rosaxma/github_repos_2/multiomic_ATAC_QC_Snakemake"
cd ${snakemake}

LOG=${project}/logs_both_filter_multiomics_custom_filters
mkdir -p ${LOG}

snakemake \
	--profile ${snakemake}/sherlock \
	--configfile ${project}/config/config_both_filter_multiomics_custom_filters.yaml \
	--snakefile ${snakemake}/workflow/Snakefile \
	--cluster "sbatch -n 1 -p engreitz,owners,normal -J {rule} -o ${LOG}/{rule}_{wildcards}.qout -e ${LOG}/{rule}_{wildcards}.e  --cpus-per-task {threads} --mem {resources.mem_gb}GB --time {resources.runtime_hr}:00:00" --rerun-incomplete
