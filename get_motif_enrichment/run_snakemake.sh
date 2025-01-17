#!/bin/bash
snakemake="/oak/stanford/groups/engreitz/Users/rosaxma/github_repos_2/get_motif_enrichment"
LOG=logs
mkdir -p ${LOG}
snakemake \
	--reason \
	--keep-going \
	--rerun-incomplete \
	--conda-frontend mamba \
	--use-conda \
	--profile sherlock \
	--configfile config/config.yaml \
	--snakefile ${snakemake}/workflow/Snakefile_2 \
	--cluster "sbatch -n 1 -p engreitz,normal,owners -J {rule} -o ${LOG}/{rule}_{wildcards}.qout -e ${LOG}/{rule}_{wildcards}.e --cpus-per-task {threads} --mem {resources.mem_gb}GB --time {resources.runtime_hr}:00:00" 


