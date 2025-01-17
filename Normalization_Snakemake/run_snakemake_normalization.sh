#!/bin/bash
project="alignment/pool_3"
config=config/config_normalization.yaml
LOG=${project}/logs
mkdir -p ${LOG}
snakemake \
	--rerun-incomplete \
	--use-conda \
	--conda-frontend mamba \
	--profile sherlock \
	--configfile ${config} \
	--snakefile workflow/Snakefile \
	--cluster "sbatch -n 1 -p engreitz,owners,normal -J {rule} -o ${LOG}/{rule}_{wildcards}.qout -e ${LOG}/{rule}_{wildcards}.e --cpus-per-task {threads} --mem {resources.mem_gb}GB --time {resources.runtime_hr}:00:00"

