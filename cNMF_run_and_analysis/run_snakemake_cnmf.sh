#!/bin/bash
LOG=${SCRATCH}/cnmf_logs
mkdir -p ${LOG}

snakemake \
	-n \
	--rerun-incomplete \
	--keep-going \
	--conda-frontend mamba \
	--profile sherlock \
	--configfile config/config.yaml \
	--snakefile workflow/Snakefile_cNMF \
	--cluster "sbatch -n 1 -p engreitz,normal,owners -J {rule} -o ${LOG}/{rule}_{wildcards}.qout -e ${LOG}/{rule}_{wildcards}.e --cpus-per-task {threads} --mem {resources.mem_gb}GB --time {resources.runtime_hr}:00:00"

