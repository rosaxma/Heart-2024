#!/bin/bash
project="/oak/stanford/groups/engreitz/Users/rosaxma/230828_heart_map/multiomic_cnmf/complete_heartmap"

LOG=${SCRATCH}/cnmf_logs
mkdir -p ${LOG}
snakemake \
	--rerun-incomplete \
	--keep-going \
	--conda-frontend mamba \
	--profile sherlock \
	--configfile ${project}/config/config_analysis.yaml \
	--snakefile workflow/Snakefile_cNMF_analysis \
	--cluster "sbatch -n 1 -p engreitz,normal,owners -J {rule} -o ${LOG}/{rule}_{wildcards}.qout -e ${LOG}/{rule}_{wildcards}.e --cpus-per-task {threads} --mem {resources.mem_gb}GB --time {resources.runtime_hr}:00:00" 
