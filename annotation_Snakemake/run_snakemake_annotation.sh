#!/bin/bash
project="clustering"
config="config/config_annotation.yaml"

LOG=${project}/logs
mkdir -p ${LOG}

snakemake \
  -n \
  --rerun-incomplete \
  --conda-frontend mamba \
  --profile sherlock \
  --configfile ${config} \
  --snakefile  workflow/Snakefile \
  --cluster "sbatch -n 1 -p {resources.partition} -J {rule} -o ${LOG}/{rule}_{wildcards}.qout -e ${LOG}/{rule}_{wildcards}.e --cpus-per-task {threads} --mem {resources.mem_gb}GB --time {resources.runtime_hr}:00:00"  
