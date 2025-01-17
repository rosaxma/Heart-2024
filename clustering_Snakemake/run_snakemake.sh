#!/bin/bash
project="clustering"
config="config/config.yaml"
LOG=${project}/logs
mkdir -p ${LOG}

snakemake \
  -n \
  --keep-going \
  --rerun-incomplete \
  --conda-frontend mamba \
  --profile ${snakemake}/sherlock \
  --configfile ${config} \
  --snakefile  ${snakemake}/workflow/Snakefile \
  --cluster "sbatch -n 1 -p {resources.partition} -J {rule} -o ${LOG}/{rule}_{wildcards}.qout -e ${LOG}/{rule}_{wildcards}.e --cpus-per-task {threads} --mem {resources.mem_gb}GB --time {resources.runtime_hr}:00:00"  
