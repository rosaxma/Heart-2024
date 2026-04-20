#!/bin/bash
project=""
snakemake="slingshot_snakemake"
config=config/config.yaml

LOG=${project}/logs
mkdir -p ${LOG}
cd ${snakemake}
snakemake \
  --keep-going \
  --use-conda \
  --rerun-incomplete \
  --conda-frontend mamba \
  --profile ${snakemake}/sherlock \
  --configfile ${config} \
  --snakefile  ${snakemake}/workflow/Snakefile_phase1 \
  --cluster "sbatch -n 1 -p engreitz,owners,normal -J {rule} -o ${LOG}/{rule}_{wildcards}.qout -e ${LOG}/{rule}_{wildcards}.e --cpus-per-task {threads} --mem {resources.mem_gb}GB --time {resources.runtime_hr}:00:00" 


snakemake \
  --keep-going \
  --use-conda \
  --rerun-incomplete \
  --conda-frontend mamba \
  --profile ${snakemake}/sherlock \
  --configfile ${config} \
  --snakefile  ${snakemake}/workflow/Snakefile_phase2 \
  --cluster "sbatch -n 1 -p engreitz,owners,normal -J {rule} -o ${LOG}/{rule}_{wildcards}.qout -e ${LOG}/{rule}_{wildcards}.e --cpus-per-task {threads} --mem {resources.mem_gb}GB --time {resources.runtime_hr}:00:00"
