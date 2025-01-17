#!/bin/bash
config=config/config.yaml

LOG=logs
mkdir -p ${LOG}
snakemake \
  --keep-going \
  --use-conda \
  --rerun-incomplete \
  --conda-frontend mamba \
  --configfile ${config} \
  --snakefile  ${snakemake}/workflow/Snakefile \
  --cluster "sbatch -n 1 -p engreitz,owners,normal -J {rule} -o ${LOG}/{rule}_{wildcards}.qout -e ${LOG}/{rule}_{wildcards}.e --cpus-per-task {threads} --mem {resources.mem_gb}GB --time {resources.runtime_hr}:00:00"
