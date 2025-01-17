#!/bin/bash
project=alignment/pool_3
config=config/config.yaml

LOG=${project}/logs
mkdir -p ${LOG}
snakemake \
  -n \
  --rerun-incomplete \
  --conda-frontend mamba \
  --profile ${snakemake}/sherlock \
  --configfile ${config} \
  --snakefile ${snakemake}/workflow/Snakefile \
  --cluster "sbatch -n 1 -p gpu,owners -J {rule} -o ${LOG}/{rule}_{wildcards}.qout -e ${LOG}/{rule}_{wildcards}.e --ntasks 1 -C GPU_MEM:{resources.gpu_mem_gb} -G 1 --time {resources.runtime_hr}:00:00 --mem {resources.cpu_mem_gb}" 

