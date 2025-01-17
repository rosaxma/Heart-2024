#!/bin/bash
project="/oak/stanford/groups/engreitz/Users/rosaxma/230828_heart_map/multiomic_evaluation/all_samples"
snakemake="/oak/stanford/groups/engreitz/Users/rosaxma/github_repos_2/evaluation_subset_Snakemake"
config=${project}/config/config_evaluation.yaml

LOG=${project}/logs
mkdir -p ${LOG}
cd ${snakemake}

snakemake \
  -n \
  --forceall \
  --rerun-incomplete \
  --conda-frontend mamba \
  --profile ${snakemake}/sherlock \
  --configfile ${config} \
  --snakefile  ${snakemake}/workflow/Snakefile \
  --cluster "sbatch -n 1 -p {resources.partition} -J {rule} -o ${LOG}/{rule}_{wildcards}.qout -e ${LOG}/{rule}_{wildcards}.e --cpus-per-task {threads} --mem {resources.mem_gb}GB --time {resources.runtime_hr}:00:00"  
