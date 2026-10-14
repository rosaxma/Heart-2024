# Fetal heart map — iPSC-VIC Perturb-seq analysis

Snakemake workflows for the iPSC-derived valve interstitial cell (VIC) high-MOI
CRISPRi Perturb-seq screen targeting congenital heart disease (CHD) genes.

**Contact:** Maggie Chen (msychen@stanford.edu), Engreitz Lab

**Stages so far:**

| Workflow | What it does |
|---|---|
| `00_counts_matrix_processing.smk` | kb counts → QC → clustering + CellRanger export |
| `01_sceptre_trans.smk` | SCEPTRE trans: calibration check, power check, discovery |
| `02_per_guide_knockdown.smk` | per-gRNA knockdown of the CHD target panel |

**Standalone analyses** (outside the Snakemake workflows; each folder has its
own README and `environment.yml`):

| Folder | What it does |
|---|---|
| `bulk_meis2_kd/` | bulk RNA-seq validation of the MEIS2 knockdown in iPSC-VICs (pydeseq2) |
| `vic_chd_developmental_expression/` | expression of the 23 VIC-nominated CHD genes across mouse and human valve development (Supplementary Fig. 24) |

---

## Stage 00 — counts-matrix processing

Picks up the raw per-sublibrary `adata.h5ad` files written by
[`tkzeng/perturb_pipeline`](https://github.com/tkzeng/perturb_pipeline) and
produces two things:

1. a clustered, annotated combined AnnData — the input to cNMF;
2. per-sublibrary CellRanger feature-barcode matrices — the input to SCEPTRE.

Scripts in `scripts/` are numbered in the order they run.

```
kb counts (GEX + guide, per sublibrary)
  │
  ├─ 01_umi_threshold_filter.py            rule: umi_filter
  │     Per-sublibrary UMI floor from the upstream barcode-rank inflection
  │     point, shared ceiling of 75,000. Intersects GEX ∩ guide barcodes.
  │     → umi_filtered/{sublibrary}_umi_filtered_gex_and_guide.h5mu
  │
  ├─ 02_qc_metrics_filter.py               rule: qc_filter
  │     MT ≤ 15%, ribo ≤ 4%, ±4 MAD outliers on log1p total counts and
  │     log1p n_genes, then per-gene filters (≥10 cells, ≥100 counts).
  │     → qc_filtered/{sublibrary}_qc_filtered_gex_and_guide.h5mu
  │
  ├─ 03a_concatenate_mudata.py             rule: concatenate_sublibraries
  │  │    Outer join across sublibraries, barcodes suffixed with the label.
  │  │    → combined/combined_all_sublibraries.h5mu
  │  │
  │  └─ 04_cluster_and_annotate.py         rule: cluster_and_annotate
  │        CP10k + log1p, cell-cycle scoring, 3,000 HVGs (seurat_v3,
  │        batch_key=sublib), 50 PCs / 30 used, k=15 neighbours, UMAP,
  │        Leiden at 0.1/0.25/0.5/1.0, marker scoring and dotplots.
  │        → clustering/combined_clustered.h5ad              [feeds cNMF]
  │
  └─ 03b_build_master_gene_list.py         rule: build_master_gene_list
     │    Union of gene IDs across all sublibraries.
     │    → cellranger_format_featurecorrected/master_gene_ids.txt
     │
     ├─ 03c_convert_to_cellranger.py       rules: convert_to_cellranger,
     │  │                                  convert_to_cellranger_featurecorrected
     │  │  Per-sublibrary barcodes/features/matrix triplets. Counts are
     │  │  rounded to int32 here. Genes stacked above guides
     │  │  ("Gene Expression" / "CRISPR Guide Capture").
     │  │  Without --master-gene-list → cellranger_format/{sublibrary}/…
     │  │  With    --master-gene-list → cellranger_format_featurecorrected/…
     │  │
     └─ 03d_validate_cellranger_features.py  rule: validate_cellranger_features
           Confirms every sublibrary emitted an identical features.tsv.gz.
           → cellranger_format_featurecorrected/…            [feeds SCEPTRE]
```

`03a`/`04` and the `03b`→`03c`→`03d` export branch are independent branches off
`02`, so they run in parallel. `03c` without `--master-gene-list` does not
depend on `03b` and can run at any point after `02`.

These map onto the numbered directories of the original analysis tree as:
`01` → `1.umifiltering`, `02` → `2.qcfiltering`, `03a` →
`3a.concat-old-forclustering`, `03c` → `3b.cellranger_format`, `03b`+`03c`+`03d`
→ `3b-i.cellranger_format_featurecorrected`, `04` → `4.clustering`.

### Why the feature-corrected export exists

`02_qc_metrics_filter.py` applies its per-gene filters one sublibrary at a time,
so the sublibraries end up with different gene sets — for the 2026-04-08 run,
25,703 / 25,451 / 25,266 / 24,522 / 25,234 / 22,960 / 17,354 genes for
sublibraries 1–7, union 26,070. SCEPTRE's `import_data_from_cellranger` reads a list of matrix
directories and indexes features **positionally**, so mismatched
`features.tsv.gz` files corrupt every result rather than raising. The
`validate_cellranger_features` rule fails loudly if they ever diverge.

---

## Stage 01 — SCEPTRE trans differential expression

Tests every perturbation against every protein-coding gene genome-wide, with a
calibration check on negative-control pairs and a power check on the TSS
positive controls. All published results come from the trans analysis.

```
CellRanger matrices (stage 00, feature-corrected)
  │
  ├─ 05_extract_sceptre_covariates.py   rule: extract_sceptre_covariates
  │     S_score, G2M_score, pct_ribo, pct_mito read out of the clustered
  │     AnnData, in CellRanger import order.
  │     → sceptre_ondisc/sceptre_extra_covariates.csv
  │
  ├─ 06c_validate_cell_ordering.py      rule: validate_cell_ordering
  │     Confirms row i of that CSV is sceptre cell i. The covariate join is
  │     positional with no barcode key, so this check is load-bearing.
  │
  ├─ 06a_import_into_sceptre.R          rule: import_into_sceptre
  │     import_data_from_cellranger(moi="high"), ondisc-backed so the full
  │     240k x 26k matrix never has to be held in memory.
  │     → sceptre_ondisc/{sceptre_object.rds, gene.odm, grna.odm}
  │
  ├─ 06b_set_analysis_parameters.R      rule: set_analysis_parameters
  │     Attaches the covariates, builds the trans pair set, restricts responses
  │     to protein-coding genes, sets the association formula explicitly.
  │     → sceptre_ondisc/sceptre_object_trans_v2.rds
  │
  ├─ nextflow sceptre-pipeline          rule: run_sceptre_pipeline
  │     Wraps `nextflow run timothy-barry/sceptre-pipeline`, which does the
  │     gRNA assignment and all the association testing. ~12 h wall clock.
  │     → sceptre/trans_v2/outputs/results_run_*.rds
  │
  └─ 06d_export_sceptre_results.R       rule: export_sceptre_results
        Adds gene symbols and pct_knockdown = 100 * (1 - 2^log2FC).
        → sceptre/trans_v2/sceptre_trans_v2_*.tsv
```

Covariate sets and pair types are declared under `sceptre.variants`;
`sceptre.run_variants` picks which ones `rule all` builds: `trans_v2`
(cell-state covariates, no Leiden cluster), the trans discovery analysis behind
every published SCEPTRE result, and `cis`, which is built only because stage 02
reads its sceptre object and gRNA assignment.

---

## Stage 02 — per-gRNA knockdown

SCEPTRE reports one fold change per target, pooling that target's gRNAs. This
stage recomputes knockdown one gRNA at a time over the CHD gene panel, by two
independent estimators — they share the counts and the gRNA assignment but
nothing else, so their agreement is evidence the knockdown is real rather than
a modelling artefact.

```
cis sceptre object + gRNA assignment matrix (stage 01)
  │
  ├─ 07a_per_guide_fold_changes_poisson.R      x2, independent
  │     Poisson GLM log2FC against sceptre's own null model
  │     (sceptre book §10.4). One script, two guide sets:
  │       --guide-set targeting      → every pair in the cis object: the
  │                                    TSS on-target pairs plus their cis
  │                                    neighbours (~18k rows); only the
  │                                    on-target rows are plotted
  │       --guide-set non-targeting  → the empirical null
  │
  ├─ 07b_target_gene_expression_per_guide.R
  │     Model-free CP10K pseudobulk vs a per-gene non-targeting pool
  │     (NT cells carrying none of that gene's TSS gRNAs).
  │
  ├─ 07c_plot_tss_knockdown.py     Mann-Whitney on log2FC → the per-gRNA dots
  ├─ 07d_plot_nt_vs_tss.py         Mann-Whitney on CP10K  → stars + ordering,
  │                                and the correlation between the two estimators
  │
  └─ 07e_plot_pct_knockdown_bar.py
        The combined figure. Bars are the mean of each gene's own dots.
        19 columns, the main-figure panel E genes in panel E's order
        (plots.panel_e_groups; 8 panel genes excluded).
        → per_guide_knockdown/fig_pct_knockdown_bar.svg
```

`07a` (targeting), `07a` (non-targeting) and `07b` are mutually independent
despite the numbering, and all three read only stage-01 outputs.

**Percent knockdown is `100 * (1 - 2^log2FC)`, a concave transform.** Do not
average percentages across gRNAs and compare that to the percentage of a mean
fold change — average log2 fold changes and convert once. `07e` asserts the
resulting Jensen inequality, along with five other guards (excluded genes must
exist, `plots.exclude_genes` and `plots.panel_e_groups` must partition the
27-gene panel, every plotted gene must have dots, each bar must equal the mean
of its dots to 1e-9, and the gene-label geometry must match the figure spec).

---

## Software

Snakemake 7.32.4 (the `--cluster` interface, not the 8+ `--executor slurm`
plugin) with Python 3.10.19 and:

| Package | Version | Package | Version |
|---|---:|---|---:|
| anndata | 0.11.4 | numpy | 2.2.6 |
| h5py | 3.15.1 | openpyxl | 3.1.5 |
| matplotlib | 3.10.8 | pandas | 2.3.3 |
| mudata | 0.3.3 | pynndescent | 0.6.0 |
| muon | 0.1.7 | python-igraph | 1.0.0 |
| natsort | 8.4.0 | PyYAML | 6.0.3 |
| numba | 0.64.0 | scanpy | 1.11.5 |
| scikit-learn | 1.7.2 | scikit-misc | 0.5.2 |
| scipy | 1.15.3 | seaborn | 0.13.2 |
| umap-learn | 0.5.11 | | |

Two environment files, so the orchestrator and the analysis stack solve
independently:

```bash
conda activate mamba-env          # see the lab manual, 03-software.md
mamba env create -f environment.yml        # -> fetalheartmap
mamba env create -f envs/snakemake.yml     # -> fetalheartmap-snakemake
```

`leidenalg` is deliberately absent: `04_cluster_and_annotate.py` uses
`sc.tl.leiden(flavor="igraph")`, which goes through `python-igraph`.

The versions above are transcribed from the environment the published run used
(`$OAK/Users/msychen/conda/envs/perturbprocessing`, which `config.conda.env`
points at). CI rebuilds `environment.yml` from scratch on every push, so a
solve that has drifted shows up as a red build rather than as a surprise months
later.

`envs/snakemake.yml` pins Python 3.11. An older Snakemake environment exists at
`~/.conda/envs/snakemake` on Python 3.8, which still parses these workflows but
emits `importlib.metadata has no attribute packages_distributions` on startup
and picks up `~/.local/lib/python3.8/site-packages`. Build the environment from
`envs/snakemake.yml` instead, or prefix commands with `PYTHONNOUSERSITE=1`.

---

## Layout

```
config/       one YAML per analysis; all paths and parameters live here
references/   small inputs that ship with the code (sample table, gene lists,
              gRNA targets, Ensembl->symbol map)
scripts/      step scripts, numbered in run order, plus shared helpers
              (pipeline_utils.py, snakemake_helpers.py, r_utils.R, figstyle.py)
workflows/    numbered Snakemake stages (*.smk)
envs/         the orchestrator environment
tests/        contract tests that run without the sequencing data
submit.sh     Slurm submission wrapper
```

Script numbers run across the whole analysis, not per stage: `01`–`04` are
stage 00, `05`–`06` stage 01, `07` stage 02.

R scripts take their arguments through `scripts/r_utils.R`, a ~100-line base-R
parser. That is deliberate: the shared R library these scripts run against has
neither `optparse` nor `getopt`, and its `argparse` is a reticulate wrapper
around Python's, which would make every R script depend on a working Python
inside R. The SCEPTRE pipeline's own scripts use plain `commandArgs()` for the
same reason.

Paths under `config/` and `references/` resolve relative to this code
directory. Runtime paths (`scratch_base`, `results_base`) resolve relative to
wherever Snakemake is launched; set them to absolute paths if you prefer.

---

## Required inputs

### 1. kb-python count matrices (from GEO)

The pipeline starts from the 14 unfiltered kb-python count objects deposited
in GEO (accession [GSE349107](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE349107)), one GEX and one guide object per sublibrary:

```
sublib{1..7}_gex_adata.h5ad     # kb nac, X = mature + nascent + ambiguous
sublib{1..7}_guide_adata.h5ad   # kb kite, 3,026 gRNAs
```

Download all 14 into one directory and point `input_paths.kb_output_base` at
it. The shipped config uses `geo_counts`, relative to the directory you launch
Snakemake from:

```bash
mkdir geo_counts            # then place the 14 files here
```

For each sublibrary the workflow reads

```
{kb_output_base}/[{counts_dir}/][{counts_subdir}/]{counts_file}
```

where `counts_file` comes from the sample table. The two directory columns are
optional; they are blank for the flat GEO layout. The GEX objects hold
fractional counts as written by kb. They stay float through QC and are rounded
(`np.round`, half to even, so 0.5 becomes 0) at the CellRanger export (step
3b-i), as in the published run. Padding each sublibrary to the master gene list
upcasts the matrix back to float64, so the feature-corrected `matrix.mtx.gz` is
written as `real`; its values are still whole numbers.

On Sherlock, the same files live in the nested kb output tree on Oak
(`…/20260407_mergedfastq_barcoderanksinflection_300kround3feb_ipscvicd9_newguidelist_analysis/merged:*/kb_*_main_raw/counts_unfiltered_modified/sublib*_adata.h5ad`).
Either point `kb_output_base` at that tree and fill in `counts_dir`
(`merged:iPSC-VICs_sublibrary_1_merged`, …) and `counts_subdir`
(`kb_all_main_raw/counts_unfiltered_modified`, …) in the sample table, or
keep the shipped table and symlink the files into one flat directory:

```bash
mkdir geo_counts
ln -s /oak/stanford/groups/engreitz/Users/msychen/fetalheartmap/20260408-FEBVICHIGHIMOI-300kchd-d9-novaseq/countsmatrices/20260407_mergedfastq_barcoderanksinflection_300kround3feb_ipscvicd9_newguidelist_analysis/merged:*/kb_*_main_raw/counts_unfiltered_modified/sublib*_adata.h5ad geo_counts/
```

Run from that layout, step 1 reproduces the published step-1 output exactly
(checked for sublib7: the same 9,011 cells, the same features in order, and
zero count difference).

### 2. Sample table — `references/sample_info.<analysis>.tsv`

One row per library. Required columns:

| Column | Meaning |
|---|---|
| `sublibrary` | Sublibrary label; the unit of work, and the barcode suffix |
| `sample_id` | Unique library identifier |
| `sample_type` | `gex` or `guide` |
| `counts_file` | Count AnnData file name, e.g. `sublib1_gex_adata.h5ad` |
| `counts_dir` | Optional sublibrary directory under `kb_output_base` |
| `counts_subdir` | Optional further subdirectory holding `counts_file` |
| `min_umi_threshold` | UMI floor (gex rows); the barcode-rank inflection point |
| `max_umi_threshold` | UMI ceiling (gex rows); blank means no ceiling |
| `exclude` | `True` to drop a library without editing the workflow |

Every sublibrary needs exactly one `gex` and one `guide` row — `tests/` enforces
this.

### 3. Ensembl → symbol map

`input_paths.ensembl_to_symbol` points at `references/ens2symbol_v43.tsv`
(GENCODE v43, columns `ensembl_id` / `gene_symbol`), shipped here. It is
identical to the pickled `ensembl_to_symbol_dict_v43.npy` the published run
read (all 67,524 entries; the same 37 MT- and 1,518 RPS/RPL genes among the
63,086 in the data), so QC and clustering are unchanged. A `.npy` dict is
still accepted. `symbol_to_ensembl` is not used by any step. Setting
`ensembl_to_symbol` to `""` skips symbol mapping, but then the MT/ribo QC
annotation and the marker panels do not resolve, so the QC step keeps a
different set of cells.

### 4. Cell-cycle gene lists

Shipped in `references/` (Tirosh et al. 2016, 42 S and 54 G2M genes). Two G2M
symbols, `FAM64A` and `HN1` (now `PIMREG` and `JPT1`), are not in GENCODE v43,
so 42 S and 52 G2M genes are scored.

### 5. gRNA target table — `references/guide_targets.*.tsv`

`grna_id` / `grna_target` (plus `chr`/`start`/`end`), 3,026 gRNAs: 2,826
targeting across 179 targets, and 200 non-targeting. `grna_target` is an
Ensembl gene ID for a TSS target, an element name for an enhancer, or the
literal `non-targeting`. Shipped in `references/`.

### 6. Clustered AnnData, for the SCEPTRE covariates

The covariates `S_score`, `G2M_score`, `pct_ribo` and `pct_mito` come from a
clustered AnnData. With `sceptre.covariates.clustered_h5ad: ""` (the default),
stage 01 reads stage 00's own output, `{results_base}/clustering/combined_clustered.h5ad`,
and declares it as an input, so no extra file is needed. The published runs
read the cell-cycle-regressed clustering instead. Those four covariates are
byte-identical between it and the no-regression output (verified 2026-09-08),
so either reproduces them. Only `leiden_res0_25` differs between the two (6
clusters vs 5), and the default variant `trans_v2` does not use it. To
reproduce `trans_v3` exactly, set `clustered_h5ad` to the published file and
`barcode_h5ad_infix` to `_300k`; that file names cells `BARCODE_300k_{sublibrary}`.

### 7. GENCODE annotation and the R toolchain

Outside Sherlock, four config values point at Stanford-local resources and
need changing:

| Setting | Change to |
|---|---|
| `sceptre.gencode_gtf` | GENCODE v43 `gencode.v43.chr_patch_hapl_scaff.annotation.gtf.gz`, from https://www.gencodegenes.org/human/release_43.html |
| `sceptre.r_env` | `modules: ""` and `r_libs_user: ""`, with `sceptre` 0.10.2 and `ondisc` 1.2.0 installed in your R library |
| `conda.env` | your environment built from `environment.yml` (or `""` to use the one on PATH) |
| `sceptre.pipeline.work_dir_base`, `resources.*.slurm_partition` | your scratch directory and Slurm partitions |
| `sceptre.pipeline.nextflow_config` | `config/nextflow.sherlock.config`: Slurm queue, per-process memory, and the module loads and R library path the Nextflow subjobs run with |


`sceptre.gencode_gtf` (GENCODE v43, 53 MB) restricts trans discovery pairs to
protein-coding responses. `sceptre.r_env` names the R modules and the shared
library providing `sceptre` 0.10.2, `ondisc` 1.2.0, `data.table` and
`rtracklayer`. `nextflow` and `java` must be on PATH for the pipeline rule.

Note that library has neither `optparse` nor `getopt`, which is why the R
scripts parse arguments through `scripts/r_utils.R` instead.

### 8. CHD gene panel — `references/chd_genes.suppt8.tsv`

27 genes (23 highly expressed in VICs, 4 not) with the flags stage 02 uses to
pick the main-figure subset. Shipped
in `references/`.

---

## Running

```bash
conda activate fetalheartmap-snakemake
export CONFIG=config/config.20260408_ipscvic_300k.yaml

./submit.sh --dry-run        # check the DAG
./submit.sh --jobs 20        # one Slurm job per rule
```

Per-rule memory, walltime, threads and partition come from the `resources:`
block of the config; `submit.sh` maps them onto `sbatch`. Cluster stdout/stderr
land in `{results_base}/logs/cluster/`, per-rule script logs in
`{results_base}/logs/{rule}/`.

`submit.sh` defaults to stage 00. Select a later stage with `WORKFLOW`, and run
them in order — each depends on the one before:

```bash
export CONFIG=config/config.20260408_ipscvic_300k.yaml

WORKFLOW=workflows/00_counts_matrix_processing.smk ./submit.sh
WORKFLOW=workflows/01_sceptre_trans.smk            ./submit.sh
WORKFLOW=workflows/02_per_guide_knockdown.smk      ./submit.sh
```

Stage 01's Nextflow rule asks for a 3-day walltime and only 10 GB, because it
mostly waits on the subjobs it submits itself. Do not assume `-resume` will
save you a re-run: changing the analysis parameters changes that rule's input
hash, and everything downstream of gRNA assignment re-runs with it.

To run everything inside a single allocation instead:

```bash
conda activate fetalheartmap-snakemake
snakemake -s workflows/00_counts_matrix_processing.smk \
  --configfile config/config.20260408_ipscvic_300k.yaml --cores 8
```

The default target builds the clustering output and the feature-corrected
export. The pre-correction export is a separate target:

```bash
./submit.sh cellranger_export_plain
```

### Tests

```bash
conda activate fetalheartmap
pytest -q
```

These check the config/sample-table/workflow contract — required keys, that
`dotplot_resolution` is one of `leiden_resolutions`, that every rule declares a
log and resources, that every script is argparse-driven and referenced by a
workflow, and that no absolute path has crept back into the code. They need no
sequencing data, so CI runs them on every push.

---

## Adapting this to a new experiment

1. Copy `config/config.20260408_ipscvic_300k.yaml`, change `analysis_name`,
   `results_base` and `input_paths.kb_output_base`.
2. Copy `references/sample_info.20260408_ipscvic_300k.tsv`, fill in your
   sublibraries and their UMI bounds, and point the new config at it.
3. Set `qc_filtering.reference` to `mouse` if applicable, and replace
   `clustering.marker_panels` with panels for your cell types.
4. `pytest -q`, then `./submit.sh --dry-run`.

For another human experiment no code changes should be needed. Mouse needs
some: `04_cluster_and_annotate.py` hard-codes the human `MT-` / `RPS` / `RPL`
prefixes, and the cell-cycle lists and `ens2symbol_v43.tsv` are human-only.

---

## Relationship to the other lab repos

| Repo | Overlap with this stage |
|---|---|
| [`tkzeng/perturb_pipeline`](https://github.com/tkzeng/perturb_pipeline) | Produces this stage's inputs. FASTQ → kallisto/bustools → cell calling → per-sublibrary `adata.h5ad`. Not vendored here. |
| [`EngreitzLab/telohaec-genomewide-perturb-seq-analysis`](https://github.com/EngreitzLab/telohaec-genomewide-perturb-seq-analysis) | Same conventions (numbered `.smk` stages, `config/`, `references/`, `scripts/`, `pipeline_utils` + `snakemake_helpers`), and this repo follows them. The *analyses* differ — see below. |


## Changes made when importing these scripts into the repo

The scripts came from `20260408-postprocessing/` as loose `sbatch` submitters
with hardcoded absolute paths. **No analysis logic, parameter or threshold was
changed.** What changed:

- **Snakemake replaces the `sbatch` loops.** The `for SUBLIB in 1 2 3 4 5 6 7`
  submitters and the `--dependency=afterok` chain that built the 3b-i master
  gene list are now rules. The chain became two rules
  (`build_master_gene_list`, then a fan-out), which removes the dependency on
  sublibrary 1 running first.
- **All paths moved into `config/` and `references/`.** Nothing in `scripts/` or
  `workflows/` contains a user- or project-specific absolute path (the only
  absolute path is the system `/etc/profile.d/modules.sh` fallback in
  `snakemake_helpers.py`); `tests/` enforces that.
- **3b and 3b-i are one script.** `03c_convert_to_cellranger.py
  --master-gene-list` is the 3b-i behaviour; without the flag it is 3b. The
  duplicated file also drifted — 3b lacked 3b-i's feature-uniqueness assert,
  which is now applied in both modes.
- **The `reference_features.txt` cross-check became its own rule.** In the
  original it was written by whichever sublibrary ran first and compared against
  by the rest, which races once the jobs run in parallel.
  `validate_cellranger_features` compares all sublibraries after the fact and
  reports per-sublibrary feature counts.
- **`1.umi_threshold_int_layercheck.py` folded into
  `01_umi_threshold_filter.py`** as `--check-layer-sums` and `--round-counts`. It
  was a diagnostic variant (asserts `nascent + ambiguous + mature == X`, rounds
  to int32); nothing downstream consumed its output. Neither flag is set by the
  workflow, so the default path is unchanged.
- **Random seeds are now explicit.** `clustering.random_seed` is passed to PCA,
  the neighbour graph, UMAP, Leiden and gene scoring. The original relied on
  scanpy's defaults.
- **`use_highly_variable=True` → `mask_var="highly_variable"`** in `sc.tl.pca`,
  which is the non-deprecated spelling of the same thing.
- **`TEST_MODE` toggles removed** from the two CellRanger scripts.

