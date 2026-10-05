# fetalheartmap — bulk MEIS2 knockdown validation

Bulk RNA-seq validation of the MEIS2 trans signature from the iPSC-VIC
high-MOI CRISPRi Perturb-seq screen, **in iPSC-derived VICs only — not primary
VICs**.

**Contact:** Charu Balamurugan (charuba@stanford.edu), Engreitz Lab

`scripts/for_FHM_MEIS2KD.py` highlights a 15-gene fetal-heart-map (FHM)
panel on the combined MEIS2 KD vs NT volcano. It is the script that produced
the figures and tables cited in the paper.

---

## What it does

Takes the **unstimulated** samples only from both experiments (6 NT vs 6 MEIS2
KD across 12 libraries), rebuilds the count matrix from the raw per-experiment
matrices, and fits with pydeseq2 at `ALPHA = 0.05` after a `total count >= 10`
filter. Doubling the sample count for the one contrast of interest is the
point: a gene genuinely down in the KD should come through tighter than either
experiment gives alone.

The script used to have a second, read-only half that highlighted the panel on
each experiment's already-written `MEIS2_KD_effect_unstim` table, producing
`wnt9b`- and `wnt5b`-tagged outputs. That was removed — nothing consumed those
outputs, and dropping it means the raw count matrices and metadata are now the
script's only inputs, so everything it needs is shipped here. The
per-experiment volcanoes are still produced by each experiment's own
`run_deseq2_v2_all_biorep.py` if they are wanted.

The combined fit is run **twice** and both are saved, so they can be compared
directly:

| Design | Tag | |
|---|---|---|
| `~ genotype` | `combined` | all 12 libraries as one cohort |
| `~ experiment + genotype` | `combined_batch` | **the design used for the paper** |

### Why the experiment term matters

The layout is balanced — 3 NT + 3 KD per experiment — so `~ genotype` cannot
*bias* the KD estimate. Batch only confounds genotype when the two are
correlated, and here they are evenly crossed. But it forces one NT level to
describe both experiments, so the between-experiment baseline gap (PBX1
`baseMean` 2742 vs 1059; FOXP2 67 vs 27) is charged to the residual, inflating
dispersion and widening every `lfcSE`.

That is severe, not cosmetic: **PBX1 goes from `lfcSE` 0.041 in the wnt9b
experiment alone to 0.346 pooled, and from `padj` 2e-25 to 0.70.** Adding the
per-experiment offset absorbs the gap into a coefficient and estimates the KD
effect from *within*-experiment differences (paired rather than unpaired).
`genotype[T.MEIS2]` means exactly the same thing under both designs and still
uses all 12 samples.

On-target anchor under `combined_batch`: MEIS2 itself `log2FC = -0.742`,
`lfcSE = 0.049`, Wald *z* = -15.04, `padj = 5.5e-48`. 1,948 genes reach
`padj < 0.05`.

---

## Inputs

### Shipped in this repo

```
data/MEIS2_wnt9b/metadata.csv                      216 B
data/MEIS2_wnt9b/9YFNNW-expression-matrix.csv.gz   3.4 MB  (16 MB raw)
data/MEIS2_wnt5b/metadata.csv                      213 B
data/MEIS2_wnt5b/2SZ9JF-expression-matrix.csv.gz   3.4 MB  (16 MB raw)
```

**This is everything the script needs** — it reads nothing else. Both matrices
are 78,986 genes x 27 columns (78,987 lines with the header).

`metadata.csv` is one row per library — `sample,genotype,treatment`. The
`2a`/`2b` suffixes in wnt9b (`2-NT2a`/`3-NT2b`, `8-M2a`/`9-M2b`) are the
technical-replicate pairs behind the pseudoreplication caveat below, and
`treatment` is what the unstimulated subset filters on (6 of 12 per experiment).

The script reads the gzipped matrices directly; no need to decompress them.

Upstream of these matrices (FASTQs, alignments) is not redistributed here.
The 24 raw FASTQs, one per library (`{accession}_{n}_{sample}.fastq.gz`,
e.g. `9YFNNW_7_7-M1.fastq.gz`), are deposited in GEO under accession
[GSE349277](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE349277). They include the 12
Wnt-treated libraries, which this script does not use.

### Layout the script expects

`BASE` is this directory (`bulk_meis2_kd/`), found relative to the script, so
it runs from any checkout without edits. (It was originally an absolute Oak
path with one project directory per experiment; outputs are unchanged by the
move, checked byte-for-byte on all six tables.) Per experiment it reads:

```
data/MEIS2_<experiment>/<accession>-expression-matrix.csv.gz
                                          genes x (gene_id, gene_name, gene_biotype,
                                          then a <sample>_cpm / <sample>_count pair
                                          per sample)
data/MEIS2_<experiment>/metadata.csv      one row per library, indexed by sample id,
                                          with `genotype` (NT | MEIS2) and
                                          `treatment` (none | Wnt) columns
```

Accessions are wired in `EXPERIMENTS`: `wnt9b` → `9YFNNW`, `wnt5b` → `2SZ9JF`.

Only `_count` columns are read; the `_cpm` columns are ignored. featureCounts
assigns multimappers fractionally, so counts are rounded to int.

---

## Outputs

All into `bulk_meis2_kd/combined_MEIS2KD_unstim/` (git-ignored), two tags — `combined` and
`combined_batch`, one per design:

```
deseq2_MEIS2KD_unstim_<tag>.csv                  genome-wide, unshrunk
deseq2_MEIS2KD_unstim_<tag>_shrunk.csv           genome-wide, shrunk LFC
highlighted_FHMpanel_MEIS2KD_unstim_<tag>.tsv    the 15-gene panel
volcano_MEIS2KD_unstim_FHMpanel_<tag>.png        panel on the volcano
```

A design-comparison table (LFC / `lfcSE` / `padj` side by side across the two
combined designs) is printed to stdout but **not saved** — it is the clearest
justification for choosing `combined_batch`, so capture the log if you need it.

### The panel TSV has two LFC columns and they are not interchangeable

```
gene_id  gene_name  baseMean  log2FoldChange  LFC_shrunk  lfcSE  pvalue  padj  gene_biotype
                              └─ unshrunk ─┘  └─ shrunk ┘
```

`log2FoldChange` is the raw Wald estimate and is the one `pvalue`/`padj`
correspond to. `LFC_shrunk` is the apeGLM-shrunk estimate, joined in from the
shrunk fit. **Quote `log2FoldChange` whenever an effect size sits next to its
p-value, and `LFC_shrunk` when reporting effect size alone.** The top hits
barely differ (JAG1 -0.612 vs -0.600) but PBX3 does: -0.382 vs -0.272, sitting
at `padj` 0.0500.

Related trap: in the genome-wide `_shrunk.csv` the `pvalue`/`padj` columns are
carried over from the Wald test, **not** recomputed from the shrunk LFC. There
is no "shrunk test".

Under `combined_batch` all 15 panel genes are present and testable (no
independent-filtering dropouts) and 11 clear `padj < 0.05`.

---

## Running

Needs a compute node — compute nodes have no bare `python3` on `PATH`, and the
loaded Lmod `py-numpy`/`python` modules leak in via `PYTHONPATH` and break
numpy's C extensions, hence `env -u PYTHONPATH`:

From `bulk_meis2_kd/`:

```bash
srun -n1 -p dev --time=00:20:00 --mem=16G --cpus-per-task=4 \
  env -u PYTHONPATH MPLBACKEND=Agg \
  <conda-env>/bin/python scripts/for_FHM_MEIS2KD.py
```

`N_CPUS = 4` at the top of the script is a **reproducibility** setting, not
just a speed knob: pydeseq2's parallel dispersion fitting is not bit-identical
across worker counts, so changing it moves genes across the `padj = ALPHA`
boundary. Keep it fixed between runs you compare.

---

## Caveats this fit carries

**Pseudoreplication.** All 12 libraries are treated as independent biological
replicates, by decision. wnt9b's `2a`/`2b` and `M2a`/`M2b` pairs are two
libraries from one biological replicate, so the true biological n is below 6 —
dispersion is underestimated and the p-values are anti-conservative.

**Effect-size heterogeneity is not modelled.** KD effects differ 2-3x between
the two experiments (WNT5A -0.50 vs -1.49; FOXP2 -0.74 vs -2.79), so the pooled
LFC is an average. A gene whose effect is real in only one experiment (LGR4,
SALL4, PBX3, RSPO3) comes out **attenuated rather than confirmed**. Testing
that heterogeneity directly needs an `experiment:genotype` term, which neither
design includes.

**Panel matching is by symbol only.** Genes are matched on `gene_name`, not
symbol-or-Ensembl-ID. Fine for a curated symbol panel, and anything not found
is reported explicitly rather than silently dropped.

---

## Environment

`environment.yml` pins the env this was actually run in (pydeseq2 0.5.4,
pandas 3.0.5, numpy 2.5.2). Recreate with:

```bash
conda env create -f environment.yml -n deseq2
```
