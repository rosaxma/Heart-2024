# VIC-nominated CHD genes across valve-lineage development — figure

Everything needed to regenerate one figure: expression of the 23 VIC-nominated
CHD genes (Fig. 2) in

- **a** mouse endocardium → AV mesenchyme (Lotto et al. 2023, SCP1897),
- **b** human valve interstitial cells by post-conception week (Bayraktar,
  Cranley et al. 2024),
- **c** VICs of our fetal heart atlas (VIC_1–5 pooled) by post-conception
  week, 6–11 pcw and 12–22 pcw pooled.

This folder is self-contained: the figure needs no imports or reads outside it.

## Run

```bash
python build_lotto_endmt_table.py   # -> tables/lotto_endmt_celltype_timepoint.tsv  (seconds)
python make_figure.py               # -> figures/combined_valve_lineage_expression_revised.{svg,png,pdf}
```

The panel b and c tables are shipped. To rebuild them from source:

```bash
python build_cranley_vic_tables.py --h5ad RNA_global_lognormalised.h5ad   # panel b, ~10 s
python build_fhm_vic_by_pcw_table.py                                     # panel c, Sherlock only
```

- **Panel b:** the `.h5ad` is the public Heart Cell Atlas foetal object, 5.8 GB,
  from https://cellgeni.cog.sanger.ac.uk/heartcellatlas/foetal_webportal/RNA_global_lognormalised.h5ad.
  On 2026-09-30 the script rewrote both shipped tables byte for byte.
- **Panel c:** the script needs the atlas's per-cell counts and donor ages
  (Engreitz lab access). It reproduces the shipped table byte for byte.

Tested with Python 3.10.19, numpy 2.2.6, pandas 2.3.3 and matplotlib 3.10.8
(matplotlib ≥ 3.10 is required). Use the **SVG**: text is live and set in
Arial. The PNG and PDF fall back to DejaVu Sans where Arial is not installed.
`make_figure.py` prints the panel-a columns and cell counts, the groups dropped
for having fewer than 50 cells, the panel-c columns with cell and donor counts,
and the per-gene E12.5 AV mesenchyme minus endocardium differences.

## Files

| file | used for | source |
|---|---|---|
| `make_figure.py` | the figure | — |
| `build_lotto_endmt_table.py` | panel a table | — |
| `build_cranley_vic_tables.py` | panel b tables, from the public `.h5ad` | port of upstream `taskC_cranley_valve_expression.py` |
| `build_fhm_vic_by_pcw_table.py` | panel c table (Sherlock only) | — |
| `FIGURE_LEGEND.md`, `METHODS.md` | legend and Methods text for this figure | — |
| `gene_set.py`, `vic_nominated_chd_genes_fig2.tsv` | the 23 genes, human→mouse symbols | `analysis/vic_chd_developmental_expression/` |
| `data/lotto_scp1897/lotto_endmt_cell_annotation.tsv.gz` | panel a cells, cell states, stages, pseudotime | `clusters`, `timepoints`, `pseudotime` columns of SCP1897 `cluster/endo_clusterfile.txt` (10,354 cells), values verbatim |
| `data/lotto_scp1897/lotto_endmt_counts_23genes.tsv.gz` | panel a | raw UMI counts of the 23 genes for the 10,354 EndMT cells, from SCP1897 `expression/total_counts2.csv` |
| `data/lotto_scp1897/lotto_endmt_total_umi.tsv.gz` | panel a | per-cell total UMI over **all** 21,862 genes of that matrix, same cells |
| `data/taskC_cranley_expr_by_celltype_week_tidy.tsv` | panel b | written by `build_cranley_vic_tables.py` |
| `data/taskC_cranley_early_pcw4to6_summary.tsv` | gene row order (human pcw4–6 VIC mean) | same |
| `data/fhm_vic_by_pcw.tsv` | panel c | written by `build_fhm_vic_by_pcw_table.py` from the atlas's 13,636 VIC cells; per-age-bin summaries only (gene, mean log1p(CP10K), % expressing, cells, donors, median UMI, weeks pooled) |
| `tables/lotto_endmt_celltype_timepoint.tsv` | panel a input | written by `build_lotto_endmt_table.py` |

The two `data/taskC_*` tables were first produced in `EngreitzLab/FetalHeartMap-paper-internal`
at commit `ad1d03d`; `build_cranley_vic_tables.py` reproduces them exactly. The
23-gene list (`vic_nominated_chd_genes_fig2.tsv`) is the one input not derived
here: it is the Fig. 2 gene set, copied verbatim.

Only processed extracts and summaries are committed; no raw data. For panel c,
per-cell values are not shipped because they identify individual donors
(donor ID and gestational week); the figure uses only the age-bin summaries.
md5 of the uncompressed inputs (gzip checksums change on every rebuild because
of the header timestamp, so compare those with `zcat`):
```
f644039d5f1d457c9c3cde0eed45083f  data/taskC_cranley_expr_by_celltype_week_tidy.tsv
e5bb1c433b4fc6be16d762de206f1a09  data/taskC_cranley_early_pcw4to6_summary.tsv
a941a72ce83eaad9743e742aeceeac21  data/fhm_vic_by_pcw.tsv
9242fdd4a614bf2b7dbcfd8ecfa33800  vic_nominated_chd_genes_fig2.tsv
c98810e5268f738dbf3eabd98a91998a  SCP1897 expression/total_counts2.csv (not committed)
```

## Rebuilding the mouse extract from the raw matrix

The raw matrix is 2.4 GB, which is over GitHub's file-size limit, so only the
23-gene extract is committed. To rebuild it from source:

1. Sign in at https://singlecell.broadinstitute.org/single_cell/study/SCP1897.
   Open **Download → Bulk download** and copy the `curl` command.
2. Run only the `curl ... -o cfg.txt` part. In `cfg.txt`, delete the
   `url=`/`output=` pairs you don't need; `total_logcounts2.csv` (22 GB) is
   never needed. Then run `curl -K cfg.txt`.
3. Run `python build_lotto_endmt_table.py --raw SCP1897`. It needs
   `cluster/endo_clusterfile.txt` and `expression/total_counts2.csv` in that
   folder, takes about 90 s, and overwrites the three `.tsv.gz` extracts. On
   2026-09-28 the rebuilt extracts, table and figure were identical to the
   committed ones.

## Methods in brief (panel a)

- **Expression:** log1p(CP10K) = log1p(count / cell total × 10⁴), computed from
  raw counts. The portal's `total_logcounts2.csv` is Seurat `ScaleData` output
  (per-gene z-scores), so it is not used.
- **Columns:** the authors' EndMT cell states (`clusters`) on the endocardium →
  AV mesenchyme path, in mean-pseudotime order (asserted in `make_figure.py`),
  at each stage (`timepoints`) with ≥ 50 cells. Dot colour is the mean
  log1p(CP10K); dot area is the % of cells with count > 0.
- **Not shown:**
  - early endothelium, which comes before endocardium;
  - VECs I/II (valvular endothelium, a side branch that stays endothelial);
  - five groups with 1–43 cells.
- **Sources of the stages:** E7.75–E9.25 are de Soysa et al. 2019 whole hearts
  (10x 3′ v2). E10.5–E12.5 are Lotto et al. dissected AV canals (mostly 10x
  3′ v3), integrated by Lotto et al. under one annotation.
- **Depth check:** thinning every cell to 10,000 UMI changes the colour value
  by ≤ 0.14 in every column and leaves the endocardium → mesenchyme increases in
  place.

## Methods in brief (panel c)

- **Cells:** all 13,636 cells of the five VIC subpopulations, from 40 donors.
  Raw counts come from the atlas's per-subtype matrices. Each cell's donor and
  post-conception week come from the atlas annotation, joined on sample ID +
  ATAC barcode + subtype. The build script asserts that the join is unique and
  complete, that each cell's UMI total matches the atlas, and that pooled counts
  reproduce the atlas pseudobulk TPM (r > 0.999).
- **Columns:** 6W–11W individually (555–4,693 cells, 2–8 donors each) and
  12–22W pooled (848 cells, 12 donors; 1–2 donors per week there).
- **Values:** as in panel a. Dot colour is the mean log1p(CP10K); dot area is
  the % of cells with count > 0.

Panels a and b share one absolute colour scale (PowerNorm γ = 0.55, 0–2.6).
Panel c uses the same ramp with its own range (0–3.4): the atlas is about ten
times shallower per cell (median 2,314 UMI), which raises log1p(CP10K) for
detected genes, so panel c colours are not comparable with a and b.

## Citations

- Lotto J, et al. *Cell diversity and plasticity during atrioventricular heart
  valve EMTs.* Nat Commun 2023;14:5567. doi:10.1038/s41467-023-41279-6
- de Soysa TY, et al. *Single-cell analysis of cardiogenesis reveals basis for
  organ-level developmental defects.* Nature 2019;572:120–124.
  doi:10.1038/s41586-019-1414-x
- Bayraktar S, Cranley J, et al. *High-resolution atlas of the developing human
  heart and the great vessels.* bioRxiv 2024.04.27.591127 (preprint).
