"""Panel b tables: the 23 genes in human valve interstitial cells by post-conception week.

Source: Bayraktar, Cranley et al. 2024 (bioRxiv 2024.04.27.591127), developing human heart atlas,
log-normalised RNA object from the Heart Cell Atlas foetal portal (CC BY 4.0, 5.8 GB):
  https://cellgeni.cog.sanger.ac.uk/heartcellatlas/foetal_webportal/RNA_global_lognormalised.h5ad
297,473 cells and nuclei x 36,601 genes. `.X` is the authors' log1p(CP10K) (expm1(.X) sums to
10,000 per cell) and is used as-is; cell types are the authors' `fine_grain` annotation; stage is
`week` (4W-20W, post-conception weeks).

Port of analysis/vic_chd_developmental_expression/scripts/taskC_cranley_valve_expression.py
(EngreitzLab/FetalHeartMap-paper-internal @ ad1d03d): same cell types, weeks, gene order and
arithmetic, without the plots, so it rewrites the two shipped tables byte for byte. The figure
uses only the ValveInterstitialCells rows; the other six types are kept for reference.

  python build_cranley_vic_tables.py --h5ad /path/to/RNA_global_lognormalised.h5ad

Outputs (overwritten): data/taskC_cranley_expr_by_celltype_week_tidy.tsv (cell type x week) and
data/taskC_cranley_early_pcw4to6_summary.tsv (pcw4-6 pooled; sets the figure's gene row order).
"""
import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from gene_set import GENES as GENES_TSV, EXPRESSION_ONLY

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"

# Same order as upstream (expression-only group first), so the tables match byte for byte.
EXPR_ONLY = ["ABL1", "ARID1B", "BMPR2", "FGFR2", "GLI3", "GPC3", "INVS", "PRKACB", "PRKD1", "RBFOX2"]
OTHER = ["CREBBP", "FBN1", "GATA4", "GATA6", "KANSL1", "LTBP2", "MEIS2", "NF1",
         "NOTCH2", "SETD5", "SMAD6", "TBX20", "TGFBR1"]
GENES = EXPR_ONLY + OTHER
assert set(GENES) == set(GENES_TSV) and set(EXPR_ONLY) == set(EXPRESSION_ONLY), \
    "gene lists drifted from vic_nominated_chd_genes_fig2.tsv"
HEADLINE = EXPR_ONLY + ["GATA4", "GATA6", "MEIS2"]

WEEK_ORDER = ["4W", "5W", "6W", "7W", "8W", "9W", "10W", "12W", "13W", "14W", "15W", "20W"]
EARLY_WEEKS = ["4W", "5W", "6W"]
CT_ORDER = [
    "ValveInterstitialCells",           # VICs -- the only type the figure uses
    "EndocardialCells",
    "MyocardialInterstitialFibroblasts",
    "SubEpicardialFibroblasts",
    "Myofibroblasts",
    "EpicardiumDerivedCells",
    "MesothelialEpicardialCells",
]

ap = argparse.ArgumentParser()
ap.add_argument("--h5ad", required=True, help="RNA_global_lognormalised.h5ad")
args = ap.parse_args()

a = ad.read_h5ad(args.h5ad, backed="r")
assert a.shape == (297473, 36601), f"unexpected object shape {a.shape}"
assert set(a.obs["week"].astype(str)) == set(WEEK_ORDER), "week values differ from WEEK_ORDER"
vn = list(map(str, a.var_names))
sub = a[:, [vn.index(g) for g in GENES]].to_memory()
X = sub.X.toarray() if sp.issparse(sub.X) else np.asarray(sub.X)
df = pd.DataFrame(X, index=sub.obs_names, columns=GENES).join(a.obs[["fine_grain", "week"]])
df = df[df["fine_grain"].isin(CT_ORDER)].copy()
df["week"] = df["week"].astype(str)


def summarise(s, **keys):
    return [dict(gene=g, **keys, n_cells=len(s),
                 mean_lognorm=float(np.mean(s[g].values)),
                 pct_expressing=float(np.mean(s[g].values > 0) * 100),   # log1p > 0 iff count > 0
                 expression_only_group=(g in EXPR_ONLY), headline_group=(g in HEADLINE))
            for g in GENES]


rows = []
for ct in CT_ORDER:
    for wk in WEEK_ORDER:
        s = df[(df["fine_grain"] == ct) & (df["week"] == wk)]
        if len(s):
            rows += summarise(s, cell_type=ct, week=wk)
pd.DataFrame(rows).to_csv(DATA / "taskC_cranley_expr_by_celltype_week_tidy.tsv", sep="\t", index=False)

early = df[df["week"].isin(EARLY_WEEKS)]
erows = []
for ct in CT_ORDER:
    s = early[early["fine_grain"] == ct]
    if len(s):
        erows += summarise(s, cell_type=ct, weeks="pcw4-6")
pd.DataFrame(erows).to_csv(DATA / "taskC_cranley_early_pcw4to6_summary.tsv", sep="\t", index=False)

vic = df[df["fine_grain"] == "ValveInterstitialCells"]
print(f"ValveInterstitialCells: {len(vic)} cells; per week "
      + ", ".join(f"{w} {int((vic.week == w).sum())}" for w in WEEK_ORDER))
print("wrote data/taskC_cranley_expr_by_celltype_week_tidy.tsv, "
      "data/taskC_cranley_early_pcw4to6_summary.tsv")
