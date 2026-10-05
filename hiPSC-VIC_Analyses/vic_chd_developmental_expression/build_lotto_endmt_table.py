"""Panel a table: mouse EndMT cell state x timepoint expression of the 23 genes.

Source: Lotto et al. 2023 (Nat Commun 14:5567) mouse valve-development atlas,
Broad Single Cell Portal SCP1897.

Default: reads the committed extract in data/lotto_scp1897/ for the 10,354
cells of the authors' EndMT subset (portal file cluster/endo_clusterfile.txt):
  lotto_endmt_cell_annotation.tsv.gz  clusters, timepoints, pseudotime (verbatim)
  lotto_endmt_counts_23genes.tsv.gz   raw UMI counts of the 23 genes
  lotto_endmt_total_umi.tsv.gz        per-cell total UMI over ALL genes

  python build_lotto_endmt_table.py

To rebuild the extract from the portal download (raw files are not committed;
see README.md), pass the SCP1897 download folder:

  python build_lotto_endmt_table.py --raw /path/to/SCP1897

Output: tables/lotto_endmt_celltype_timepoint.tsv -- per gene x cell state x
timepoint: n_cells, mean log1p(CP10K), % cells with count > 0, and the cell
state's mean pseudotime.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from gene_set import GENES_MOUSE

HERE = Path(__file__).resolve().parent
LOTTO = HERE / "data" / "lotto_scp1897"
ANNOT = LOTTO / "lotto_endmt_cell_annotation.tsv.gz"
COUNTS = LOTTO / "lotto_endmt_counts_23genes.tsv.gz"
TOTAL = LOTTO / "lotto_endmt_total_umi.tsv.gz"
OUT = HERE / "tables" / "lotto_endmt_celltype_timepoint.tsv"
TP_ORDER = ["E7.75", "E8.25", "E9.25", "E10.5", "E12.5"]


def extract(scp):
    """Annotation from cluster/endo_clusterfile.txt (row 2 is the portal's TYPE
    line), then one streaming pass over expression/total_counts2.csv (genes x
    cells): per-cell total over all genes, and the 23 genes' counts, kept for
    the EndMT cells only."""
    scp = Path(scp)
    endo = pd.read_csv(scp / "cluster" / "endo_clusterfile.txt", sep="\t",
                       skiprows=[1], index_col="NAME")
    endo["pseudotime"] = pd.to_numeric(endo["pseudotime"])
    endo[["clusters", "timepoints", "pseudotime"]].rename_axis("cell").to_csv(ANNOT, sep="\t")
    raw = scp / "expression" / "total_counts2.csv"
    want = set(GENES_MOUSE)
    with open(raw) as fh:
        barcodes = fh.readline().rstrip("\n").split(",")[1:]
        total = np.zeros(len(barcodes), dtype=np.int64)
        keep = {}
        for line in fh:
            gene, rest = line.rstrip("\n").split(",", 1)
            gene = gene.strip('"')
            v = np.array(rest.split(","), dtype=np.int64)
            assert v.size == len(barcodes), f"{gene}: {v.size} values"
            total += v
            if gene in want:
                assert gene not in keep, f"duplicated gene row {gene}"
                keep[gene] = v
    assert set(keep) == want, f"missing genes: {want - set(keep)}"
    counts = pd.DataFrame({g: keep[g] for g in GENES_MOUSE}, index=barcodes)
    total = pd.Series(total, index=barcodes, name="total_umi")
    assert endo.index.isin(counts.index).all(), "EndMT cells missing from the matrix"
    counts.loc[endo.index].rename_axis("cell").to_csv(COUNTS, sep="\t")
    total.loc[endo.index].rename_axis("cell").to_csv(TOTAL, sep="\t")
    print(f"extracted {len(endo)} EndMT cells x {len(GENES_MOUSE)} genes from {raw}")


ap = argparse.ArgumentParser()
ap.add_argument("--raw", help="SCP1897 download folder (holds cluster/ and expression/)")
args = ap.parse_args()
if args.raw:
    extract(args.raw)

endo = pd.read_csv(ANNOT, sep="\t", index_col="cell")
counts = pd.read_csv(COUNTS, sep="\t", index_col=0).loc[endo.index, GENES_MOUSE]
total = pd.read_csv(TOTAL, sep="\t", index_col=0)["total_umi"].loc[endo.index]
assert (total > 0).all(), "cell with zero total UMI"
lognorm = np.log1p(counts.div(total, axis=0) * 1e4)          # log1p(CP10K)

grp = [endo["clusters"], endo["timepoints"]]
mean = lognorm.groupby(grp).mean().stack().rename("mean_lognorm")
pct = ((counts > 0).groupby(grp).mean() * 100).stack().rename("pct_expressing")
n = lognorm.groupby(grp).size().rename("n_cells")
tab = pd.concat([mean, pct], axis=1).reset_index()
tab.columns = ["cell_type", "timepoint", "gene", "mean_lognorm", "pct_expressing"]
tab = tab.merge(n.rename_axis(["cell_type", "timepoint"]).reset_index())
pt = endo.groupby("clusters")["pseudotime"].mean().rename("cell_type_mean_pseudotime")
tab = tab.merge(pt, left_on="cell_type", right_index=True)
tab["_t"] = tab.timepoint.map(TP_ORDER.index)
tab = (tab.sort_values(["cell_type_mean_pseudotime", "_t", "gene"])
       [["gene", "cell_type", "timepoint", "n_cells", "mean_lognorm", "pct_expressing",
         "cell_type_mean_pseudotime"]])
OUT.parent.mkdir(exist_ok=True)
tab.to_csv(OUT, sep="\t", index=False)
print(f"wrote {OUT.relative_to(HERE)} ({len(endo)} EndMT cells, {len(tab)} rows)")
