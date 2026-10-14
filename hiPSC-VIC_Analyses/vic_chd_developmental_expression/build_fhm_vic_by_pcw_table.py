"""Panel c table: our fetal heart atlas VICs (VIC_1-5 pooled) by post-conception week.

Same metrics as panels a and b: per cell log1p(count / cell total x 1e4), with the cell total summed
over every gene in the matrix; then per age bin, the mean over cells and the % of cells with
count > 0. Weeks >= POOL_FROM (12-22W, 1-2 donors each) are pooled into one bin.

Sources (read-only):
  counts  Rosa's per-subtype raw count matrices (10x mtx, full transcriptome, the matrices behind
          GEX_TPM / TPM_long.tsv), barcodes "<sample_id>_<ATAC barcode>":
            rosaxma/230828_heart_map/multiomic_pseudobulk/240516_subset_rna_mtx/<subtype>/matrix/
  ages    cell -> donor (`source`) -> `Rounded.PCW` from the atlas object obs:
            msychen/fetalheartmap/20260312-fhm-obj-seurat-andata/fhm_slim_2k.h5ad

Needs read access to both sources (Engreitz lab, Sherlock). The figure does not: its output,
data/fhm_vic_by_pcw.tsv, is shipped in this folder. Per-cell data stay in memory and are not
written here, because they identify individual donors (donor ID + gestational week).

Checks, all asserted: every VIC barcode joins to exactly one atlas cell of the same subtype; each
cell's total UMI equals the atlas nCount_RNA; pooled pseudobulk CPM tracks Rosa's TPM (r > 0.999).

Run (on Sherlock):
  python build_fhm_vic_by_pcw_table.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from gene_set import GENES  # noqa: E402  (human symbols, 23)

ROSA = Path("/oak/stanford/groups/engreitz/Users/rosaxma/230828_heart_map/"
            "multiomic_pseudobulk/240516_subset_rna_mtx")
OBJ = Path("/oak/stanford/groups/engreitz/Users/msychen/fetalheartmap/"
           "20260312-fhm-obj-seurat-andata/fhm_slim_2k.h5ad")
OUT = HERE / "data" / "fhm_vic_by_pcw.tsv"
SUBTYPES = {"VIC_free_segments_1": "VIC_1", "VIC_free_segments_2": "VIC_2",
            "VIC_free_segments_3": "VIC_3", "VIC_NC_derived_mescenchyme": "VIC_4",
            "VIC_valve_roots": "VIC_5"}
POOL_FROM = 12


def extract():
    import anndata as ad
    import scipy.io

    obs = ad.read_h5ad(OBJ, backed="r").obs
    v = obs[obs.celltype.astype(str).isin(SUBTYPES)].copy()
    v["key"] = (v.sample_id.astype(str) + "_" + v.ATAC_barcode.astype(str) + "|"
                + v.celltype.astype(str))
    assert not v.key.duplicated().any(), "barcode key not unique in the atlas"
    v = v.set_index("key")

    blocks = []
    for ct in SUBTYPES:
        d = ROSA / ct / "matrix"
        bc = pd.read_csv(d / "barcodes.tsv.gz", header=None)[0]
        feat = pd.read_csv(d / "features.tsv.gz", header=None)[0]
        X = scipy.io.mmread(d / "matrix.mtx.gz").tocsc()            # genes x cells
        assert X.shape == (len(feat), len(bc))
        assert np.allclose(X.data, np.round(X.data)), f"{ct}: non-integer counts"
        meta = v.reindex(bc + "|" + ct)
        assert meta.celltype.notna().all(), f"{ct}: barcodes missing from the atlas"
        total = np.asarray(X.sum(axis=0)).ravel()
        assert np.array_equal(total, meta.nCount_RNA.to_numpy()), f"{ct}: totals != nCount_RNA"
        gi = pd.Index(feat).get_indexer(GENES)
        assert (gi >= 0).all(), f"missing genes: {[g for g, i in zip(GENES, gi) if i < 0]}"
        b = pd.DataFrame(X[gi, :].T.toarray().astype(int), columns=GENES)
        b.insert(0, "cell", bc.values)
        b.insert(1, "subtype", SUBTYPES[ct])
        b.insert(2, "donor", meta.source.astype(str).values)
        b.insert(3, "pcw", meta["Rounded.PCW"].astype(int).values)
        b.insert(4, "total_umi", total.astype(int))
        blocks.append(b)
        print(f"{ct}: {len(bc)} cells, all joined")
    cells = pd.concat(blocks, ignore_index=True)

    # pseudobulk CPM per subtype vs Rosa's TPM for the 23 genes
    tpm = pd.read_csv(ROSA / "TPM_2" / "TPM_long.tsv", sep="\t")
    tpm = tpm[tpm.Celltype.isin(SUBTYPES) & tpm.Gene.isin(GENES)].pivot(
        index="Gene", columns="Celltype", values="TPM").rename(columns=SUBTYPES)
    cpm = (cells.groupby("subtype")[GENES].sum()
           .div(cells.groupby("subtype").total_umi.sum(), axis=0).T * 1e6)
    for s in cpm.columns:
        r = np.corrcoef(np.log1p(cpm[s]), np.log1p(tpm.loc[cpm.index, s]))[0, 1]
        assert r > 0.999, f"{s}: pseudobulk vs TPM r = {r:.4f}"
    print("pseudobulk CPM vs Rosa TPM: r > 0.999 for all 5 subtypes")

    return cells


if not (ROSA.exists() and OBJ.exists()):
    sys.exit("Atlas sources not readable here (Sherlock / Engreitz lab only). "
             "data/fhm_vic_by_pcw.tsv is already shipped; the figure does not need this step.")
cells = extract()
assert len(cells) == 13636 and (cells.total_umi > 0).all()

counts = cells[GENES]
ln = np.log1p(counts.div(cells.total_umi, axis=0) * 1e4)
late = cells.pcw >= POOL_FROM
late_label = f"{POOL_FROM}–{cells.pcw.max()}W"
cells["bin"] = np.where(late, late_label, cells.pcw.astype(str) + "W")
order = [f"{p}W" for p in sorted(cells.pcw[~late].unique())] + [late_label]

tab = pd.concat({"mean_lognorm": ln.groupby(cells.bin).mean().stack(),
                 "pct_expressing": (counts > 0).groupby(cells.bin).mean().stack() * 100}, axis=1)
tab.index.names = ["bin", "gene"]
meta = cells.groupby("bin").agg(n_cells=("cell", "size"), n_donors=("donor", "nunique"),
                                median_umi=("total_umi", "median"),
                                weeks=("pcw", lambda x: ",".join(map(str, sorted(x.unique())))))
tab = tab.reset_index().merge(meta.reset_index())
tab["bin_order"] = tab["bin"].map(order.index)
tab = tab.sort_values(["bin_order", "gene"]).drop(columns="bin_order")
tab.to_csv(OUT, sep="\t", index=False)
print(meta.reindex(order).to_string())
print(f"wrote {OUT.relative_to(HERE)}")
