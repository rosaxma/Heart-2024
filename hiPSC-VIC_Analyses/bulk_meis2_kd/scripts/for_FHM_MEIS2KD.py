"""FHM gene panel on the combined MEIS2 KD vs NT volcano.

Takes the UNSTIMULATED samples only from both experiments -- 6 NT vs 6 MEIS2 across 12 libraries --
and fits from the raw count matrices, writing only to combined_MEIS2KD_unstim/. This doubles the
sample count for the one contrast of interest, so a gene genuinely down in the KD should come
through tighter than either experiment gives alone.

The per-experiment half was removed: it read each experiment's already-written
`deseq2_v2_MEIS2_KD_effect_unstim.csv` and highlighted the panel on it, producing the `wnt9b` and
`wnt5b` tagged outputs. Nothing consumed those, and dropping it means this script no longer depends
on anything in either deseq2_results_all_biorep/ -- the raw count matrices and metadata are now its
only inputs. The per-experiment volcanoes still exist in each experiment's own
run_deseq2_v2_all_biorep.py if they are wanted.

Design notes:

  All 12 libraries are treated as independent biological replicates, by decision. wnt9b's 2a/2b and
  M2a/M2b pairs are two libraries from one biological replicate, so the true biological n is lower
  than 6; dispersion is therefore underestimated and the p-values are anti-conservative. Same
  caveat as the header of run_deseq2_v2_all_biorep.py.

  The combined fit is run TWICE, once per design in DESIGNS, and both are saved so they can be
  compared directly:

    ~ genotype              all 12 libraries as one cohort, ignoring which experiment they came
                            from. The layout is balanced (3 NT + 3 MEIS2 per experiment) so this
                            cannot BIAS the KD estimate -- batch only confounds genotype when the
                            two are correlated, and here they are crossed evenly. But it forces one
                            NT level to describe both experiments, so the between-experiment
                            baseline gap (PBX1 baseMean 2742 vs 1059; FOXP2 67 vs 27) is charged to
                            the residual, inflating dispersion and widening every lfcSE. Severely
                            conservative here: PBX1 goes from lfcSE 0.041 in wnt9b alone to 0.346
                            pooled, and from padj 2e-25 to 0.70.

    ~ experiment + genotype adds a per-experiment offset, so the baseline gap is absorbed by a
                            coefficient instead of the error term and the KD effect is estimated
                            from the WITHIN-experiment differences (paired rather than unpaired).
                            genotype[T.MEIS2] means exactly the same thing and still uses all 12
                            samples. Recommended reading of the two.

  Neither design fixes the effect-size heterogeneity below; both assume one common KD effect.

  The KD effect sizes themselves differ 2-3x between experiments (WNT5A -0.50 vs -1.49; FOXP2 -0.74
  vs -2.79), so the pooled LFC is an average. A gene whose effect is real in only one experiment
  (LGR4, SALL4, PBX3, RSPO3) comes out attenuated rather than confirmed. Testing that heterogeneity
  directly would need an `experiment:genotype` term, which is not in scope here.

Inputs are read from bulk_meis2_kd/data/ and outputs written to bulk_meis2_kd/combined_MEIS2KD_unstim/,
both located relative to this script, so it runs from any checkout with no path edits. Run it on a
dev node with the deseq2 conda env (compute nodes have no bare python3 on PATH):

  srun -n1 -p dev --time=00:20:00 --mem=16G --cpus-per-task=4 \
    env -u PYTHONPATH MPLBACKEND=Agg <deseq2-env>/bin/python scripts/for_FHM_MEIS2KD.py
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # batch script, no display on a compute node

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from adjustText import adjust_text
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

# bulk_meis2_kd/, the directory holding scripts/ and data/. Originally an absolute Oak path
# (/oak/stanford/groups/engreitz/Users/charuba/20260811_iPSC_VIC_stimscreens) with one project
# directory per experiment; now the shipped data/ layout, so no edit is needed to run it elsewhere.
BASE = Path(__file__).resolve().parents[1]
DATA_DIR = BASE / "data"
OUTPUT_DIR = BASE / "combined_MEIS2KD_unstim"
OUTPUT_DIR.mkdir(exist_ok=True)

ALPHA = 0.05
N_CPUS = 4

# (design, filename tag, figure label). Both are fitted and saved; see the docstring.
DESIGNS = [
    ("~ genotype", "combined", "combined 6v6, no batch term"),
    ("~ experiment + genotype", "combined_batch", "combined 6v6, + experiment term"),
]

NS_COLOR = "#B0B0B0"     # recessive gray for the rest of the transcriptome
PANEL_COLOR = "#9B2D3E"  # deep maroon, same as the volcanoes in run_deseq2_v2_all_biorep.py

FHM_PANEL = [
    "PBX1", "PBX3", "FOXP2", "FOXP1", "KAT6B", "SALL4", "JAG1", "GPC3",
    "WNT5A", "ROR1", "ROR2", "PRICKLE1", "KIF26B", "RSPO3", "LGR4",
]

PANELS = {"FHMpanel": FHM_PANEL}

# Per panel: whether to highlight AND label only the FDR-significant genes. Off for the 15-gene FHM
# panel, where every gene is worth naming. Worth turning on for any larger panel added later -- once
# a panel runs to a few dozen genes the non-significant points pile up at x~0 and adjustText drags
# connector lines across the whole figure.
SIG_ONLY = {"FHMpanel": False}

# experiment -> the accession prefixing its expression matrix
EXPERIMENTS = {"wnt9b": "9YFNNW", "wnt5b": "2SZ9JF"}

PANEL_COLUMNS = ["gene_name", "baseMean", "log2FoldChange", "LFC_shrunk", "lfcSE", "pvalue", "padj",
                 "gene_biotype"]


def panel_outputs(res, label, tag, panel_name):
    """Panel TSV + volcano for one result table. `res` must already carry an LFC_shrunk column.

    Same plotting convention as the program-9/16 volcanoes in run_deseq2_v2_all_biorep.py --
    everything gray except the panel genes, which are maroon and labeled -- so every figure this
    script writes is directly comparable, whichever panel is highlighted.
    """
    genes = PANELS[panel_name]

    # log2FoldChange is the raw Wald estimate that pvalue/padj correspond to; LFC_shrunk is the one
    # to read for effect size. Sorted by padj, like the other highlighted_* TSVs.
    panel_values = res[res["gene_name"].isin(genes)][PANEL_COLUMNS].sort_values("padj")
    out_tsv = OUTPUT_DIR / f"highlighted_{panel_name}_MEIS2KD_unstim_{tag}.tsv"
    panel_values.to_csv(out_tsv, sep="\t")

    plot_data = res.dropna(subset=["padj", "log2FoldChange"]).copy()
    plot_data["neg_log10_padj"] = -np.log10(plot_data["padj"].clip(lower=1e-300))
    hits = plot_data[plot_data["gene_name"].isin(genes)]

    fig, ax = plt.subplots(figsize=(8, 7))
    ax.scatter(
        plot_data["log2FoldChange"], plot_data["neg_log10_padj"],
        color=NS_COLOR, s=10, alpha=0.6, edgecolor="none", zorder=1,
    )
    ax.axhline(-np.log10(ALPHA), color="black", linewidth=0.75, linestyle="--")
    ax.axvline(0, color="black", linewidth=0.75)
    # SIG_ONLY decides whether non-significant panel genes are highlighted at all. When it is on,
    # they stay in the gray background instead of being coloured -- same convention as the CHD/CVD
    # volcanoes in run_deseq2_v2_all_biorep.py, where a maroon point always means an FDR hit.
    n_sig = int((hits["padj"] < ALPHA).sum())
    shown = hits[hits["padj"] < ALPHA] if SIG_ONLY[panel_name] else hits

    ax.scatter(
        shown["log2FoldChange"], shown["neg_log10_padj"],
        color=PANEL_COLOR, s=40, alpha=0.9, edgecolor="none", zorder=3,
        label=f"{panel_name}, padj < {ALPHA}" if SIG_ONLY[panel_name] else panel_name,
    )
    texts = [
        ax.text(row["log2FoldChange"], row["neg_log10_padj"], row["gene_name"], fontsize=8)
        for _, row in shown.iterrows()
    ]

    ax.set_xlabel(f"log2 fold change (MEIS2 KD vs NT, unstimulated) -- {label}")
    ax.set_ylabel("-log10(padj)")
    title = f"{label}: MEIS2 KD vs NT, unstimulated\n"
    if SIG_ONLY[panel_name]:
        title += (f"{n_sig}/{len(hits)} {panel_name} genes found have padj < {ALPHA}; "
                  f"only these are highlighted")
    else:
        title += f"{len(hits)}/{len(genes)} {panel_name} genes plotted; {n_sig} with padj < {ALPHA}"
    ax.set_title(title)
    ax.legend(loc="upper right", fontsize=7, markerscale=0.7, framealpha=0.9)

    # After set_title -- adjustText needs the final axis limits. Guarded because it errors on empty
    # coordinate arrays.
    if texts:
        adjust_text(
            texts, ax=ax,
            x=shown["log2FoldChange"].to_numpy(), y=shown["neg_log10_padj"].to_numpy(),
            expand=(1.35, 1.6), force_text=(0.5, 0.9), force_static=(0.3, 0.5),
            max_move=None, min_arrow_len=3,
            arrowprops=dict(arrowstyle="-", color="black", lw=0.5),
        )

    out_png = OUTPUT_DIR / f"volcano_MEIS2KD_unstim_{panel_name}_{tag}.png"
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)

    # padj is NA for genes dropped by independent filtering; those cannot go on a -log10(padj) axis,
    # so say which panel genes that removed rather than letting them vanish.
    dropped = sorted(set(genes) - set(hits["gene_name"]))
    print(f"\n{label} / {panel_name}: {len(hits)}/{len(genes)} genes plotted, "
          f"{n_sig} padj < {ALPHA}")
    if dropped:
        print(f"  not plotted (padj NA from independent filtering, or absent): {dropped}")
    print(f"  wrote {out_png}")
    print(f"  wrote {out_tsv}")
    print(panel_values.to_string(float_format=lambda v: f"{v:10.4g}"))
    return panel_values


def load_counts(experiment):
    """Raw counts (samples x genes) and metadata for one experiment, plus its gene annotation.

    Same loading idiom as run_deseq2_v2_all_biorep.py: the matrix is genes x (gene_id, gene_name,
    gene_biotype, then a <sample>_cpm / <sample>_count pair per sample), so keep only the _count
    columns; featureCounts assigns multimappers fractionally so they need rounding to int.
    """
    project_dir = DATA_DIR / f"MEIS2_{experiment}"
    raw = pd.read_csv(project_dir / f"{EXPERIMENTS[experiment]}-expression-matrix.csv.gz")
    gene_info = raw[["gene_id", "gene_name", "gene_biotype"]].set_index("gene_id")

    count_cols = [c for c in raw.columns if c.endswith("_count")]
    counts = raw.set_index("gene_id")[count_cols]
    counts.columns = [c[: -len("_count")] for c in counts.columns]
    counts = counts.round().astype(int).T
    assert not counts.isna().any().any(), "NaNs remain in the count matrix"

    meta = pd.read_csv(project_dir / "metadata.csv", index_col=0).loc[counts.index]
    return counts, meta, gene_info


# Combined fit: unstimulated samples from both experiments, 6 NT vs 6 MEIS2.
print("\n" + "=" * 100)
print("Combined fit: unstimulated samples only, both experiments")

count_blocks, meta_blocks, gene_infos = [], [], []
for experiment in EXPERIMENTS:
    counts, meta, gene_info = load_counts(experiment)

    unstim = meta.index[meta["treatment"] == "none"]
    assert len(unstim) == 6, f"{experiment}: expected 6 unstimulated libraries, got {len(unstim)}"
    counts, meta = counts.loc[unstim], meta.loc[unstim]

    # Both metadata files use the same sample names (1-NT1, 7-M1, ...), so the experiment has to go
    # into the index or the two blocks would silently collide on concat.
    prefixed = [f"{experiment}_{s}" for s in counts.index]
    counts.index, meta.index = prefixed, prefixed
    meta["experiment"] = experiment

    count_blocks.append(counts)
    meta_blocks.append(meta)
    gene_infos.append(gene_info)
    print(f"  {experiment}: {list(counts.index)}")

# Both matrices come from the same annotation; fail loudly rather than silently inner-joining if not.
assert gene_infos[0].index.equals(gene_infos[1].index), "gene_id sets differ between experiments"
gene_info = gene_infos[0]

counts = pd.concat(count_blocks)[gene_infos[0].index]
counts.index.name = "sample_id"
meta = pd.concat(meta_blocks)

# Explicit reference level so genotype[T.MEIS2] is KD-vs-NT and a positive LFC means up in the KD.
meta["genotype"] = pd.Categorical(meta["genotype"], categories=["NT", "MEIS2"])
assert not meta["genotype"].isna().any(), "unexpected genotype levels in metadata.csv"

# The experiment reference is arbitrary -- with no interaction term, swapping it only flips the sign
# of experiment[T.*] and changes what the Intercept means; genotype[T.MEIS2] is identical either way.
# Set explicitly anyway so it cannot change with pandas' alphabetical default (which picks wnt5b).
meta["experiment"] = pd.Categorical(meta["experiment"], categories=list(EXPERIMENTS))

# Unbiased low-count filter, same threshold as the per-experiment scripts.
keep = counts.sum(axis=0) >= 10
counts = counts.loc[:, keep]
print(f"\n  kept {keep.sum()} / {len(keep)} genes (total count >= 10)")
print(f"  {counts.shape} (samples, genes)")
print(meta.groupby(["genotype", "experiment"], observed=True).size())

panels = {}
for design, tag, label in DESIGNS:
    print("\n" + "-" * 100)
    print(f"Fitting {design}")

    dds = DeseqDataSet(
        counts=counts, metadata=meta, design=design,
        refit_cooks=True, inference=DefaultInference(n_cpus=N_CPUS),
    )
    dds.deseq2()
    print(f"\n  design columns: {dds.obsm['design_matrix'].columns.tolist()}")

    ds = DeseqStats(dds, contrast=["genotype", "MEIS2", "NT"], alpha=ALPHA)
    ds.summary()
    res = ds.results_df.join(gene_info)
    res.to_csv(OUTPUT_DIR / f"deseq2_MEIS2KD_unstim_{tag}.csv")

    # lfc_shrink replaces ds.results_df in place, hence taking the unshrunk copy first.
    ds.lfc_shrink(coeff="genotype[T.MEIS2]")
    res_shrunk = ds.results_df.join(gene_info)
    res_shrunk.to_csv(OUTPUT_DIR / f"deseq2_MEIS2KD_unstim_{tag}_shrunk.csv")

    print(f"\n  {(res['padj'] < ALPHA).sum()} genes with padj < {ALPHA} ({design})")
    print(f"  wrote {OUTPUT_DIR / f'deseq2_MEIS2KD_unstim_{tag}.csv'}")
    print(f"  wrote {OUTPUT_DIR / f'deseq2_MEIS2KD_unstim_{tag}_shrunk.csv'}")

    res = res.join(res_shrunk["log2FoldChange"].rename("LFC_shrunk"))
    for panel_name in PANELS:
        panels[(panel_name, tag)] = panel_outputs(res, label, tag, panel_name)

    # MEIS2 itself is not in the panel but is the on-target anchor -- if the pooled KD effect on
    # MEIS2 is not strongly negative, the problem is the pooling, not the panel.
    print("\n  On-target anchor (not part of the panel):")
    print(res[res["gene_name"] == "MEIS2"][PANEL_COLUMNS].to_string(float_format=lambda v: f"{v:10.4g}"))

# Side-by-side: the point of fitting both designs is the lfcSE column, which is where the batch term
# earns its place. Indexed on gene_name because the two tables share gene_ids.
for panel_name, genes in PANELS.items():
    print("\n" + "=" * 100)
    print(f"{panel_name}: no batch term vs + experiment term (LFC / lfcSE / padj)")
    compare = pd.concat(
        {tag: panels[(panel_name, tag)].set_index("gene_name")[["log2FoldChange", "lfcSE", "padj"]]
         for _, tag, _ in DESIGNS},
        axis=1,
    ).reindex(genes).dropna(how="all")
    print(compare.to_string(float_format=lambda v: f"{v:9.3g}"))
