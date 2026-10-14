#!/usr/bin/env python
"""Step 7e -- combined percent-knockdown bar figure.

Merges the per-gene summary from 07d (asterisks, headline knockdown) with the
per-gRNA table from 07c (the dots), converts each gRNA's log2 fold change to
percent knockdown, and draws bars = mean of that gene's own dots with SEM
whiskers. Positive = expression went down, so the bars point up.

Reads only those two CSVs, so it runs in seconds and can be re-styled without
recomputing anything. This is the published version of the figure (revised
2026-09-15, 6.sceptre/perguideknockdown/20260903-tssknockdown-pct/
07_plot_pct_knockdown_combined_bar.py).

Gene set and column order are set by main-figure panel E, not by the data:
19 columns, in panel E's left-to-right order (``plots.panel_e_groups``), so the
panel stacks column-for-column under that heatmap. ``plots.exclude_genes`` and
``plots.panel_e_groups`` must partition the 27-gene CHD panel exactly. q_bh is
BH-corrected across the whole panel in 07d, so asterisks do not move when the
exclusions change.

Percent knockdown is ``100 * (1 - 2^log2fc)``, a **concave** transform. That
has two consequences this script enforces rather than assumes:

  * the bar must equal the mean of the dots actually plotted, not the percent
    knockdown of the mean log2 fold change -- both are reported so the
    difference is visible;
  * by Jensen's inequality ``mean(%KD) <= %KD(mean log2FC)``, checked below.

The y axis is set from the data and not clipped, so every gRNA that went the
wrong way stays visible.

Six runtime guards: excluded genes must exist, the exclusions and the panel-E
order must partition the panel, every plotted gene must have per-gRNA rows, the
bar must match the mean of its dots, the Jensen inequality must hold, and the
gene-label geometry must match the figure spec.
"""

import argparse
import os

import numpy as np
import pandas as pd

import figstyle as fs  # noqa: F401  -- sets the Agg backend on import
import matplotlib.pyplot as plt  # noqa: E402

from scripts.pipeline_utils import load_config  # noqa: E402

TOL = 1e-9


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", required=True, help="07d per-gene summary CSV")
    parser.add_argument("--per-guide", required=True, help="07c per-gRNA CSV")
    parser.add_argument("--config", required=True)
    parser.add_argument("--output-figure", required=True, help="Figure stem")
    parser.add_argument("--output-source-data", required=True)
    return parser.parse_args()


def pct_knockdown(log2fc):
    return 100.0 * (1.0 - 2.0**np.asarray(log2fc, dtype=float))


def main():
    args = parse_args()
    config = load_config(args.config)
    params = config["per_guide_knockdown"]
    plot_params = params["plots"]

    summary = pd.read_csv(args.summary)
    summary["gene_symbol"] = summary["gene_symbol"].astype(str).str.strip()
    per_guide = pd.read_csv(args.per_guide)
    per_guide["gene_symbol"] = per_guide["gene_symbol"].astype(str).str.strip()

    # --- Guard 1: excluded genes must exist -------------------------------
    exclude = list(plot_params["exclude_genes"])
    missing = [g for g in exclude if g not in set(summary["gene_symbol"])]
    if missing:
        raise ValueError(
            f"plots.exclude_genes names genes absent from the summary: {missing}. "
            "Either the gene set changed or the names are stale."
        )
    n_panel = len(summary)
    summary = summary[~summary["gene_symbol"].isin(exclude)].copy()
    print(f"Gene set: {len(summary)} of the {n_panel} panel genes "
          f"(excluded {', '.join(exclude)})")

    # --- Guard 2: exclusions + panel-E order partition the panel ----------
    groups = plot_params["panel_e_groups"]
    genes = [g for group in groups for g in group["genes"]]
    kept = set(summary["gene_symbol"])
    if len(set(genes)) != len(genes):
        raise ValueError("plots.panel_e_groups contains a duplicate gene")
    if set(genes) != kept:
        raise ValueError(
            "plots.panel_e_groups does not match the kept gene set: "
            f"missing from it {sorted(kept - set(genes))}, "
            f"not in the kept set {sorted(set(genes) - kept)}"
        )
    print(f"Column order: main-figure panel E, {len(genes)} columns  "
          + " | ".join(f"{g['name']}: {', '.join(g['genes'])}" for g in groups))

    # Not filtered on in_main_figure -- that flag selects a different gene set.
    per_guide = per_guide[per_guide["gene_symbol"].isin(kept)].copy()
    per_guide["pct_kd"] = pct_knockdown(per_guide["log2_fold_change"])
    tss = per_guide[per_guide["series"] == "TSS"]
    ntc = per_guide[per_guide["series"] == "NTC"]

    # --- Guard 3: every plotted gene needs dots ---------------------------
    no_dots = sorted(kept - set(tss["gene_symbol"]))
    if no_dots:
        raise ValueError(f"no per-gRNA rows for: {no_dots}")

    lo_q, hi_q = plot_params["ntc_percentiles"]
    by_tss = tss.groupby("gene_symbol")
    by_ntc = ntc.groupby("gene_symbol")["pct_kd"]
    stat = pd.DataFrame({
        "n_tss_guides_plotted": by_tss["pct_kd"].size(),
        "perguide_mean_pct_kd": by_tss["pct_kd"].mean(),
        "perguide_sem_pct_kd": by_tss["pct_kd"].sem(),
        "perguide_mean_log2fc": by_tss["log2_fold_change"].mean(),
    })
    stat["pct_kd_of_mean_log2fc"] = pct_knockdown(stat["perguide_mean_log2fc"])
    stat[f"ntc_p{lo_q:02d}_pct_kd"] = by_ntc.quantile(lo_q / 100)
    stat[f"ntc_p{hi_q:02d}_pct_kd"] = by_ntc.quantile(hi_q / 100)
    # The summary already carries n_ntc_guides; a new name avoids the collision.
    stat["n_ntc_guides_plotted"] = by_ntc.size()

    # Panel-E order, not effect-size order: reindex, because the sequence is
    # editorial. Guard 2 guarantees every name is present.
    df = summary.set_index("gene_symbol").join(stat, how="inner").reindex(genes)

    # --- Guard 4: bar height == mean of the plotted dots ------------------
    dot_mean = tss.groupby("gene_symbol")["pct_kd"].mean().reindex(df.index)
    diff = (df["perguide_mean_pct_kd"] - dot_mean).abs().max()
    if diff > TOL:
        raise AssertionError(f"bar height != mean of plotted dots (max diff {diff:.3g})")

    # --- Guard 5: Jensen's inequality -------------------------------------
    bad = df.index[df["perguide_mean_pct_kd"] > df["pct_kd_of_mean_log2fc"] + TOL]
    if len(bad):
        raise AssertionError(
            f"mean(%KD) > %KD(mean log2FC) for {list(bad)}, which a concave "
            "transform cannot do -- check the inputs."
        )

    # --- Figure -----------------------------------------------------------
    # Font size and jitter seed are the published figure's.
    fs.set_style(base_font_size=12)
    rng = np.random.RandomState(0)

    n = len(df)
    xs = np.arange(n)
    fig, ax = plt.subplots(figsize=(0.50 * n + 2.1, 4.9))
    ax.axhline(0, color=fs.INK_SOFT, lw=0.8, zorder=2)
    ax.bar(xs, df["perguide_mean_pct_kd"], width=0.66, color=fs.BAR_GREY,
           edgecolor="none", zorder=3)
    ax.errorbar(xs, df["perguide_mean_pct_kd"], yerr=df["perguide_sem_pct_kd"],
                fmt="none", ecolor=fs.DOT_INK, elinewidth=1.0, capsize=2.2,
                capthick=1.0, zorder=5)
    for i, gene in enumerate(df.index):
        vals = tss.loc[tss["gene_symbol"] == gene, "pct_kd"].to_numpy()
        ax.scatter(i + rng.uniform(-0.19, 0.19, len(vals)), vals, s=7,
                   color=fs.DOT_INK, alpha=0.75, linewidths=0.3,
                   edgecolors=fs.SURFACE, zorder=6)

    # Significance, in a row above the tallest dot.
    lo = float(tss["pct_kd"].min())
    top = float(tss["pct_kd"].max())
    span = top - lo
    for i, q in enumerate(df["q_bh"]):
        mark = fs.stars(q)
        ax.text(i, top + 0.025 * span, mark if mark else "n.s.", ha="center",
                va="bottom", fontsize=10 if mark else 8,
                color=fs.INK if mark else fs.INK_MUTED)

    ax.set_xlim(-0.75, n - 0.25)
    ax.set_ylim(lo - 0.06 * span, top + 0.14 * span)
    ax.set_xticks(xs)
    ax.set_xticklabels(df.index)
    ax.set_ylabel("Knockdown (%)")
    ax.tick_params(axis="x", length=5, width=1.0, pad=3)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_position(("outward", 2))

    # Every text property of the gene labels is applied here, after all other
    # x-axis calls, and none as set_xticklabels kwargs: in this matplotlib those
    # kwargs are silently dropped and the rcParams values win.
    label_pt = plot_params["gene_label_pt"]
    plt.setp(ax.get_xticklabels(), rotation=90, ha="center", va="top",
             fontsize=label_pt, fontstyle="italic", color=fs.INK)
    # --- Guard 6: label geometry is what the figure spec asks for --------
    got = {(t.get_rotation(), t.get_fontsize()) for t in ax.get_xticklabels()}
    if got != {(90.0, float(label_pt))}:
        raise AssertionError(f"gene labels not (90 deg, {label_pt}pt): {got}")

    os.makedirs(os.path.dirname(os.path.abspath(args.output_figure)), exist_ok=True)
    fs.save(fig, args.output_figure)
    print(f"Wrote {args.output_figure}.svg / .png")

    cols = ["n_tss_guides_plotted", "n_ntc_guides_plotted",
            "perguide_mean_pct_kd", "perguide_sem_pct_kd",
            "pct_kd_of_mean_log2fc", "pct_knockdown",
            "perguide_mean_log2fc", "log2fc_tss_vs_ntc_guides",
            f"ntc_p{lo_q:02d}_pct_kd", f"ntc_p{hi_q:02d}_pct_kd",
            "p_mannwhitney", "q_bh"]
    out = df[cols].copy()
    out.insert(0, "gene_symbol", out.index)
    os.makedirs(
        os.path.dirname(os.path.abspath(args.output_source_data)), exist_ok=True
    )
    out.to_csv(args.output_source_data, index=False)
    print(f"Wrote {args.output_source_data}")

    fdr = plot_params["fdr_threshold"]
    print(f"\n{int((df['q_bh'] < fdr).sum())} / {n} genes significant at q < {fdr}")
    print(f"Bars (per-gRNA mean %KD) {df['perguide_mean_pct_kd'].min():.1f}% .. "
          f"{df['perguide_mean_pct_kd'].max():.1f}%; dots span {lo:.1f}% .. "
          f"{top:.1f}% ({int((tss['pct_kd'] < 0).sum())}/{len(tss)} gRNAs below 0)")
    gap = df["pct_kd_of_mean_log2fc"] - df["perguide_mean_pct_kd"]
    print(f"Jensen gap %KD(mean log2FC) - mean(%KD): {gap.min():.2f} .. "
          f"{gap.max():.2f} points (largest: {gap.idxmax()})")
    print("All guards passed.")


if __name__ == "__main__":
    main()
