"""Combined figure: the 23 VIC-nominated CHD genes across valve-lineage development.

Supports this sentence in the text:

  "Finally, we also examined the expression levels of these 23 genes across time
   in other transcriptomics datasets in mouse and human heart development. Many
   of these genes are indeed activated during the endocardial-to-mesenchymal
   transition (a), are expressed and generally up-regulated in mouse mesenchymal
   cells from E10.5-E12.5, and are expressed in human valvular interstitial
   cells/mesenchymal at earlier timepoints (p.c.w. 4-6)."

Self-contained: panel a reads tables/lotto_endmt_celltype_timepoint.tsv (built by
build_lotto_endmt_table.py from data/lotto_scp1897/); panel b reads the Cranley
tables in data/; panel c reads data/fhm_vic_by_pcw.tsv (built on Sherlock by
build_fhm_vic_by_pcw_table.py from the atlas's per-cell VIC counts). See README.md.

Three panels, one shared set of 23 gene rows:
  a  mouse endocardium -> AV mesenchyme (Lotto 2023, SCP1897): dot plot per
     EndMT cell type x timepoint, using the authors' endo_clusterfile labels,
     at the timepoints where each cell type has >= MIN_CELLS cells (the
     "activated during EndMT" and "up-regulated in mesenchyme E10.5-E12.5"
     claims). Replaces the 2026-09-15 panels a (row-scaled pseudotime) and
     b (AV mesenchyme I+II pooled); that version is archived in
     figures/archive_4panel_20260915/.
  b  human VIC by post-conception week, pcw4-6 shaded      (the "expressed at
     earlier timepoints" claim). Dataset = Bayraktar, Cranley et al. 2024,
     bioRxiv 2024.04.27.591127; NOT the Development 2025 paper the upstream
     METHODS_NOTES cites -- see FIGURE_LEGEND.md.
  c  our fetal heart atlas, VICs (VIC_1-5 pooled) by post-conception week,
     6W-11W plus 12-22W pooled (1-2 donors per week there). Replaces the
     2026-09-28 panel c (pseudobulk TPM per VIC subtype, ages pooled).

Substantive differences from the upstream version
-------------------------------------------------
1. Panel b shows ALL 12 available post-conception weeks.
   Upstream showed 6 hand-picked weeks, dropping 7W (n=2,577 VICs, the
   best-powered week in the atlas) while keeping 9W (n=182) and 20W (n=112):
   half of all VIC cells were not shown, and the thinnest weeks were.
2. Panels a and b share one absolute log1p(CP10K) colour scale with a monotone
   power-law ramp (gamma<1) ticked in absolute units, so the low end is
   readable. The upstream shared 0-3.5 linear cap was set by human MEIS2.
3. Panel c uses the same metrics and dot-size key as a and b but its own
   colour scale (0-3.4): the atlas is ~10x shallower (median ~2,200 UMI per
   cell), which raises log1p(CP10K) of detected genes, so on the shared 2.6
   cap its highest genes (MEIS2, TBX20) would saturate and hide the change
   with age. Colours are comparable within panel c, not with a/b.

Outputs: figures/combined_valve_lineage_expression_revised.{svg,png,pdf}

NOTE ON FONTS: Arial is requested, but no Arial (or metric-compatible
Liberation Sans / Arimo) is installed on Sherlock. Because svg.fonttype is
"none", the SVG carries font-family "Arial" as live text and renders in real
Arial when opened anywhere that has it -- that is the file to use. The PNG and
PDF written here fall back to DejaVu Sans for rasterisation.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import Normalize, PowerNorm, LogNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

HERE = Path(__file__).resolve().parent
TABLES = HERE / "data"
OUT = HERE / "figures"
OUT.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(HERE))
from gene_set import MOUSE_TO_HUMAN, EXPRESSION_ONLY, TSV as GENESET_TSV  # noqa: E402

# Set True to bring back the 4-column V2G / cNMF-27 / Perturb-seq-KD / bulk-KD
# evidence matrix at the left (an empty row = one of the 10 expression-only
# genes). Off by default: the sentence this figure supports is about expression
# across development, not about the evidence breakdown.
SHOW_EVIDENCE = False

# ---------------------------------------------------------------- style -----
INK = "#1A1A1A"
INK_SOFT = "#5A5F66"
RULE = "#D3D8DD"
POS = "#B2182B"
NEG = "#2166AC"
ACCENT = "#2166AC"
WARN = "#B2182B"

plt.rcParams.update({
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",          # SVG keeps live text, so it picks up real Arial
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "Arimo", "DejaVu Sans"],
    "mathtext.fontset": "custom",
    "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold",
    "font.size": 8,
    "axes.labelcolor": INK,
    "axes.edgecolor": INK_SOFT,
    "axes.linewidth": 0.6,
    "text.color": INK,
    "xtick.color": INK_SOFT,
    "ytick.color": INK_SOFT,
    "legend.frameon": False,
    "figure.dpi": 150,
    "savefig.dpi": 400,
    "savefig.transparent": False,
})

CMAP_EXPR = "viridis"    # absolute log1p(CP10K), panels a and b
VMAX_EXPR = 2.6          # shared cap for log1p(CP10K); over-colour flags saturation
VMAX_C = 3.4             # panel c's own cap (shallower data; see docstring note 3)
GAMMA = 0.55             # power-law norm: monotone, expands the low end
MIN_CELLS = 50           # panel a: drop cell type x timepoint groups smaller than this

# Panel a: the endocardium -> AV mesenchyme path through the authors' EndMT
# labels (endo_clusterfile `clusters`), in mean-pseudotime order (asserted
# below). Early endothelium (upstream of endocardium) and VECs I/II (valvular
# endothelium, which stays endothelial) are in the same subset but are not on
# this path.
ENDMT_PATH = ["endocardium", "AV endocardium", "early AV mesenchyme",
              "proliferating mesenchyme", "AV mesenchyme"]
ENDMT_LABELS = {"endocardium": "Endocardium", "AV endocardium": "AV\nendocardium",
                "early AV mesenchyme": "Early AV\nmesenchyme",
                "proliferating mesenchyme": "Proliferating\nmesenchyme",
                "AV mesenchyme": "AV\nmesenchyme"}
TP_MOUSE = ["E7.75", "E8.25", "E9.25", "E10.5", "E12.5"]
NORM_EXPR = PowerNorm(GAMMA, vmin=0, vmax=VMAX_EXPR)
NORM_C = PowerNorm(GAMMA, vmin=0, vmax=VMAX_C)
DOT_BASE, DOT_SPAN = 8.0, 200.0   # scatter area: s = DOT_BASE + pct/100 * DOT_SPAN

# ------------------------------------------------------------------ data ----
early = pd.read_csv(TABLES / "taskC_cranley_early_pcw4to6_summary.tsv", sep="\t")
early = early[early.cell_type == "ValveInterstitialCells"]
GENES = early.sort_values("mean_lognorm", ascending=False, kind="stable")["gene"].tolist()
assert len(GENES) == 23, f"expected 23 genes, got {len(GENES)}"
yorder = list(reversed(GENES))                 # bottom -> top, so GENES[0] is on top
ypos = {g: i for i, g in enumerate(yorder)}
NG = len(yorder)

endmt = pd.read_csv(HERE / "tables" / "lotto_endmt_celltype_timepoint.tsv", sep="\t")
endmt["gene_h"] = endmt["gene"].map(MOUSE_TO_HUMAN)
assert endmt["gene_h"].notna().all(), "unmapped mouse symbol in EndMT table"
assert set(endmt.gene_h) == set(GENES), "panel a gene set mismatch"
_pt = endmt.drop_duplicates("cell_type").set_index("cell_type")["cell_type_mean_pseudotime"]
assert list(_pt.reindex(ENDMT_PATH).sort_values().index) == ENDMT_PATH, \
    "ENDMT_PATH is not in mean-pseudotime order"
groups = endmt.drop_duplicates(["cell_type", "timepoint"])[["cell_type", "timepoint", "n_cells"]]
groups = groups[groups.cell_type.isin(ENDMT_PATH)]
dropped = groups[groups.n_cells < MIN_CELLS]
groups = groups[groups.n_cells >= MIN_CELLS].assign(
    _c=lambda d: d.cell_type.map(ENDMT_PATH.index), _t=lambda d: d.timepoint.map(TP_MOUSE.index)
).sort_values(["_c", "_t"]).reset_index(drop=True)
A_COLS = list(zip(groups.cell_type, groups.timepoint))
xa = {k: i for i, k in enumerate(A_COLS)}
endmt = endmt[[k in xa for k in zip(endmt.cell_type, endmt.timepoint)]]
assert len(endmt) == len(A_COLS) * 23

cr = pd.read_csv(TABLES / "taskC_cranley_expr_by_celltype_week_tidy.tsv", sep="\t")
cr = cr[cr.cell_type == "ValveInterstitialCells"].copy()
WEEK_ORDER = ["4W", "5W", "6W", "7W", "8W", "9W", "10W", "12W", "13W", "14W", "15W", "20W"]
WEEKS = [w for w in WEEK_ORDER if w in set(cr.week)]
assert set(cr.week) == set(WEEKS), f"unexpected weeks: {set(cr.week) - set(WEEKS)}"
N_EARLY = 3                                    # 4W/5W/6W = the pcw4-6 block

fhm = pd.read_csv(TABLES / "fhm_vic_by_pcw.tsv", sep="\t")
assert set(fhm.gene) == set(GENES), "panel c gene set mismatch"
C_BINS = (fhm.drop_duplicates("bin").assign(_w=lambda d: d["bin"].str.extract(r"^(\d+)")[0].astype(int))
          .sort_values("_w")["bin"].tolist())                   # 6W ... 11W, 12-22W
assert len(fhm) == len(C_BINS) * 23 and fhm.mean_lognorm.max() <= VMAX_C
C_META = fhm.drop_duplicates("bin").set_index("bin").loc[C_BINS,
                                                         ["n_cells", "n_donors", "weeks"]]

ev = pd.read_csv(GENESET_TSV, sep="\t").set_index("Gene")
EV_COLS = ["V2G_VIC_nominated", "in_VIC_program27", "perturbseq_KD_target", "bulk_KD_target"]
EV_LABELS = ["V2G", "cNMF 27", "Perturb-seq KD", "bulk KD"]


# ---------------------------------------------------------------- layout ----
# dot panels a and b get the same width per column, so dots line up in size
ratios = ([0.52] if SHOW_EVIDENCE else []) + [0.232 * (len(A_COLS) + 0.3),
                                              0.232 * (12 + 0.3),
                                              0.232 * (len(C_BINS) + 0.3)]
fig = plt.figure(figsize=(11.4 if SHOW_EVIDENCE else 10.8, 8.2))
outer = gridspec.GridSpec(
    1, len(ratios), figure=fig,
    width_ratios=ratios, wspace=0.14,
    left=0.078, right=0.885, top=0.815, bottom=0.115,
)
col = {name: i for i, name in enumerate((["ev"] if SHOW_EVIDENCE else [])
                                        + ["a", "b", "c"])}
ax_ev = fig.add_subplot(outer[0, col["ev"]]) if SHOW_EVIDENCE else None
ax_a = fig.add_subplot(outer[0, col["a"]])
ax_b = fig.add_subplot(outer[0, col["b"]])
ax_c = fig.add_subplot(outer[0, col["c"]])


def bare(ax, keep=()):
    for name, spine in ax.spines.items():
        spine.set_visible(name in keep)
    ax.tick_params(length=0, pad=2.5)


def rows(ax, labels=False, bands=False):
    ax.set_ylim(-0.6, NG - 0.4)
    ax.set_yticks(range(NG))
    ax.set_yticklabels(yorder if labels else [], fontsize=7.8, color=INK)
    if bands:                     # faint row guides so the eye tracks across gaps
        for i in range(NG):
            ax.axhline(i, color=RULE, lw=0.35, zorder=0)


def title_above(ax, letter, text):
    """Panel title anchored in figure space, so titles can never collide with
    each other. The panel letter is real bold text rather than mathtext, so
    every glyph in the file stays in one font."""
    x0 = ax.get_position().x0
    if letter:
        fig.text(x0, 0.832, letter, fontsize=8.8, color=INK, ha="left", va="bottom",
                 fontweight="bold")
        x0 += 0.011
    fig.text(x0, 0.832, text, fontsize=8.8, color=INK, ha="left", va="bottom",
             linespacing=1.5)


# ------------------------------------------- evidence matrix (optional) -----
if SHOW_EVIDENCE:
    for g in yorder:
        for j, c_ in enumerate(EV_COLS):
            on = str(ev.loc[g, c_]).strip().upper() == "TRUE"
            ax_ev.add_patch(Rectangle((j - 0.33, ypos[g] - 0.33), 0.66, 0.66,
                                      facecolor=ACCENT if on else "#FFFFFF",
                                      edgecolor=ACCENT if on else RULE, linewidth=0.5))
    ax_ev.set_xlim(-0.62, len(EV_COLS) - 0.38)
    rows(ax_ev, labels=True)
    ax_ev.set_xticks(range(len(EV_COLS)))
    ax_ev.set_xticklabels(EV_LABELS, rotation=90, fontsize=6.3, color=INK_SOFT)
    bare(ax_ev)
    title_above(ax_ev, "", "Evidence\nroutes")

# ------------------------------------ panel a: mouse endocardium -> AV mes --
rows(ax_a, labels=not SHOW_EVIDENCE, bands=True)
for _, r in endmt.iterrows():
    ax_a.scatter(xa[(r.cell_type, r.timepoint)], ypos[r.gene_h],
                 s=DOT_BASE + r.pct_expressing / 100 * DOT_SPAN,
                 c=[r.mean_lognorm], cmap=CMAP_EXPR, norm=NORM_EXPR,
                 edgecolors="#FFFFFF", linewidths=0.35, zorder=3)
ax_a.set_xlim(-0.65, len(A_COLS) - 0.35)
ax_a.set_xticks(range(len(A_COLS)))
ax_a.set_xticklabels([t for _, t in A_COLS], rotation=90, fontsize=7.4, color=INK)
bare(ax_a)
# cell-type groups: a rule between groups, and a bracket + label under the ticks
_bounds = [(ct, min(i for i, k in enumerate(A_COLS) if k[0] == ct),
            max(i for i, k in enumerate(A_COLS) if k[0] == ct)) for ct in ENDMT_PATH
           if any(k[0] == ct for k in A_COLS)]
for _, i0, i1 in _bounds[1:]:
    ax_a.axvline(i0 - 0.5, color=RULE, lw=0.8, zorder=1)
_tr = ax_a.get_xaxis_transform()
for ct, i0, i1 in _bounds:
    ax_a.plot([i0 - 0.32, i1 + 0.32], [-0.075, -0.075], transform=_tr, color=INK_SOFT,
              lw=0.8, clip_on=False, solid_capstyle="butt")
    ax_a.text((i0 + i1) / 2, -0.088, ENDMT_LABELS[ct], transform=_tr, rotation=90,
              ha="center", va="top", fontsize=7.4, color=INK, linespacing=1.0,
              multialignment="right")
title_above(ax_a, "a", "Mouse endocardium → AV mesenchyme\n(Lotto et al. 2023)")

# ------------------------------------------------- panel b: human VIC -------
xc = {w: i for i, w in enumerate(WEEKS)}
ax_b.axvspan(-0.5, N_EARLY - 0.5, color="#FBE8CC", zorder=0, linewidth=0)
rows(ax_b, bands=True)
for _, r in cr.iterrows():
    ax_b.scatter(xc[r.week], ypos[r.gene],
                 s=DOT_BASE + r.pct_expressing / 100 * DOT_SPAN,
                 c=[r.mean_lognorm], cmap=CMAP_EXPR, norm=NORM_EXPR,
                 edgecolors="#FFFFFF", linewidths=0.35, zorder=3)
ax_b.set_xlim(-0.65, len(WEEKS) - 0.35)
ax_b.set_xticks(range(len(WEEKS)))
ax_b.set_xticklabels(WEEKS, rotation=90, fontsize=7.4, color=INK)
bare(ax_b)
ax_b.text((N_EARLY - 1) / 2, NG - 0.3, "pcw4–6", ha="center", va="bottom", fontsize=7.0,
          color="#9A6A1E")
title_above(ax_b, "b", "Human VIC\n(Bayraktar & Cranley 2024)")

# ------------------------------------ panel c: our atlas VIC by pcw --------
xd = {b: i for i, b in enumerate(C_BINS)}
rows(ax_c, bands=True)
for _, r in fhm.iterrows():
    ax_c.scatter(xd[r["bin"]], ypos[r.gene],
                 s=DOT_BASE + r.pct_expressing / 100 * DOT_SPAN,
                 c=[r.mean_lognorm], cmap=CMAP_EXPR, norm=NORM_C,
                 edgecolors="#FFFFFF", linewidths=0.35, zorder=3)
ax_c.set_xlim(-0.65, len(C_BINS) - 0.35)
ax_c.set_xticks(range(len(C_BINS)))
ax_c.set_xticklabels(C_BINS, rotation=90, fontsize=7.4, color=INK)
bare(ax_c)
title_above(ax_c, "c", "Our atlas VIC\n(fetal heart map)")

# ---------------------------------------------------------------- legends ---
cax = fig.add_axes([0.905, 0.560, 0.0105, 0.200])
cb = fig.colorbar(plt.cm.ScalarMappable(norm=NORM_EXPR, cmap=CMAP_EXPR), cax=cax,
                  extend="max", ticks=[0, 0.25, 0.5, 1.0, 1.5, 2.0, 2.5])
cb.set_label("mean log1p(CP10K)\npanels a, b", fontsize=7.0, color=INK, labelpad=4)
cb.ax.tick_params(labelsize=6.4)
cb.outline.set_linewidth(0)

cax2 = fig.add_axes([0.905, 0.330, 0.0105, 0.160])
cb2 = fig.colorbar(plt.cm.ScalarMappable(norm=NORM_C, cmap=CMAP_EXPR), cax=cax2,
                   ticks=[0, 0.25, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
cb2.set_label("mean log1p(CP10K)\npanel c", fontsize=7.0, color=INK, labelpad=4)
cb2.ax.tick_params(labelsize=6.4)
cb2.outline.set_linewidth(0)

lax = fig.add_axes([0.899, 0.130, 0.095, 0.130])
lax.axis("off")
lax.legend(handles=[Line2D([0], [0], marker="o", color="none", markerfacecolor="#7C8792",
                           markeredgecolor="none",
                           markersize=np.sqrt(DOT_BASE + p / 100 * DOT_SPAN), label=f"{p}%")
                    for p in [25, 50, 75, 100]],
           title="% cells expressing\n(panels a–c)", loc="upper left",
           labelspacing=1.0, fontsize=6.6, title_fontsize=7.0, borderpad=0.2,
           handletextpad=0.9)

# ------------------------------------------------------------ annotation ----
fig.suptitle("The 23 VIC-nominated CHD genes are activated during EndMT and "
             "sustained through the valve-interstitial-cell stage",
             fontsize=11.5, x=0.078, ha="left", y=0.945)
fig.text(0.078, 0.905,
         "Gene rows shared across panels, ordered by human pcw4–6 VIC expression "
         "(highest at top).",
         fontsize=7.6, color=INK_SOFT, ha="left")

base = OUT / "combined_valve_lineage_expression_revised"
for ext in ("svg", "png", "pdf"):
    fig.savefig(f"{base}.{ext}", bbox_inches="tight")
plt.close(fig)

print("wrote", base, "(.svg/.png/.pdf)")
print(f"expression norm PowerNorm(gamma={GAMMA}, 0, {VMAX_EXPR}) for panels a, b; "
      f"panel c PowerNorm(gamma={GAMMA}, 0, {VMAX_C}), data max {fhm.mean_lognorm.max():.2f}")
print("\npanel c columns (n cells, n donors, weeks):")
for b, r in C_META.iterrows():
    print(f"  {b:7s} n={r.n_cells:5d}  donors={r.n_donors:2d}  weeks {r.weeks}")
print("\ngene order (top -> bottom):", GENES)
print(f"\npanel a columns (cell type, timepoint, n), MIN_CELLS={MIN_CELLS}:")
for _, r in groups.iterrows():
    print(f"  {r.cell_type:26s} {r.timepoint:6s} n={r.n_cells}")
print("dropped (< MIN_CELLS): " + ", ".join(
    f"{r.cell_type} {r.timepoint} (n={r.n_cells})" for _, r in dropped.iterrows()))
over = endmt[endmt.mean_lognorm > VMAX_EXPR]
print(f"panel a max log1p(CP10K) {endmt.mean_lognorm.max():.2f}; "
      f"{len(over)} dots above the {VMAX_EXPR} cap"
      + "".join(f"; {r.gene_h} {r.cell_type} {r.timepoint}" for _, r in over.iterrows()))

# Same-stage lineage contrast (no confound with developmental time):
# E12.5 AV mesenchyme minus E12.5 endocardium, absolute log1p(CP10K).
_m = endmt.set_index(["gene_h", "cell_type", "timepoint"])["mean_lognorm"]
print("\nE12.5 AV mesenchyme - E12.5 endocardium (log1p CP10K):")
for g in GENES:
    d = _m[(g, "AV mesenchyme", "E12.5")] - _m[(g, "endocardium", "E12.5")]
    print(f"  {g:8s} {d:+.2f}{'   [expression-only]' if g in EXPRESSION_ONLY else ''}")
