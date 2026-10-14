# Figure legend — Supplementary Fig. 24

## Full version

**Supplementary Figure 24. Expression of the 23 VIC-nominated CHD genes across
mouse and human valve-lineage development.**
Rows are the 23 CHD causal genes assigned to valve interstitial cells (VICs) in
Fig. 2 on the basis of cell-type-specific expression, in the same order in all
panels: by mean expression in human VICs at 4–6 post-conception weeks (pcw),
highest at top. Mouse genes are the 1:1 orthologs. In every panel, dot colour is
the mean log1p(CP10K) of a gene across the cells in a column, and dot area is
the percentage of those cells with at least one count. Cell types and
developmental stages are those annotated by the authors of each dataset; no
cells were re-clustered and no values are z-scored.
(**a**) Mouse endocardium-to-atrioventricular (AV) mesenchyme transition in the
valve-development atlas of Lotto et al. (2023; Broad Single Cell Portal
SCP1897). Columns are the authors' endocardial-to-mesenchymal transition (EndMT)
cell states, in order of their mean pseudotime, at each embryonic stage where
the state has at least 50 cells:
- endocardium: E7.75, E8.25, E9.25 and E12.5 (n = 260, 220, 417 and 153 cells);
- AV endocardium: E8.25 to E12.5 (n = 143, 762, 345 and 213);
- early AV mesenchyme: E10.5 (n = 2,408);
- proliferating mesenchyme: E12.5 (n = 987);
- AV mesenchyme: E12.5 (n = 2,409).

That is 8,317 cells in total. E7.75–E9.25 are whole hearts from de Soysa et al.
(2019); E10.5 and E12.5 are dissected AV canals from Lotto et al.
(**b**) Human VICs (annotation `ValveInterstitialCells`) in the developing heart
atlas of Bayraktar, Cranley et al. (2024), at the 12 pcw with VIC data
(4–10, 12–15 and 20 pcw; 31–2,577 cells per week, 7,398 in total). Shading marks
4–6 pcw (2,853 cells from 5 hearts), when the valves are still endocardial
cushions. 4 and 5 pcw are single-nucleus data from one heart each; 6 pcw is
whole-cell data from three hearts.
(**c**) VICs of the fetal heart atlas presented here (subpopulations VIC_1–VIC_5
pooled; 13,636 cells from 40 donors) by pcw. 6–11 pcw are shown individually
(555–4,693 cells, 2–8 donors per week). 12–22 pcw are pooled (848 cells from 12
donors) because each week has one or two donors.
Panels a and b share one colour scale, capped at 2.6, with values above the cap
marked on the colour bar. Panel c has its own scale (0–3.4), because the atlas is
sequenced about ten times less deeply per cell (median 2,314 UMI), which raises
log1p(CP10K) for detected genes. Colours are therefore comparable within panel
c but not with a and b, and dot areas are comparable within a panel only. The
23-gene set was defined on the population in panel c, so that panel is the
endpoint of the gene selection rather than an independent test of it.
Expression supports cell-type assignment but does not by itself establish a
causal role in valve development.

## Short version

**Supplementary Figure 24. Expression of the 23 VIC-nominated CHD genes across
mouse and human valve-lineage development.** Rows: the 23 CHD genes assigned to
VICs in Fig. 2, ordered by mean expression in human 4–6 pcw VICs. Dot colour,
mean log1p(CP10K); dot area, % of cells expressing; author-provided annotations
throughout.
(**a**) Mouse endocardium → AV mesenchyme (Lotto et al. 2023): EndMT cell states
in pseudotime order, at each stage with ≥ 50 cells (8,317 cells).
(**b**) Human VICs by pcw (Bayraktar, Cranley et al. 2024; 7,398 cells), 4–6 pcw
shaded.
(**c**) VICs of this atlas by pcw, 12–22 pcw pooled (13,636 cells, 40 donors).
Panels a and b share a colour scale; panel c has its own, because its sequencing
depth is lower.

---

## Notes for the text (not for the legend)

- **Panel a: "activated during EndMT".** The cleanest comparison is E12.5 AV
  mesenchyme vs E12.5 endocardium, the same stage. The largest gains are GPC3
  +1.94, TBX20 +1.78, MEIS2 +1.66 and GATA6 +1.63 (log1p CP10K). PRKD1, INVS and
  LTBP2 barely change (+0.02, +0.02, +0.01), so "many of these genes" is the
  right hedge.
- **Panel a: mesenchyme starts at E10.5.** Only 18 early AV mesenchyme cells
  are annotated at E9.25, below the 50-cell threshold, so "E10.5–E12.5" matches
  the panel.
- **Panel b: 4–6 pcw rests on 5 hearts and two modalities.** 4 and 5 pcw are
  one heart each, both single-nucleus; 6 pcw is three hearts, whole cells. Some
  of the week-to-week differences may reflect modality rather than age. At these
  stages the cells are cushion mesenchyme (VIC precursors), so "VICs/cushion
  mesenchyme at 4–6 pcw" is the precise wording.
- **Panel c: most genes decline with age.** Expression is highest at 6–7 pcw
  for most of the 23 (e.g. MEIS2 3.18 at 6 pcw vs 2.39 at 12–22 pcw). LTBP2
  (0.11 → 0.50) and FBN1 (1.63 → 1.93) rise. All stay detected throughout. The
  decline holds within each VIC subtype and after depth-matching, so it is not a
  change in subtype mix or sequencing depth.
- **LTBP2 is the clearest exception** to early expression: near-absent in
  mouse and in human 4–6 pcw VICs, but rising later in both human datasets. It
  behaves like a late valve-maturation ECM gene.
