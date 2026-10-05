# Methods — developmental expression of the VIC-nominated CHD genes (Supplementary Fig. 24)

## Gene set

We examined the 23 CHD causal genes assigned to valve interstitial cells (VICs)
in Fig. 2: CHD causal genes with an expression z-score > 1 across cell types and
> 1 TPM in at least one of the five VIC subpopulations. For mouse we used the 1:1
orthologs (checked against MGI nomenclature). Genes are shown in the same order
in every panel, by mean expression in human VICs at 4–6 post-conception weeks
(pcw).

## Expression summary used in all panels

For each cell, expression was log1p(CP10K), the natural log of
1 + count / total UMI of that cell × 10,000. The total was taken over every
gene in the dataset's matrix. For each group of cells (a column in the figure)
we report two numbers for each gene:
- the mean log1p(CP10K) across the cells;
- the percentage of cells with at least one count.

Cell-type annotations and developmental stages are those provided by the
authors of each dataset. No cells were re-clustered and no values were z-scored.

## Mouse endocardium to AV mesenchyme (Supplementary Fig. 24a)

We used the atrioventricular (AV) valve development atlas of Lotto et al.
(Nat. Commun. 2023), downloaded from the Broad Single Cell Portal (SCP1897). It
combines whole hearts at E7.75, E8.25 and E9.25 from de Soysa et al. (2019) with
AV canals dissected at E10.5 and E12.5, under a single annotation by the
original authors.

**Normalisation.** The portal's log-normalised matrix is per-gene scaled
(Seurat `ScaleData`) and so reports relative expression. We therefore computed
log1p(CP10K) from the raw count matrix (21,862 genes × 54,561 cells), with each
cell's total summed over all genes (median 20,159 UMI per cell).

**Cells and columns.** We used the authors' endocardial-to-mesenchymal
transition (EndMT) subset (10,354 cells), with its cell-state annotation, stage
and pseudotime. We kept the five states on the path from endocardium to AV
mesenchyme: endocardium, AV endocardium, early AV mesenchyme, proliferating
mesenchyme and AV mesenchyme. They are ordered by their mean pseudotime (0.13,
0.31, 0.80, 0.85 and 0.95). Early endothelium, which precedes endocardium, and
valvular endothelial cells I and II, which remain endothelial, were not shown.
Each state is shown at every stage where it has at least 50 cells, 8,317 cells
in all. Five groups with 1–43 cells were omitted.

**Checks.**
- A separate implementation, reading the raw count matrix independently,
  reproduced every value exactly.
- Thinning every cell to 10,000 UMI changed the mean log1p(CP10K) by at most
  0.14 in any column.

## Human valve interstitial cells (Supplementary Fig. 24b)

We used the developing human heart atlas of Bayraktar, Cranley et al. (bioRxiv
2024.04.27.591127): 297,473 cells and nuclei × 36,601 genes from 21 hearts at
4–20 pcw. The log-normalised object (`RNA_global_lognormalised.h5ad`) was
downloaded from the Heart Cell Atlas foetal portal (CC BY 4.0).

**Normalisation.** Its expression values are the authors' log1p(CP10K): the
per-cell sum of expm1(values) is 10,000. They were used as provided.

**Cells and columns.** VICs are the authors' `ValveInterstitialCells`
annotation at the finest of their three levels. They are summarised for each of
the 12 post-conception weeks with VIC data (4–10, 12–15 and 20 pcw;
31–2,577 cells per week, 7,398 in total).

**Earliest weeks.** The 4–6 pcw group, shaded in the figure and used to order
the genes, is 2,853 cells from 5 hearts. The 4 and 5 pcw data are single-nucleus
profiles from one heart each; the 6 pcw data are whole-cell profiles from three
hearts. At these stages the valves are endocardial cushions, so these cells are
cushion mesenchyme, the precursors of VICs.

## VICs of this atlas by post-conception week (Supplementary Fig. 24c)

We used the raw gene-expression counts of all 13,636 cells assigned to the five
VIC subpopulations (VIC_1–VIC_5) of our fetal heart atlas, from 40 donors.

**Linking cells to age.** Each cell was linked to its donor and that donor's
post-conception week through the atlas cell annotation (sample identifier plus
ATAC barcode). Every cell joined uniquely to an atlas cell of the same
subpopulation. Each cell's total UMI matched the atlas record exactly, and the
pooled counts reproduced the atlas pseudobulk TPM (r > 0.999 for each
subpopulation).

**Columns.** The five subpopulations were pooled and summarised per week as
above. 6–11 pcw are shown individually (555–4,693 cells and 2–8 donors per
week). 12–22 pcw have one or two donors each, so they were pooled into one group
(848 cells, 12 donors).

**Checks.** The age trends held within each subpopulation analysed separately,
and after every cell was thinned to 1,500 UMI.

The 23-gene set was defined on these same subpopulations, so this panel is the
endpoint of the gene selection and not an independent test of it.

## Visualisation

Dot colour encodes mean log1p(CP10K) and dot area the percentage of cells with
non-zero counts. Colour scales use a monotone power-law ramp (γ = 0.55), ticked
in absolute units, so the low-expression range is legible without changing rank
or magnitude.
- **Panels a and b** share one scale from 0 to 2.6; values above the cap are
  flagged on the colour bar.
- **Panel c** uses the same ramp from 0 to 3.4. The atlas is sequenced about ten
  times less deeply per cell (median 2,314 UMI), which raises log1p(CP10K) for
  detected genes, so panel c colours are not comparable with panels a and b.

These are descriptive expression summaries; no statistical tests were performed.

## Software and availability

Analyses used Python 3.10 with numpy 2.2, pandas 2.3, scipy, anndata and
matplotlib 3.10. Code, the processed inputs and scripts that rebuild every table
from its source are available at [repository URL / DOI to be inserted]. Public
datasets: Broad Single Cell Portal SCP1897 (mouse) and the Heart Cell Atlas
foetal portal (human).
