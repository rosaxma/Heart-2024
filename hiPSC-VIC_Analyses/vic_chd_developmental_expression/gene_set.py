"""Single source of truth for the 23 VIC-nominated CHD genes (Fig 2).

Loads the gene set and the 10-gene "expression-only" subset (genes assigned to
VICs by expression specificity alone -- no variant-to-gene, cNMF program-27, or
knockdown support) from vic_nominated_chd_genes_fig2.tsv, so the list lives in
exactly one place. Every figure/analysis script imports from here instead of
re-declaring the genes by hand.
"""
from pathlib import Path
import csv

TSV = Path(__file__).resolve().parent / "vic_nominated_chd_genes_fig2.tsv"
OTHER_EVIDENCE_COLUMNS = (
    "V2G_VIC_nominated",
    "in_VIC_program27",
    "perturbseq_KD_target",
    "bulk_KD_target",
)


def to_mouse(human_symbol):
    """Human -> mouse 1:1 ortholog (Title-case). Valid for all 23 genes in this set."""
    return human_symbol[0] + human_symbol[1:].lower()


def load_gene_set(tsv=TSV):
    """Return (genes, expression_only) as human symbols, read from the canonical TSV."""
    genes, expression_only = [], []
    with open(tsv) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            gene = row["Gene"]
            genes.append(gene)
            if not any(row[col] == "TRUE" for col in OTHER_EVIDENCE_COLUMNS):
                expression_only.append(gene)
    return genes, expression_only


GENES, EXPRESSION_ONLY_LIST = load_gene_set()          # human symbols (23, 10)
EXPRESSION_ONLY = set(EXPRESSION_ONLY_LIST)
GENES_MOUSE = [to_mouse(g) for g in GENES]
EXPRESSION_ONLY_MOUSE = {to_mouse(g) for g in EXPRESSION_ONLY}
MOUSE_TO_HUMAN = {to_mouse(g): g for g in GENES}       # explicit ortholog join key
