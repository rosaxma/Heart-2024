import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os
import random
import pandas as pd

random.seed(42)

parser = argparse.ArgumentParser()
parser.add_argument('--mtx', type=str, required=True)
parser.add_argument('--gene_tsv', type=str, required=True)
parser.add_argument('--barcodes', type=str, required=True)
parser.add_argument('--sample', type=str, required=True)
parser.add_argument('--outdir', type=str, required=True)
parser.add_argument('--threshold', type=float, required=True)
args = parser.parse_args()

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42

mtx=args.mtx
gene_tsv=args.gene_tsv
sample=args.sample
counts_matrix = scipy.io.mmread(mtx).T.tocsc()
genes = np.array(scr.load_genes(gene_tsv, delimiter='\t', column=0))
print(counts_matrix.shape)
#cells as rows and genes as columns.
print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
print('Number of genes in gene list: {}'.format(len(genes)))
#######################################################
nCell=counts_matrix.shape[0]
#https://docs.google.com/spreadsheets/d/17UXVWsSNuMrF_RMg_ewBGmtEACm3lwmYu2wxPB7IYZo/edit?usp=sharing
mulitplet_rate = (nCell*0.000779-0.0667)/100
####################################################### 
scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=mulitplet_rate)
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_gene_variability_pctl=85, n_prin_comps=30)
scrub.call_doublets(threshold=args.threshold)
#######################################################
scrub.plot_histogram()
plt.suptitle(args.sample)
plt.savefig(args.outdir + "/" + sample +"_line_plot.pdf", format="pdf")
print('Running UMAP...')
"""
n_neighbors: This determines the number of neighboring points used in local approximations of manifold structure. Larger values will result in more global structure being preserved at the loss of detailed local structure. In general this parameter should often be in the range 5 to 50, with a choice of 10 to 15 being a sensible default.
min_dist: This controls how tightly the embedding is allowed compress points together. Larger values ensure embedded points are more evenly distributed, while smaller values allow the algorithm to optimise more accurately with regard to local structure. Sensible values are in the range 0.001 to 0.5, with 0.1 being a reasonable default.
metric: This determines the choice of metric used to measure distance in the input space. A wide variety of metrics are already coded, and a user defined function can be passed as long as it has been JITd by numba
"""
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
scrub.plot_embedding('UMAP', order_points=True)
plt.suptitle(args.sample)
plt.savefig(args.outdir + "/" + sample + "_Doublet_UMAP.pdf", format="pdf")
#######################################################
df = pd.DataFrame({
    'doublet_score': scrub.doublet_scores_obs_,
    'predicted_doublet': scrub.predicted_doublets_
})

barcodes=pd.read_csv(args.barcodes,compression='gzip',sep="\t", header=None)
print(str(barcodes))
df["cell"]=barcodes[0]
print(df.shape)
df.to_csv(args.outdir + "/"+ sample + "_scrublet_output_table.csv", index=False)
