import os
import argparse

parser = argparse.ArgumentParser()

parser = argparse.ArgumentParser()
parser.add_argument("--mtx", type=str, required=True)
parser.add_argument("--filtered_barcode_gz", type=str, required=True)
parser.add_argument("--filtered_gene_gz", type=str, required=True)
parser.add_argument("--output_h5ad", type=str, required=True)
args = parser.parse_args()

import anndata as ad
import scipy.sparse as sparse
from scipy.io import mmread
import pandas as pd

mtx = mmread(args.mtx)
mtx = sparse.csr_matrix(mtx.transpose())
obs = pd.read_csv(args.filtered_barcode_gz, header=None)
obs.columns = ["cbc"]
obs = obs.set_index("cbc", drop=False)
obs.columns = ["cellbarcode"]
var = pd.read_csv(args.filtered_gene_gz, header=None, index_col=0)
var.index.names = ["Gene"]
adata = ad.AnnData(mtx, obs=obs, var=var)
adata.write_h5ad(args.output_h5ad)
