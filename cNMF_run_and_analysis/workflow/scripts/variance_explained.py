import os
import argparse
import numpy as np
import pandas as pd
from scipy import sparse
import sklearn.preprocessing as pp
import anndata as ad

parser = argparse.ArgumentParser()

parser.add_argument( "--tpm", type=str, required=True)
parser.add_argument("--spectra_tpm", type=str, required=True)
parser.add_argument("--usage", type=str, required=True)
parser.add_argument("-K", type=str, required=True)
parser.add_argument("--output", type=str, required=True)
args = parser.parse_args()


def load_df_from_npz(filename):
    with np.load(filename, allow_pickle=True) as f:
        obj = pd.DataFrame(**f)
    return obj


# target
tpm = ad.read_h5ad(args.tpm)

# usage
usage = load_df_from_npz(args.usage)
spectra_tpm = load_df_from_npz(args.spectra_tpm)


recovered_mtx = usage.dot(spectra_tpm)
variance_unexplained = np.subtract(tpm.X.todense(), recovered_mtx).var(axis=0).sum()

tpm_variance = tpm.X.todense().var(axis=0).sum()
variance_explained = 1 - (variance_unexplained / tpm_variance)
print(variance_explained)
output = pd.DataFrame([variance_explained], columns=["VarianceExplained"])
output["K"] = args.K
output["TPMVariance"] = tpm_variance
output["UnexplainedVariance"] = variance_unexplained
output.to_csv(args.output, sep="\t", index=False)
