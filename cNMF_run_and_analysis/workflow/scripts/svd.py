import os
import argparse

parser = argparse.ArgumentParser()

parser = argparse.ArgumentParser()
parser.add_argument("--threads", type=str, required=True)
parser.add_argument("--tpm_mtx", type=str, required=True)
parser.add_argument("--k", type=int, required=True)
parser.add_argument("--outdir", type=str, required=True)
args = parser.parse_args()

os.environ["OMP_NUM_THREADS"] = args.threads

import numpy as np
from scipy.sparse.linalg import svds
from scipy.sparse import issparse
from scipy.io import mmread
import sklearn.preprocessing as pp
import pickle
import matplotlib.pyplot as plt
import anndata as ad

tpm = ad.read_h5ad(args.tpm_mtx)

mtx = tpm.X
if issparse(mtx):
    norm_mtx = pp.scale(mtx, axis=0, with_mean=False, copy=True)
else:
    norm_mtx = mtx / mtx.std(axis=0, ddof=1)


k = args.k
s = svds(norm_mtx, k=k, which="LM", return_singular_vectors=False, random_state=42)

with open(os.path.join(args.outdir, "singular_values.pickle"), "wb") as handle:
    pickle.dump(s, handle)


print(s.flatten().shape)
s_sorted = np.sort(s)[::-1]
index = np.array(range(1, k + 1))
print(s_sorted)
fig = plt.figure(figsize=(20, 10))
plt.scatter(index, s_sorted)
plt.xticks(range(1, k + 1))

plt.savefig(os.path.join(args.outdir, "singular_value.pdf"))
s_sorted.tofile(os.path.join(args.outdir, "singular_values.csv"), sep=",")
