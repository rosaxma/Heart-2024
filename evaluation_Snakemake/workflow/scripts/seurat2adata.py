import pandas as pd
import argparse
import os
import sys

parser = argparse.ArgumentParser()
parser.add_argument("--seurat", type=str, default="", required=True),
parser.add_argument(
    "--output_file",
    type=str,
    required=True,
),
parser.add_argument(
    "--rlibPath",
    type=str,
    required=True,
),

args = parser.parse_args()

os.environ["R_LIBS"] = args.rlibPath
import scib
import scanpy as sc

adata = scib.preprocessing.read_seurat(args.seurat)
print(adata.X.shape)
adata.write(args.output_file)
