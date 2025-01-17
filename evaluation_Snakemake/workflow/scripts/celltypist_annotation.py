import os
import argparse

parser = argparse.ArgumentParser()

parser = argparse.ArgumentParser()
parser.add_argument("--trained_model", type=str, required=True)
parser.add_argument("--model_name", type=str, required=True)
parser.add_argument("--mtx", type=str, required=True)
parser.add_argument("--gene_list", type=str, required=True)
parser.add_argument("--barcode_list", type=str, required=True)
parser.add_argument("--transpose_input", type=bool, required=True)
parser.add_argument("--outdir", type=str, required=True)
args = parser.parse_args()
print(args)
import celltypist
from celltypist import models
from scipy.io import mmread, mmwrite

models.models_path= os.path.dirname(args.trained_model)
model=models.Model.load(os.path.basename(args.trained_model))
print(args.transpose_input)

predictions=celltypist.annotate(args.mtx, gene_file=args.gene_list, cell_file=args.barcode_list, model=model, majority_voting=True, transpose_input=args.transpose_input)

df = predictions.predicted_labels
df.to_csv(os.path.join(args.outdir, args.model_name+"_annotation.tsv"))
