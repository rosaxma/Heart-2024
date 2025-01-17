import pandas as pd
import argparse
import pickle
from Bio.Seq import Seq
import os

parser = argparse.ArgumentParser()
parser.add_argument('--barcode_dict', type=str, required=True)
parser.add_argument('--GEX_barcodes', type=str, required=True)
parser.add_argument('--outdir', type=str, required=True)
parser.add_argument('--revcom', type=bool, required=True)
parser.add_argument('--sample', type=str, required=True)
parser.add_argument("--frag", type=str, required=True)
args = parser.parse_args()

barcode_dict_pickle_path=args.barcode_dict

barcode_dict_pickle = open(barcode_dict_pickle_path,'rb')
barcode_dict = pickle.load(barcode_dict_pickle)
barcode_dict_pickle.close()


barcode_dict_gex={v: k for k, v in barcode_dict.items()}

GexInfo_df =pd.read_csv(args.GEX_barcodes, delimiter="\t", header=None)
fragment_file=pd.read_csv(args.frag, delimiter="\t", header=None)

print(GexInfo_df.head())
print(fragment_file.head())

GexInfo_df_keep =GexInfo_df.iloc[:, [0]]
GexInfo_df_keep.columns=["GEX_cbc"]
print(GexInfo_df_keep.shape)
ATAC_barcode_keep=fragment_file.iloc[:, [3]].drop_duplicates()
ATAC_barcode_keep.columns=["ATAC_cbc"]
print(GexInfo_df_keep.head())
output_df = GexInfo_df_keep.loc[:, ["GEX_cbc"]]

print(ATAC_barcode_keep.shape)
print(ATAC_barcode_keep.head())
if args.revcom:
     output_df["ATAC_cbc"]=output_df["GEX_cbc"].apply(lambda x: str(Seq(barcode_dict_gex[x]).reverse_complement()))
else:
     output_df["ATAC_cbc"]=output_df["GEX_cbc"].apply(lambda x: barcode_dict_gex[x])
print(output_df.shape)
output_df = output_df.merge(ATAC_barcode_keep, how="inner", on="ATAC_cbc")
print(output_df.shape)
print(output_df.head())
output_df.to_csv(os.path.join(args.outdir, "retained_GEX_cells_for_archR.tsv"), sep='\t', index=False)
