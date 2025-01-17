import pandas as pd
import argparse
import pickle
from Bio.Seq import Seq
import os

parser = argparse.ArgumentParser()
parser.add_argument('--barcode_dict', type=str, required=True)
parser.add_argument('--GEX_barcodes_df', type=str, required=True)
parser.add_argument('--outdir', type=str, required=True)
parser.add_argument('--revcom', type=bool, required=True)
parser.add_argument('--sample', type=str, required=True)
args = parser.parse_args()

barcode_dict_pickle_path=args.barcode_dict

barcode_dict_pickle = open(barcode_dict_pickle_path,'rb')
barcode_dict = pickle.load(barcode_dict_pickle )
barcode_dict_pickle.close()


barcode_dict_gex={v: k for k, v in barcode_dict.items()}

GexInfo_df =pd.read_csv(args.GEX_barcodes_df, delimiter="\t", header=0)
print(GexInfo_df.shape)
print(GexInfo_df.head().to_string())
GexInfo_df = GexInfo_df[GexInfo_df.EffectiveCells == True]
print(GexInfo_df.shape)
print(GexInfo_df.head().to_string())
output_df = GexInfo_df.loc[:, ["cbc"]]
print(output_df.head().to_string())
if args.revcom:
     output_df["ATAC_Barcodes"]=output_df["cbc"].apply(lambda x: Seq(barcode_dict_gex[x]).reverse_complement())
else:
     output_df["ATAC_Barcodes"]=output_df["cbc"].apply(lambda x: barcode_dict_gex[x])

output_df["ATAC_Barcodes_achr"] = args.sample + "#" + output_df['ATAC_Barcodes']


output_df.to_csv(os.path.join(args.outdir, "ATAC_barcodes_for_archr.tsv"), sep='\t', index=False)
