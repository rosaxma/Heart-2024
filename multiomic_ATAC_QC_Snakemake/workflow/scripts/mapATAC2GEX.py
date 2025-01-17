import pandas as pd
import argparse
import pickle
from Bio.Seq import Seq
import os

parser = argparse.ArgumentParser()
parser.add_argument('--barcode_dict', type=str, required=True)
parser.add_argument('--ATACMultipletsBarcodes', type=str, required=True)
parser.add_argument('--ATACDoubletsBarcodes', type=str, required=True)
parser.add_argument('--ATACCellInfoTable', type=str, required=True)
parser.add_argument('--outdir', type=str, required=True)
parser.add_argument('--revcom', type=bool, required=True)
args = parser.parse_args()

barcode_dict_pickle_path=args.barcode_dict
ATAC_multiplet_barcodes=args.ATACMultipletsBarcodes
ATAC_doublet_barcodes=args.ATACDoubletsBarcodes
ATAC_cell_info_table=args.ATACCellInfoTable

barcode_dict_pickle = open(barcode_dict_pickle_path,'rb')
barcode_dict = pickle.load(barcode_dict_pickle )
barcode_dict_pickle.close()

ATAC_cell_info_df=pd.read_csv(ATAC_cell_info_table, delimiter=",", header=0)


ATAC_multiplet_barcodes_df=pd.read_csv(ATAC_multiplet_barcodes, header=None, delimiter=",")
list_multiplets=ATAC_multiplet_barcodes_df[0].tolist()

if os.stat(ATAC_doublet_barcodes).st_size == 0:
	list_doublets=[]
else:
	ATAC_doublet_barcodes_df=pd.read_csv(ATAC_doublet_barcodes, header=None, delimiter=",")
	list_doublets=ATAC_doublet_barcodes_df[0].tolist()

print(ATAC_multiplet_barcodes_df.head(n=10))
if args.revcom:
    ATAC_cell_info_df["GEX_Barcodes"]=ATAC_cell_info_df["barcode"].apply(lambda x: barcode_dict[Seq(x).reverse_complement()])
else:
    ATAC_cell_info_df["GEX_Barcodes"]=ATAC_cell_info_df["barcode"].apply(lambda x: barcode_dict[x])

ATAC_cell_info_df["AmuletMultiplets"]=ATAC_cell_info_df["barcode"].apply(lambda x:  True if (x in list_multiplets) else False)
ATAC_cell_info_df["ArchRDoublets"]=ATAC_cell_info_df["barcode"].apply(lambda x:  True if (x in list_doublets) else False)
ATAC_cell_info_df.to_csv(os.path.join(args.outdir, "ATACBarcode2GEXmapping.tsv"), sep='\t', index=False)
