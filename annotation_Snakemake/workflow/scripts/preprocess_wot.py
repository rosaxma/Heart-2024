import os
import argparse
import wot
import os
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--h5ad", type=str, required=True)
parser.add_argument("--hvg", type=str, required=True)
parser.add_argument("--geneset", type=str, required=True)
parser.add_argument("--outdir", type=str, required=True)
parser.add_argument("--sample", type=str, required=True)
args = parser.parse_args()
##############################################################################################################
def convert2gmx(df, save_path, name):
    all_lists=[]
    for i in df.iloc[:,0].unique():
        sub_list=["cluster_"+str(i), ""]
        subset_df=cell_annotation[cell_annotation.iloc[:,0]==i]
        sub_list.extend(subset_df.index)
        all_lists.append(sub_list)
    df = pd.DataFrame(all_lists).transpose()
    df.to_csv(os.path.join(save_path, name+".gmx"), sep="\t", header=False, index=False)
##############################################################################################################
# subset for hvg h5ad
adata = wot.io.read_dataset(args.h5ad)
hvg=pd.read_csv(args.hvg, header=None)
adata=adata[:, hvg[0].tolist()]
new_filename=args.h5ad.replace(".h5ad", ".var.gene.h5ad")
adata.write(os.path.join(new_filename))
##############################################################################################################
# calculate gene sets scores
gs= wot.io.read_sets(args.geneset, adata.var.index.values)
gene_set_scores_df = pd.DataFrame(index=adata.obs.index)
for j in range(gs.shape[1]):
    gene_set_name = str(gs.var.index.values[j])
    if gene_set_name in ["Cell.cycle", "Apoptosis"]:
        result = wot.score_gene_sets(ds=adata, gs=gs[:, [j]], permutations=0, method='mean_z_score')
        gene_set_scores_df[gene_set_name] = result['score']
gene_set_scores_df.to_csv(os.path.join(args.outdir, "proliferation_apop_geneset_scores.csv"), index_label="id")
##############################################################################################################
#create cell days
adata.obs["day"]=adata.obs["Rounded.PCW"]
df_days=adata.obs["day"].to_frame().reset_index()
df_days=df_days.rename(columns={"index":"id"})
df_days=df_days.set_index("id")
df_days.to_csv(os.path.join(args.outdir, "cell_days.txt"), sep="\t", header=True, index_label="id")
##############################################################################################################
# create cell sets
feature="annotation"
cell_annotation=adata.obs[feature].to_frame()
cell_annotation_gmx=convert2gmx(cell_annotation,args.outdir, feature)

proliferation=gene_set_scores_df["Cell.cycle"]
apoptosis = gene_set_scores_df['Apoptosis']
def logistic(x, L, k, x0=0):
    f = L / (1 + np.exp(-k * (x - x0)))
    return f
def gen_logistic(p, beta_max, beta_min, pmax, pmin, center, width):
    return beta_min + logistic(p, L=beta_max - beta_min, k=4 / width, x0=center)

def beta(p, beta_max=1.7, beta_min=0.3, pmax=1.0, pmin=-0.5, center=0.25):
    return gen_logistic(p, beta_max, beta_min, pmax, pmin, center, width=0.5)

def delta(a, delta_max=1.7, delta_min=0.3, amax=0.5, amin=-0.4, center=0.1):
    return gen_logistic(a, delta_max, delta_min, amax, amin, center,
                          width=0.2)

birth = beta(proliferation)
death = delta(apoptosis)

# growth rate is given by 
gr = np.exp(birth-death)
growth_rates_df = pd.DataFrame(index=gene_set_scores_df.index, data={'cell_growth_rate':gr})
growth_rates_df.to_csv(os.path.join(args.outdir, "growth_gs_init.txt"), index_label="id")
