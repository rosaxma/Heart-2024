from kneed import KneeLocator
import pandas as pd

import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--betaScoreTable", type=str, required=True)
parser.add_argument("--k", type=int, required=True)
parser.add_argument("--outdir", type=str, required=True)
args = parser.parse_args()


beta_df = pd.read_table(args.betaScoreTable, sep="\t", index_col=0, skiprows=0)
k = args.k
col_names = ["Knee_X"]
knee_df = pd.DataFrame(columns=col_names)

for i in range(1, k + 1):
    print(i)
    colname = "GEP_" + str(i) + "_beta_score"
    sub_df = beta_df.loc[:, [colname]]
    sub_df = sub_df[sub_df[colname] > 0]
    sub_df["rank"] = sub_df[colname].rank(
        ascending=False, method="first", na_option="bottom"
    )
    try:
        kn = KneeLocator(
            sub_df["rank"],
            sub_df[colname],
            S=5,
            curve="convex",
            direction="decreasing",
            online=True,
            interp_method="polynomial",
        )
        knee_df.loc[colname, "Knee_X"] = kn.knee
    except IndexError:
        knee_df.loc[colname, "Knee_X"] = len(sub_df[colname])

knee_df.to_csv(os.path.join(args.outdir, "cNMF_GEP_knee_K" + str(k) + ".tsv"), sep="\t")
