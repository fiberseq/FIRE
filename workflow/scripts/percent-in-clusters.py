import pandas as pd
import sys

def my_groupby(gdf):
    real=gdf[gdf[5] == "Real"]
    null=gdf[gdf[5] == "Null"]
    return real.length.sum() - null.length.sum()

bed=sys.argv[1]
df = pd.read_csv(bed, sep="\t", header=None)
df["length"] = df[2] - df[1]
print(df)
df = df.groupby([4]).apply(my_groupby).reset_index()
print(df)
