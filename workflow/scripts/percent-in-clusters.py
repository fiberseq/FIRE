import pandas as pd
import sys

def my_groupby(gdf):
    real=gdf[gdf.case == "Real"]
    null=gdf[gdf.case == "Null"]
    return real.length.sum() - null.length.sum()

bed=sys.argv[1]
df = pd.read_csv(bed, sep="\t", header=None)
df.columns = c["ct", "st", "en", "cov", "case"]
df["length"] = df.en - df.st
print(df)
df = df.groupby("cov").apply(my_groupby).reset_index()
print(df)
