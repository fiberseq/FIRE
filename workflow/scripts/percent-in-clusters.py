import pandas as pd
import sys

def my_groupby(gdf):
    real=gdf[gdf.case == "Real"]
    null=gdf[gdf.case == "Null"]
    return (real.length*real.cov).sum() - (null.length*real.cov).sum()

bed=sys.argv[1]
df = pd.read_csv(bed, sep="\t", header=None)
df.columns = ["ct", "st", "en", "cov", "case"]
df["length"] = df.en - df.st
real_bp = (df[df.case == "Real"].length * df[df.case == "Real"].cov).sum()
print(df)
df = df.groupby("cov").apply(my_groupby).reset_index()
print(df)
over_expected = df[0][df[0]>0].sum()
print(f"{over_expected/real_bp:%}\n")
