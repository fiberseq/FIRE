import pandas as pd
import sys
import numpy as np


def my_groupby(gdf):
    real = gdf[gdf.case == "Real"]
    null = gdf[gdf.case == "Null"]
    return (real.length * real["cov"]).sum() - (null.length * null["cov"]).sum()


bed = sys.argv[1]
df = pd.read_csv(bed, sep="\t", header=None)
df.columns = ["ct", "st", "en", "cov", "case"]
df["length"] = df.en - df.st
real_bp = (df[df.case == "Real"]["length"] * df[df.case == "Real"]["cov"]).sum()
df = df.groupby("cov").apply(my_groupby).reset_index()
over_expected = df[0][df[0] > 0].sum()
print(f"{over_expected/real_bp:%}\n")

# chr1    0       1       12      1       1       1
cov = sys.argv[2]
cov = pd.read_csv(cov, sep="\t", header=None)
cov.columns = ["ct", "st", "en", "fdr", "acc", "link", "nuc"]
n_tests = cov.shape[0]
min_fdr = -10.0 * np.log10(0.01 / n_tests)
n_peaks = cov[cov.fdr >= min_fdr].shape[0]


rtn = f"""percent-of-MSPs-preferentially-clustered-along-the-genome\tmin_fdr
{over_expected/real_bp:%}\t{min_fdr}
"""
out = sys.argv[3]
open(out, "w").write(rtn)
