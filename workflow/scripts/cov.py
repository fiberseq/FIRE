import pandas as pd
import sys
import os
import math

coverage_within_n_sd = snakemake.params.coverage_within_n_sd
MIN_coverage = snakemake.params.min_coverage


def get_min_coverage(median):
    sd = math.sqrt(median)
    mmin = median - coverage_within_n_sd * sd
    return max(mmin, MIN_coverage)


def get_max_coverage(median):
    sd = math.sqrt(median)
    return median + coverage_within_n_sd * sd


def weighted_median(df, val, weight):
    df.sort_values(val, inplace=True)
    cumsum = df[weight].cumsum()
    cutoff = df[weight].sum() / 2.0
    return df[cumsum >= cutoff][val].iloc[0]


df = pd.read_csv(
    snakemake.input.bg,
    sep="\t",
    header=None,
    names=["chr", "start", "end", "coverage"],
)
df = df.loc[df["coverage"] > 0]
df = df[df["chr"].isin(snakemake.params.chroms)]
df = df[~df["chr"].isin(["chrX", "chrY", "chrM", "chrEBV"])]
df["weight"] = df["end"] - df["start"]
print(df, file=sys.stderr)
coverage = weighted_median(df, "coverage", "weight")

min_coverage = get_min_coverage(coverage)
max_coverage = get_max_coverage(coverage)
print(coverage, file=sys.stderr)
if coverage <= 1:
    raise ValueError(
        f"Median coverage is {coverage}! Did you use the correct reference, or is data missing from most of your genome. If so consider the keep_chromosomes parameter in config.yaml"
    )
open(snakemake.output.cov, "w").write(str(round(coverage)) + "\n")
open(snakemake.output.minimum, "w").write(str(round(min_coverage)) + "\n")
open(snakemake.output.maximum, "w").write(str(round(max_coverage)) + "\n")
