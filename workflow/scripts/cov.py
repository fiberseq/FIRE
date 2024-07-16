import pandas as pd
import numpy as np
import sys
import os
import math
import polars as pl

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
    # group by value and sum the weights
    gdf = df.groupby(val)[weight].sum().reset_index().sort_values(val)
    print(gdf, file=sys.stderr)

    gdf["cumsum"] = gdf[weight].cumsum()
    gdf["cutoff"] = gdf[weight].sum() / 2.0
    print(gdf, file=sys.stderr)
    comparison = gdf[gdf["cumsum"] >= gdf["cutoff"]][val]
    # print(comparison, file=sys.stderr)
    return comparison.iloc[0]


def pandas_read():
    df = pd.read_csv(
        snakemake.input.bg,
        sep="\t",
        header=None,
        names=["chr", "start", "end", "coverage"],
    )
    df = df.loc[df["coverage"] > 0]
    df = df.loc[df["chr"].isin(snakemake.params.chroms)]
    df["weight"] = df["end"] - df["start"]
    return df


def polars_read():
    # Reading the CSV file using the lazy API
    df = pl.scan_csv(
        snakemake.input.bg,
        separator="\t",
        has_header=False,
        #new_columns=["chr", "start", "end", "coverage"],
    )
    # Applying filters and adding the new column using the lazy API
    df = (
        df.filter(pl.col("coverage") > 0)
        .filter(pl.col("chr").is_in(snakemake.params.chroms))
        .with_column((pl.col("end") - pl.col("start")).alias("weight"))
    )
    return df.collect().to_pandas() 

df = polars_read()
print(df, file=sys.stderr)
coverage = weighted_median(df, "coverage", "weight")

min_coverage = get_min_coverage(coverage)
max_coverage = get_max_coverage(coverage)
mean = (df["coverage"] * df["weight"]).sum() / df["weight"].sum()
print(f"\nmean coverage: {mean}", file=sys.stderr)
print(f"median coverage: {coverage}\n", file=sys.stderr)

if coverage <= 1:
    raise ValueError(
        f"Median coverage is {coverage}! Did you use the correct reference, or is data missing from most of your genome. If so consider the keep_chromosomes parameter in config.yaml"
    )
open(snakemake.output.cov, "w").write(str(round(coverage)) + "\n")
open(snakemake.output.minimum, "w").write(str(round(min_coverage)) + "\n")
open(snakemake.output.maximum, "w").write(str(round(max_coverage)) + "\n")
