import pandas as pd


def weighted_median(df, val, weight):
    df_sorted = df.sort_values(val)
    cumsum = df_sorted[weight].cumsum()
    cutoff = df_sorted[weight].sum() / 2.0
    return df_sorted[cumsum >= cutoff][val].iloc[0]


df = pd.read_csv(
    snakemake.input.median,
    sep="\t",
    header=None,
    names=["chr", "start", "end", "coverage"],
)
df = df[df.coverage > 0]
df = df[df["chr"].isin(get_chroms())]
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
open(output.cov, "w").write(str(round(coverage)) + "\n")
open(output.minimum, "w").write(str(round(min_coverage)) + "\n")
open(output.maximum, "w").write(str(round(max_coverage)) + "\n")
