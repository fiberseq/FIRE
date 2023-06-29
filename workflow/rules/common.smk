import re


def get_chroms():
    chroms = fai["chr"]
    chroms = sorted([chrom for chrom in fai["chr"] if "_" not in chrom])
    chroms = [chrom for chrom in chroms if re.fullmatch(keep_chrs, chrom)]
    return chroms


def get_mem_mb(wildcards, attempt):
    if attempt < 3:
        return attempt * 1024 * 32
    return attempt * 1024 * 48


def get_large_mem_mb(wildcards, attempt):
    return attempt * 1024 * 64


def get_load(wc):
    if "all" in wc.sm:
        return 100
    return 50


def find_median_coverage(file, outfile=None):
    df = pd.read_csv(
        file, sep="\t", header=None, names=["chr", "start", "end", "coverage"]
    )
    df = df[df.coverage > 0]
    df = df[df["chr"].isin(get_chroms())]
    total = (df.end - df.start).sum()
    coverage = (df.coverage * (df.end - df.start)).sum() / total
    if outfile is not None:
        open(outfile, "w").write(str(round(coverage)) + "\n")
    return round(coverage)
