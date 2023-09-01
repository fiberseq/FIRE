import re


def get_chroms():
    chroms = fai["chr"]
    chroms = sorted([chrom for chrom in fai["chr"] if "chrUn_" not in chrom])
    chroms = [chrom for chrom in chroms if "_random" not in chrom]
    chroms = [chrom for chrom in chroms if re.fullmatch(keep_chrs, chrom)]
    if len(chroms) == 0:
        raise ValueError(
            f"No chromosomes left after filtering. Check your keep_chromosomes parameter in config.yaml. "
            f"Your fai file contains the following chromosomes: {fai['chr']}"
        )
    return chroms


def get_mem_mb(wildcards, attempt):
    if attempt < 3:
        return attempt * 1024 * 32
    return attempt * 1024 * 48


def get_large_mem_mb(wildcards, attempt):
    return attempt * 1024 * 64


def get_mem_mb_small(wildcards, attempt):
    return attempt * 1024 * 4


def get_load(wc):
    if "all" in wc.sm:
        return 100
    return 50


def find_median_coverage(file, outfile=None):
    if force_coverage is not None:
        coverage = force_coverage
    else:
        df = pd.read_csv(
            file, sep="\t", header=None, names=["chr", "start", "end", "coverage"]
        )
        df = df[df.coverage > 0]
        df = df[df["chr"].isin(get_chroms())]
        total = (df.end - df.start).sum()
        coverage = (df.coverage * (df.end - df.start)).sum() / total

    if coverage <= 1:
        raise ValueError(
            f"Median coverage is {coverage}! Did you use the correct reference, or is data missing from most of your genome. If so consider the keep_chromosomes parameter in config.yaml"
        )

    if outfile is not None:
        open(outfile, "w").write(str(round(coverage)) + "\n")
    return round(coverage)


def get_median_coverage(wc):
    if force_coverage is not None:
        return force_coverage
    median_coverages = expand(rules.genome_bedgraph.output.median, sm=wc.sm)
    print(median_coverages)
    return find_median_coverage(median_coverages)


def get_min_coverage(wc):
    median = get_median_coverage(wc)
    sd = math.sqrt(median)
    mmin = median - coverage_within_n_sd * sd
    return max(mmin, min_coverage)


def get_max_coverage(wc):
    median = get_median_coverage(wc)
    sd = math.sqrt(median)
    return median + coverage_within_n_sd * sd
