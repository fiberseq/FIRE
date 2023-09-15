import re
import logging


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


def get_mem_mb_xl(wildcards, attempt):
    return attempt * 1024 * 92


def get_mem_mb_small(wildcards, attempt):
    return attempt * 1024 * 4


def get_load(wc):
    if "all" in wc.sm:
        return 100
    return 50


def find_median_coverage(file, outfile=None, min_out=None, max_out=None):
    if force_coverage is not None:
        coverage = force_coverage
    else:
        df = pd.read_csv(
            file, sep="\t", header=None, names=["chr", "start", "end", "coverage"]
        )
        logging.info(f"Calculating median coverage from {file}\n{df}")
        df = df[df.coverage > 0]
        df = df[df["chr"].isin(get_chroms())]
        total = (df.end - df.start).sum()
        coverage = (df.coverage * (df.end - df.start)).sum() / total

    min_coverage = get_min_coverage(coverage)
    max_coverage = get_max_coverage(coverage)

    if coverage <= 1:
        raise ValueError(
            f"Median coverage is {coverage}! Did you use the correct reference, or is data missing from most of your genome. If so consider the keep_chromosomes parameter in config.yaml"
        )

    if outfile is not None:
        open(outfile, "w").write(str(round(coverage)) + "\n")
    if min_out is not None:
        open(min_out, "w").write(str(round(min_coverage)) + "\n")
    if max_out is not None:
        open(max_out, "w").write(str(round(max_coverage)) + "\n")
    return round(coverage)


def get_median_coverage(wc):
    if force_coverage is not None:
        return force_coverage
    median_coverages = expand(rules.genome_bedgraph.output.median, sm=wc.sm)
    return find_median_coverage(median_coverages[0])


def get_min_coverage(median):
    sd = math.sqrt(median)
    mmin = median - coverage_within_n_sd * sd
    return max(mmin, min_coverage)


def get_max_coverage(median):
    sd = math.sqrt(median)
    return median + coverage_within_n_sd * sd


def grep_command_for_el_type(wc):
    if wc.el_type == "nucleosome":
        return "(rg '230,230,230' || true)"
    elif wc.el_type == "linker":
        return f"(rg -v '230,230,230' || true) | awk '$10>{min_fire_fdr}'"
    elif wc.el_type == "fire":
        return f"awk '$10<={min_fire_fdr}'"
    else:
        raise ValueError(f"Unknown element type {wc.el_type}")


def hap_grep_term(wc):
    if wc.hp == "all":
        return '""'
    elif wc.hp == "hap1":
        return "H1"
    elif wc.hp == "hap2":
        return "H2"
    else:
        raise ValueError(f"Unknown haplotype {wc.hp}")


def hap_hck_columns(wc):
    if wc.hp == "all":
        return "-F fire_coverage -F coverage"
    elif wc.hp == "hap1":
        return "-F fire_coverage_H1 -F coverage_H1"
    elif wc.hp == "hap2":
        return "-F fire_coverage_H2 -F coverage_H2"
    else:
        raise ValueError(f"Unknown haplotype {wc.hp}")
