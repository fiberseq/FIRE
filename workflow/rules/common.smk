def get_chroms():
    import re

    chroms = fai["fai"]
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
