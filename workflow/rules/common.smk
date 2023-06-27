def get_chroms():
    import re
    chroms = fai["fai"]
    chroms = sorted([chrom for chrom in fai["chr"] if "_" not in chrom])
    chroms = [chrom for chrom in chroms if re.fullmatch(keep_chrs, chrom)]
    return chroms
    