import re
import logging
import sys


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
