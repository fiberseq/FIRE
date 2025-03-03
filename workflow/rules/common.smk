import re
import logging
import sys

FIRST_REPORT = True


def get_ref():
    if "ref" not in config:
        raise ValueError("FIRE: ref parameter is missing in config.yaml")
    ref = config["ref"]
    if not os.path.exists(ref):
        raise ValueError(f"FIRE: reference file {ref} does not exist")
    return os.path.abspath(ref)


def get_fai():
    fai = f"{get_ref()}.fai"
    if not os.path.exists(fai):
        raise ValueError(f"FIRE: reference index file {fai} does not exist")
    return fai


def get_excludes():
    excludes = config.get("excludes", [])
    if REF_NAME == "hg38" or REF_NAME == "GRCh38":
        files = [
            "../annotations/hg38.gap.bed.gz",
            "../annotations/hg38.blacklist.ENCFF356LFX.bed.gz",
            "../annotations/SDs.merged.hg38.bed.gz",
        ]
        excludes += [workflow.source_path(file) for file in files]
    return excludes


def get_fai_df():
    fai = get_fai()
    return pd.read_csv(fai, sep="\t", names=["chr", "length", "x", "y", "z"])


def get_chroms():
    global FIRST_REPORT
    min_contig_length = config.get("min_contig_length", 0)
    skipped_contigs = FAI_DF["chr"][FAI_DF["length"] < min_contig_length]
    if len(skipped_contigs) > 0 and FIRST_REPORT:
        print(
            f"WARNING: Skipping contigs with length < {min_contig_length:,}: {skipped_contigs}",
            file=sys.stderr,
        )

    chroms = FAI_DF["chr"][FAI_DF["length"] >= min_contig_length]
    chroms = sorted([chrom for chrom in chroms if "chrUn_" not in chrom])
    chroms = [chrom for chrom in chroms if "_random" not in chrom]
    chroms = [chrom for chrom in chroms if re.fullmatch(KEEP_CHRS, chrom)]

    if FIRST_REPORT:
        FIRST_REPORT = False
        print(f"INFO: Using N chromosomes: {len(chroms)}", file=sys.stderr)

    if len(chroms) == 0:
        raise ValueError(
            f"No chromosomes left after filtering. Check your keep_chromosomes parameter in config.yaml. "
            f"Your fai file contains the following chromosomes: {FAI_DF['chr']}"
        )
    return chroms


def get_manifest():
    manifest = config.get("manifest")
    if manifest is None:
        raise ValueError("manifest parameter is missing in config.yaml")
    if not os.path.exists(manifest):
        raise ValueError(f"Manifest file {manifest} does not exist")
    manifest = pd.read_csv(config["manifest"], sep=r"\s+", comment="#").set_index(
        "sample"
    )
    return manifest


def get_input_bam(wc):
    return MANIFEST.loc[wc.sm, "bam"]


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


def get_hap_col_suffix(wc):
    if wc.hp == "all":
        return ""
    elif wc.hp == "hap1":
        return "_H1"
    elif wc.hp == "hap2":
        return "_H2"
    else:
        raise ValueError(f"Unknown haplotype {wc.hp}")


def pileup_cut_cmd(wc):
    if wc.hp == "all":
        tail = ""
    elif wc.hp == "hap1":
        tail = "_H1"
    elif wc.hp == "hap2":
        tail = "_H2"
    else:
        raise ValueError(f"Unknown haplotype {wc.hp}")
    if wc.el_type == "nucleosome":
        col = f"$nuc_coverage{tail}"
    elif wc.el_type == "linker":
        col = f"$msp_coverage{tail}-$fire_coverage{tail}"
    elif wc.el_type == "fire":
        col = f"$fire_coverage{tail}"
    else:
        raise ValueError(f"Unknown element type {wc.el_type}")
    return f"bioawk -tc hdr '{{print $1,$2,$3,{col}}}'"
