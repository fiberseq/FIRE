import re
import logging
import sys
import warnings
from functools import lru_cache

import pysam

# marks a manifest cell that falls back to the config value
MANIFEST_NA = "."


@lru_cache(maxsize=None)
def _bam_contigs(bam):
    """Read (name, length) pairs from the header of a BAM/CRAM, once per file."""
    verbosity = pysam.set_verbosity(0)
    try:
        with pysam.AlignmentFile(bam, check_sq=False, require_index=False) as f:
            contigs = tuple(zip(f.references, f.lengths))
    except (OSError, ValueError) as e:
        raise ValueError(f"FIRE: cannot read input bam {bam}: {e}") from e
    finally:
        pysam.set_verbosity(verbosity)
    if not contigs:
        raise ValueError(
            f"FIRE: input bam {bam} has no reference sequences (@SQ) in its header; "
            "FIRE requires an aligned bam"
        )
    return contigs


@lru_cache(maxsize=None)
def _sample_chroms(sm):
    """Filtered chromosome names for one sample, in BAM header order.

    Header order is the sort order of the data (mosdepth output, coordinate
    sorted reads), so it must be preserved for the sorted bedtools
    operations downstream.
    """
    min_contig_length = config.get("min_contig_length", 0)
    try:
        contigs = _bam_contigs(MANIFEST.loc[sm, "bam"])
    except ValueError as e:
        raise ValueError(f"{e} (sample '{sm}')") from e
    skipped = [name for name, length in contigs if length < min_contig_length]
    if skipped:
        print(
            f"WARNING: {sm}: skipping contigs with length < {min_contig_length:,}: {skipped}",
            file=sys.stderr,
        )
    chroms = tuple(
        name
        for name, length in contigs
        if length >= min_contig_length
        and "chrUn_" not in name
        and "_random" not in name
        and re.fullmatch(KEEP_CHRS, name)
    )
    print(f"INFO: {sm}: using {len(chroms)} chromosomes", file=sys.stderr)
    if not chroms:
        raise ValueError(
            f"FIRE: no chromosomes left for sample '{sm}' after filtering. "
            "Check the keep_chromosomes and min_contig_length options in config.yaml. "
            f"The bam header contains: {[name for name, _ in contigs]}"
        )
    return chroms


def get_ref(wc):
    return MANIFEST.loc[wc.sm, "ref"]


def get_fai(wc):
    return f"{get_ref(wc)}.fai"


def get_ref_name(wc):
    return MANIFEST.loc[wc.sm, "ref_name"]


def get_chroms(wc):
    return list(_sample_chroms(wc.sm))


def all_chroms():
    # sorted only for a deterministic wildcard-constraint regex; the
    # alternation order has no effect on matching
    return sorted({chrom for sm in MANIFEST.index for chrom in _sample_chroms(sm)})


def get_excludes(wc):
    excludes = list(config.get("excludes", []))
    if get_ref_name(wc) in ["hg38", "GRCh38"]:
        files = [
            "../annotations/hg38.gap.bed.gz",
            "../annotations/hg38.blacklist.ENCFF356LFX.bed.gz",
            "../annotations/SDs.merged.hg38.bed.gz",
        ]
        excludes += [workflow.source_path(file) for file in files]
    return excludes


def _config_ref_value(col):
    """The config value for ref/ref_name, or None when absent or empty."""
    value = config.get(col)
    if value is None or str(value) == "":
        return None
    return str(value)


def _fill_manifest_refs(manifest):
    """Fill and validate the ref and ref_name manifest columns."""
    empty_in_config = [
        col
        for col in ["ref", "ref_name"]
        if col in config and _config_ref_value(col) is None
    ]
    if empty_in_config:
        raise ValueError(
            f"FIRE: config options {empty_in_config} are empty in config.yaml; "
            "set a value or remove the key"
        )
    in_manifest = [col for col in ["ref", "ref_name"] if col in manifest.columns]
    in_config = [
        col for col in ["ref", "ref_name"] if _config_ref_value(col) is not None
    ]
    if len(in_manifest) == 1:
        raise ValueError(
            "FIRE: manifest columns 'ref' and 'ref_name' must be provided together "
            f"(found only '{in_manifest[0]}')"
        )
    if len(in_config) == 1:
        raise ValueError(
            "FIRE: config options 'ref' and 'ref_name' must be provided together "
            f"(found only '{in_config[0]}')"
        )
    if not in_manifest and not in_config:
        raise ValueError(
            "FIRE: no reference specified: add 'ref' and 'ref_name' columns to the "
            "manifest, or set 'ref' and 'ref_name' in config.yaml"
        )
    if in_manifest and in_config:
        print(
            "INFO: manifest ref/ref_name columns override config-level values",
            file=sys.stderr,
        )
    for col in ["ref", "ref_name"]:
        if col not in manifest.columns:
            manifest[col] = _config_ref_value(col)
            continue
        sentinel = manifest[col] == MANIFEST_NA
        if sentinel.any():
            if _config_ref_value(col) is None:
                missing = manifest.index[sentinel].tolist()
                raise ValueError(
                    f"FIRE: samples {missing} use '{MANIFEST_NA}' for '{col}' in the "
                    f"manifest, but '{col}' is not set in config.yaml"
                )
            manifest.loc[sentinel, col] = _config_ref_value(col)
    return manifest


def get_manifest():
    manifest_path = config.get("manifest")
    if manifest_path is None:
        raise ValueError("FIRE: manifest parameter is missing in config.yaml")
    if not os.path.exists(manifest_path):
        raise ValueError(f"FIRE: manifest file {manifest_path} does not exist")
    try:
        # dtype=str + keep_default_na=False keep every cell as literal text:
        # numeric sample names stay strings and a sample named NA stays "NA"
        # (missing trailing cells still parse as NaN; the malformed-row check
        # below catches both NaN and ""). index_col=False stops pandas from
        # silently treating the first field as an index when every data row
        # has one extra column; promoting ParserWarning to an error turns
        # the resulting silent field drop into a loud failure
        with warnings.catch_warnings():
            warnings.simplefilter("error", pd.errors.ParserWarning)
            manifest = pd.read_csv(
                manifest_path,
                sep=r"\s+",
                comment="#",
                dtype=str,
                keep_default_na=False,
                index_col=False,
                engine="python",
            )
    except (
        pd.errors.ParserError,
        pd.errors.EmptyDataError,
        pd.errors.ParserWarning,
    ) as e:
        raise ValueError(
            f"FIRE: cannot parse manifest {manifest_path}: {e} (check that every "
            "row has the same number of whitespace-separated fields as the header)"
        ) from e
    for col in ["sample", "bam"]:
        if col not in manifest.columns:
            raise ValueError(
                f"FIRE: manifest must have 'sample' and 'bam' columns; "
                f"found: {list(manifest.columns)}"
            )
    if len(manifest) == 0:
        raise ValueError(f"FIRE: manifest {manifest_path} has no samples")
    dups = manifest["sample"][manifest["sample"].duplicated()].tolist()
    if dups:
        raise ValueError(f"FIRE: duplicate sample names in manifest: {dups}")
    for sm in manifest["sample"]:
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", sm):
            raise ValueError(f"FIRE: sample name '{sm}' must match [A-Za-z0-9_.-]+")
    manifest = manifest.set_index("sample")
    manifest = _fill_manifest_refs(manifest)
    ref_cols = manifest[["bam", "ref", "ref_name"]]
    malformed = manifest.index[(ref_cols.isna() | ref_cols.eq("")).any(axis=1)]
    if len(malformed) > 0:
        raise ValueError(
            f"FIRE: samples {malformed.tolist()} have missing or malformed manifest "
            "fields; every row must fill all manifest columns "
            f"(use '{MANIFEST_NA}' in ref/ref_name to fall back to the config value)"
        )
    manifest["ref"] = manifest["ref"].map(os.path.abspath)
    for ref in manifest["ref"].unique():
        if not os.path.isfile(ref):
            raise ValueError(f"FIRE: reference file {ref} does not exist")
        if not os.path.isfile(f"{ref}.fai"):
            raise ValueError(f"FIRE: reference index file {ref}.fai does not exist")
    return manifest


def get_input_bam(wc):
    return MANIFEST.loc[wc.sm, "bam"]


def genome_file_content(sm):
    """The genome (chrom sizes) file text for one sample: all bam header
    contigs, in header order."""
    return "".join(
        f"{name}\t{length}\n" for name, length in _bam_contigs(MANIFEST.loc[sm, "bam"])
    )


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
