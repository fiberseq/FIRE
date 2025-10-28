import re
import logging
import sys
import pysam

FIRST_REPORT = True


def get_ref_orig():
    if "ref" not in config:
        raise ValueError("FIRE: ref parameter is missing in config.yaml")
    ref = config["ref"]
    if not os.path.exists(ref):
        raise ValueError(f"FIRE: reference file {ref} does not exist")
    return os.path.abspath(ref)

def get_ref_old(wc):
    if "ref" in config:
        ref = config["ref"]
    if "ref" not in config:
        if "ref" in MANIFEST.columns:
            ref = get_input_ref()
        else: 
            raise ValueError("FIRE: ref parameter is missing in config.yaml and no ref column in manifest")
    #ref = config["ref"]
    if not os.path.exists(ref):
        raise ValueError(f"FIRE: reference file {ref} does not exist")
    return os.path.abspath(ref)

def get_ref(wc):
    ref = MANIFEST.loc[wc.sm, "ref"]
    if not os.path.exists(ref):
        raise ValueError(f"FIRE: reference file {ref} does not exist")
    return os.path.abspath(ref)

def get_ref_name_old(wc):
    if "ref_name" in config:
        ref_name=config["ref_name"]
    else:
        if "ref_name" in MANIFEST.columns:
            ref_name= get_input_ref_name()
        else:
            raise ValueError("FIRE: ref_name parameter is missing in config.yaml and no ref_name column in manifest")
    #ref=get_ref()
    #temp_split=ref.strip().split('/')[-1]
    #fa_split=temp_split.split('.fa')
    #ref_name = '.fa'.join(fa_split[:-1])
    return(ref_name)

def get_ref_name(wc):
    return MANIFEST.loc[wc.sm, "ref_name"]

def get_fai_orig(wc):
    fai = f"{get_ref(wc)}.fai"
    if not os.path.exists(fai):
        raise ValueError(f"FIRE: reference index file {fai} does not exist")
    return fai

def get_fai(wc):
    ref = MANIFEST.loc[wc.sm, "ref"]
    #fai= ref + ".fai"
    fai = str(ref) + ".fai"
    if not os.path.exists(fai):
        raise ValueError(f"FIRE: reference index file {fai} does not exist")
    return fai

def get_excludes(wc):
    excludes = config.get("excludes", [])
    ref_name = get_ref_name(wc)
    if ref_name == "hg38" or ref_name == "GRCh38":
        files = [
            "../annotations/hg38.gap.bed.gz",
            "../annotations/hg38.blacklist.ENCFF356LFX.bed.gz",
            "../annotations/SDs.merged.hg38.bed.gz",
        ]
        excludes += [workflow.source_path(file) for file in files]
    return excludes


def get_fai_df(wc):
    fai = get_fai(wc)
    return pd.read_csv(fai, sep="\t", names=["chr", "length", "x", "y", "z"])


def get_chroms_orig(wc):
    global FIRST_REPORT
    min_contig_length = config.get("min_contig_length", 0)
    fai_df=get_fai_df(wc)
    skipped_contigs = fai_df["chr"][fai_df["length"] < min_contig_length]
    if len(skipped_contigs) > 0 and FIRST_REPORT:
        print(
            f"WARNING: Skipping contigs with length < {min_contig_length:,}: {skipped_contigs}",
            file=sys.stderr,
        )

    chroms = fai_df["chr"][fai_df["length"] >= min_contig_length]
    chroms = sorted([chrom for chrom in chroms if "chrUn_" not in chrom])
    chroms = [chrom for chrom in chroms if "_random" not in chrom]
    chroms = [chrom for chrom in chroms if re.fullmatch(KEEP_CHRS, chrom)]

    if FIRST_REPORT:
        FIRST_REPORT = False
        print(f"INFO: Using N chromosomes: {len(chroms)}", file=sys.stderr)

    if len(chroms) == 0:
        raise ValueError(
            f"No chromosomes left after filtering. Check your keep_chromosomes parameter in config.yaml. "
            f"Your fai file contains the following chromosomes: {fai_df['chr']}"
        )
    return chroms

def get_chroms(wc):
    global FIRST_REPORT
    min_contig_length = config.get("min_contig_length", 0)
    fai_df=get_fai_df(wc)
    skipped_contigs = fai_df["chr"][fai_df["length"] < min_contig_length]
    if len(skipped_contigs) > 0 and FIRST_REPORT:
        print(
            f"WARNING: Skipping contigs with length < {min_contig_length:,}: {skipped_contigs}",
            file=sys.stderr,
        )

    bam_chr_list=[]
    input_bam_path=get_input_bam(wc)
    input_bam = pysam.AlignmentFile(input_bam_path, "rc", threads=MAX_THREADS)
    bam_header_dict = input_bam.header.to_dict()

    for line in bam_header_dict['SQ']:
        chr_name=line['SN']
        bam_chr_list.append(chr_name)

    input_bam.close()

    chroms = fai_df["chr"][fai_df["length"] >= min_contig_length]
    chroms = sorted([chrom for chrom in chroms if "chrUn_" not in chrom])
    chroms = [chrom for chrom in chroms if "_random" not in chrom]
    chroms = [chrom for chrom in chroms if re.fullmatch(KEEP_CHRS, chrom)]
    chroms = [chrom for chrom in chroms if chrom in bam_chr_list]

    if FIRST_REPORT:
        FIRST_REPORT = False
        print(f"INFO: Using N chromosomes: {len(chroms)}", file=sys.stderr)

    if len(chroms) == 0:
        raise ValueError(
            f"No chromosomes left after filtering. Check your keep_chromosomes parameter in config.yaml. "
            f"Your fai file contains the following chromosomes: {fai_df['chr']}"
        )
    return chroms

def get_chroms_first_element(wc):
    global FIRST_REPORT
    min_contig_length = config.get("min_contig_length", 0)
    fai_df=get_fai_df(wc)
    skipped_contigs = fai_df["chr"][fai_df["length"] < min_contig_length]
    if len(skipped_contigs) > 0 and FIRST_REPORT:
        print(
            f"WARNING: Skipping contigs with length < {min_contig_length:,}: {skipped_contigs}",
            file=sys.stderr,
        )

    bam_chr_list=[]
    input_bam_path=get_input_bam(wc)
    input_bam = pysam.AlignmentFile(input_bam_path, "rc", threads=MAX_THREADS)
    bam_header_dict = input_bam.header.to_dict()

    for line in bam_header_dict['SQ']:
        chr_name=line['SN']
        bam_chr_list.append(chr_name)

    input_bam.close()

    chroms = fai_df["chr"][fai_df["length"] >= min_contig_length]
    chroms = sorted([chrom for chrom in chroms if "chrUn_" not in chrom])
    chroms = [chrom for chrom in chroms if "_random" not in chrom]
    chroms = [chrom for chrom in chroms if re.fullmatch(KEEP_CHRS, chrom)]
    chroms = [chrom for chrom in chroms if chrom in bam_chr_list]

    if FIRST_REPORT:
        FIRST_REPORT = False
        print(f"INFO: Using N chromosomes: {len(chroms)}", file=sys.stderr)

    if len(chroms) == 0:
        raise ValueError(
            f"No chromosomes left after filtering. Check your keep_chromosomes parameter in config.yaml. "
            f"Your fai file contains the following chromosomes: {fai_df['chr']}"
        )
    return chroms[0]

def get_chroms_first_element_orig(wc):
    global FIRST_REPORT
    min_contig_length = config.get("min_contig_length", 0)
    fai_df=get_fai_df(wc)
    skipped_contigs = fai_df["chr"][fai_df["length"] < min_contig_length]
    if len(skipped_contigs) > 0 and FIRST_REPORT:
        print(
            f"WARNING: Skipping contigs with length < {min_contig_length:,}: {skipped_contigs}",
            file=sys.stderr,
        )

    chroms = fai_df["chr"][fai_df["length"] >= min_contig_length]
    chroms = sorted([chrom for chrom in chroms if "chrUn_" not in chrom])
    chroms = [chrom for chrom in chroms if "_random" not in chrom]
    chroms = [chrom for chrom in chroms if re.fullmatch(KEEP_CHRS, chrom)]

    if FIRST_REPORT:
        FIRST_REPORT = False
        print(f"INFO: Using N chromosomes: {len(chroms)}", file=sys.stderr)

    if len(chroms) == 0:
        raise ValueError(
            f"No chromosomes left after filtering. Check your keep_chromosomes parameter in config.yaml. "
            f"Your fai file contains the following chromosomes: {fai_df['chr']}"
        )
    return chroms[0]


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

def get_input_ref(wc):
    return MANIFEST.loc[wc.sm, "ref"]

def get_input_ref_name(wc):
    return MANIFEST.loc[wc.sm, "ref_name"]

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
