#!/usr/bin/env python
import defopt
import sys
import gc
import logging
from pathlib import Path
from typing import Optional
import pandas as pd
import polars as pl
import numpy as np
from numba import njit
import io

# from pyinstrument import Profiler


def chunker(seq, size):
    return (seq[pos : pos + size] for pos in range(0, len(seq), size))


def make_decorator(ct, fiber, score, strand, color, el_type, hp, st, en, starts, ends):
    start = starts[0]
    end = ends[-1]
    lengths = ",".join(map(str, ends - starts))
    offsets = ",".join(map(str, starts - start))
    block_count = len(starts)
    # chr1 12 9985 block 1000 + 12 9985 200,0,150
    # 382, 1 , 1
    # chr1:1-10000:LongRead
    # block 255,0,0,180 Ignored TypeA
    # 1 is transparent
    prime_color = "200,200,200,1"
    return (
        # bed9
        f"{ct}\t{start}\t{end}\t{el_type}\t{score}\t{strand}\t{start}\t{end}\t{prime_color}\t"
        # bed12
        f"{block_count}\t{lengths}\t{offsets}\t"
        # read tag for the decorator
        f"{ct}:{st}-{en}:{fiber}\t"
        # decorator
        f"block\t{color},0\tIgnored\t{el_type}"
    )


def subgroup(df, ct, fiber, strand, hp):
    st = df["st"].min()
    en = df["en"].max()
    # tmp = df.filter(pl.col("color") != "230,230,230")
    # linker = df.filter(pl.col("color") == "147,112,219")
    # fire = df.filter(pl.col("color") == "255,0,0")
    # for el_type, tdf in zip(["Linker", "FIRE"], [linker, fire]):
    for (color, score), gdf in df.group_by(["color", "score"]):
        if gdf.shape[0] == 0 or color == "230,230,230":
            continue
        elif color == "147,112,219":
            el_type = "Linker"
        else:
            el_type = "FIRE"
        decorator = make_decorator(
            ct,
            fiber,
            score,
            strand,
            color,
            el_type,
            hp,
            st,
            en,
            gdf["st"],
            gdf["en"],
        )
        print(decorator)

    return (
        ct,
        st,
        en,
        fiber,
        1,
        strand,
        0,
        0,
        "0,0,0,200",
        2,
        "1,1",
        f"0,{en-st-1}",
        hp,
    )


def process(df, outfile, group_size=5_000):
    data = []
    fibers = df["fiber"].unique()
    n_fibers = len(fibers)
    n = 0
    mode = "w"
    for (ct, fiber, hp), gdf in df.group_by(
        ["#ct", "fiber", "HP"], maintain_order=True
    ):
        strand = "."
        data.append(subgroup(gdf, ct, fiber, strand, hp))
        n += 1
        if n % group_size == 0 or n == n_fibers:
            logging.info(f"processed {n:,} fibers of {n_fibers:,}")
            pd.DataFrame(data).sort_values([0, 1, 2]).to_csv(
                outfile, sep="\t", header=False, index=False, mode=mode
            )
            mode = "a"
            data = []
            gc.collect()


def main(
    infile: str,
    outfile: Optional[Path],
    *,
    verbose: int = 0,
):
    """
    Author Mitchell R. Vollger
    :param infile: Input file, stdin by default
    :param outfile: Output file, stdout by default
    :param verbose: Set the logging level of the function
    """
    if infile == "-":
        infile = io.StringIO(sys.stdin.read())

    logger = logging.getLogger()
    log_format = "[%(levelname)s][Time elapsed (ms) %(relativeCreated)d]: %(message)s"
    log_level = 10 * (3 - verbose)
    logging.basicConfig(format=log_format)
    logger.setLevel(log_level)

    df = pl.read_csv(
        infile,
        separator="\t",
        low_memory=True,
        columns=[
            "#ct",
            "st",
            "en",
            "fiber",
            "score",
            # "strand",
            "HP",
            "color",
        ],
    )
    if df.shape[0] == 0:
        outfile = open(outfile, "w")
        return 0
    # df = df.filter(pl.col("color") != "230,230,230")

    # with Profiler(interval=0.1) as profiler:
    logging.info(f"{df}")
    process(df, outfile)
    # profiler.print()
    # profiler.open_in_browser()
    return 0


if __name__ == "__main__":
    defopt.run(main, show_types=True, version="0.0.1")
