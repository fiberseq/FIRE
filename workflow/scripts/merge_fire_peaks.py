#!/usr/bin/env python
import os
import defopt
import logging
from pathlib import Path
import numpy as np
from typing import Optional
import polars as pl
import io
import sys
from numba import njit


def is_grouped_with_previous(
    list_of_lists,
    starts,
    ends,
    min_frac_overlap=0.5,
    min_reciprocal_overlap=0.75,
):
    condition = []
    pre_st = 0
    pre_en = 1
    pre_list = set([])
    for cur_list, cur_st, cur_en in zip(list_of_lists, starts, ends):
        cur_list = set(cur_list)
        overlap = len(cur_list.intersection(pre_list))
        overlap_frac = overlap / max(len(cur_list), len(pre_list))
        overlap_bp = pre_en - cur_st
        reciprocal_overlap = min(
            overlap_bp / (pre_en - pre_st), overlap_bp / (cur_en - cur_st)
        )
        if (
            overlap_frac >= min_frac_overlap
            and reciprocal_overlap >= min_reciprocal_overlap
        ):
            condition.append(True)
        else:
            condition.append(False)
        pre_list = cur_list
        pre_st = cur_st
        pre_en = cur_en
    return condition


def group_peaks(df, min_frac_overlap=0.5, min_reciprocal_overlap=0.75):
    df = (
        df.sort(["#chrom", "peak_start"])
        .with_columns(
            pl.Series(
                name="shares_FIREs",
                values=is_grouped_with_previous(
                    df["FIRE_IDs"],
                    df["peak_start"],
                    df["peak_end"],
                    min_frac_overlap=min_frac_overlap,
                    min_reciprocal_overlap=min_reciprocal_overlap,
                ),
            ),
        )
        .with_columns(
            (~pl.col("shares_FIREs")).cumsum().alias("group"),
        )
        .with_columns(
            pl.col("score").max().over("group").suffix("_max"),
            peak_start=pl.col("peak_start").min().over("group").cast(pl.UInt32),
            peak_end=pl.col("peak_end").max().over("group").cast(pl.UInt32),
            local_max_count=pl.col("group").len().over("group"),
        )
        .filter(pl.col("score") == pl.col("score_max"))
        .with_columns(
            peak_length=pl.col("peak_end") - pl.col("peak_start"),
        )
    )
    return df


def main(
    *,
    max_score_every: int = 200,
    min_frac_overlap: float = 0.5,
    min_reciprocal_overlap: float = 0.0,
    max_grouping_iterations: int = 10,
    verbose: int = 0,
):
    """
    Author Mitchell R. Vollger
    """
    logger = logging.getLogger()
    log_format = "[%(levelname)s][Time elapsed (ms) %(relativeCreated)d]: %(message)s"
    log_level = 10 * (3 - verbose)
    logging.basicConfig(format=log_format)
    logger.setLevel(log_level)

    inf = io.StringIO(sys.stdin.read())
    df = (
        pl.read_csv(inf, separator="\t", null_values=".")
        # group into sliding X bp windows and only keep the highest score
        .with_columns(roll_start=pl.col("start"))
        .sort(["#chrom", "roll_start"])
        .groupby_rolling("roll_start", period=f"{max_score_every}i", by="#chrom")
        .agg([pl.exclude("roll_start").sort_by("score").last()])
        .drop("roll_start")
        # remove any peaks that are the highest score for multiple X bp windows
        .unique(subset=["#chrom", "peak_start", "peak_end", "start", "end"])
        .sort(["#chrom", "peak_start"])
        # convert the FIRE IDs strings to lists of ints
        .with_columns(
            FIRE_IDs=pl.col("FIRE_IDs").str.split(",").cast(pl.List(pl.UInt32)),
        )
    )
    logging.info(
        f"Dynamic window merging over {max_score_every} bp is done: {df.shape[0]:,}"
    )
    # group data 5 times
    n_row = None
    i = 0
    while i < max_grouping_iterations:
        df = group_peaks(
            df,
            min_frac_overlap=min_frac_overlap,
            min_reciprocal_overlap=min_reciprocal_overlap,
        )
        if n_row == df.shape[0]:
            break
        n_row = df.shape[0]
        i += 1
        logging.info(f"Iterative merging round {i} complete. Rows:{n_row:,}")

    (
        df.sort(["#chrom", "peak_start"])
        .drop("score_max", "group", "FIRE_IDs", "shares_FIREs", "is_local_max")
        .write_csv("/dev/stdout", separator="\t")
    )
    return 0


if __name__ == "__main__":
    defopt.run(main, show_types=True, version="0.0.1")
