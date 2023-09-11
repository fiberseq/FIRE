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


def is_grouped_with_previous(list_of_lists, min_frac_overlap=0.5):
    condition = []
    pre_list = set([])
    for cur_list in list_of_lists:
        cur_list = set(cur_list)
        overlap = len(cur_list.intersection(pre_list))
        overlap_frac = overlap / max(len(cur_list), len(pre_list))
        if overlap_frac >= min_frac_overlap:
            condition.append(True)
        else:
            condition.append(False)
        pre_list = cur_list
    return condition


def main(
    *,
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
    df = pl.read_csv(inf, separator="\t").with_columns(
        FIRE_IDs=pl.col("FIRE_IDs").str.split(",").cast(pl.List(pl.UInt32)),
    )

    df = (
        df.with_columns(
            pl.Series(
                name="shares_FIREs", values=is_grouped_with_previous(df["FIRE_IDs"])
            ),
        )
        .with_columns(
            (~pl.col("shares_FIREs")).cumsum().alias("group"),
        )
        .sort(["#chrom", "start"])
    )
    logging.info(f"\n{df.columns}")
    logging.info(f"\n{df}")
    df = (
        df.with_columns(
            pl.col("score").max().over("group").suffix("_max"),
            start=pl.col("peak_start").min().over("group").cast(pl.UInt32),
            end=pl.col("peak_end").max().over("group").cast(pl.UInt32),
            local_max_count=pl.col("group").len().over("group"),
        )
        .filter(pl.col("score") == pl.col("score_max"))
        .with_columns(
            peak_length=pl.col("end") - pl.col("start"),
        )
    ).drop("score_max", "group", "FIRE_IDs", "shares_FIREs", "is_local_max")
    logging.info(f"\n{df}")
    df.write_csv("/dev/stdout", separator="\t")
    return 0


if __name__ == "__main__":
    defopt.run(main, show_types=True, version="0.0.1")
