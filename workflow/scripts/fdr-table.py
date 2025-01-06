#!/usr/bin/env python
import defopt
import logging
from pathlib import Path
from typing import Optional
import pandas as pd
import polars as pl
import numpy as np
import polars.selectors as cs
import gzip
import sys

# from numba import njit
ROLLING_FIRE_SCORE_WINDOW_SIZE = 200


def is_gzipped(path):
    with open(path, "rb") as f:
        return f.read(2) == b"\x1f\x8b"


def find_nearest(array, value):
    idx = np.searchsorted(array, value, side="left")
    idx[idx < 0] = 0
    idx[idx >= len(array)] = len(array) - 1
    return idx


def my_read_csv(*args, **kwargs):
    try:
        result = pl.read_csv(*args, **kwargs)
    # do some transformation with the dataframe
    except pl.exceptions.NoDataError as e:
        print(
            "No data is found in the input file. Check the input file and make sure it is not empty. It is likely that the input data was not generated correctly or that it was impossible to find peaks at the specified FDR value.",
            file=sys.stderr,
        )
        print(e, file=sys.stderr)
        sys.exit(1)
    return result


# ['#chrom', 'start', 'end', 'coverage', 'fire_coverage', 'score', 'nuc_coverage', 'msp_coverage',
# 'coverage_H1', 'fire_coverage_H1', 'score_H1', 'nuc_coverage_H1', 'msp_coverage_H1',
# 'coverage_H2', 'fire_coverage_H2', 'score_H2', 'nuc_coverage_H2', 'msp_coverage_H2']
def read_pileup_file(infile, nrows):
    # get the header from the first line of the file
    header = my_read_csv(infile, separator="\t", n_rows=1).columns

    # check that there is at least two lines
    open_infile = gzip.open if is_gzipped(infile) else open
    with open_infile(infile) as f:
        for i, _ in enumerate(f):
            if i > 1:
                break
        if i < 2:
            return None

    # add scema overrides for the score columns
    schema_overrides = {}
    for n in ["score", "score_H1", "score_H2", "score_shuffled"]:
        if n in header:
            schema_overrides[n] = float

    logging.info(f"Header of the pileup file:\n{header}")
    logging.info(f"Schema overrides for the pileup file:\n{schema_overrides}")

    # read the file
    pileup = my_read_csv(
        infile,
        separator="\t",
        has_header=False,
        new_columns=header,
        comment_prefix="#",
        n_rows=nrows,
        infer_schema_length=100000,
        schema_overrides=schema_overrides,
    )
    logging.info(f"Done reading pileup file:\n{pileup}")
    return pileup


# @njit
def fdr_from_fire_scores(fire_scores):
    Vs = []
    Rs = []
    Ts = []
    cur_R = 0.0
    cur_V = 0.0
    pre_score = -1.0
    first = True
    for score, is_real, bp in fire_scores.iter_rows():
        # save the counts and thresholds as long as we have counts
        if score != pre_score and cur_R > 0 and not first:
            Rs.append(cur_R)
            Vs.append(cur_V)
            Ts.append(pre_score)
        # don't add negative scores to the fdr data, since they have no coverage.
        if score < 0.0:
            break
        # update the counts
        if is_real:
            cur_R += bp
        else:
            cur_V += bp
        # prepare for next iteration
        pre_score = score
        first = False

    # add the last threshold with an FDR of 1
    Rs.append(1)
    Vs.append(1)
    Ts.append(-1.0)

    # set up return values
    Vs = np.array(Vs)
    Rs = np.array(Rs)
    Ts = np.array(Ts)
    FDRs = Vs / Rs
    FDRs[FDRs > 1] = 1.0

    return (Ts, FDRs, Vs, Rs)


def fdr_table_from_scores(fire_scores):
    # Calculate FDR thresholds
    Ts, FDRs, Vs, Rs = fdr_from_fire_scores(fire_scores)
    results = pd.DataFrame(
        {
            "threshold": Ts,
            "FDR": FDRs,
            "shuffled_bp": Vs,
            "real_bp": Rs,
        }
    )
    # simplify the results a little, don't want 100,000s of thresholds
    results = results.groupby("FDR", sort=False).tail(1).reset_index(drop=True)
    results = results.groupby("shuffled_bp", sort=False).tail(1).reset_index(drop=True)
    results = results.groupby("real_bp", sort=False).tail(1).reset_index(drop=True)
    # limit the number of thresholds that can be in the table
    results["threshold"] = results["threshold"].round(2)
    results = results.groupby("threshold", sort=False).tail(1).reset_index(drop=True)
    # sort the results by threshold so that they are now acceding
    # which is needed for the find_nearest function
    results = results.sort_values("threshold")
    logging.info(f"FDR results\n{results}")
    return results


def make_fdr_table(infile, outfile, nrows, max_cov=None, min_cov=None, max_fdr=0.05):
    # read the pileup file
    pileup = read_pileup_file(infile, nrows)
    # filter on coverages if needed
    if max_cov is not None:
        pileup = pileup.filter(
            pl.col("coverage") <= max_cov, pl.col("coverage_shuffled") <= max_cov
        )
    if min_cov is not None:
        pileup = pileup.filter(
            pl.col("coverage") >= min_cov, pl.col("coverage_shuffled") >= min_cov
        )

    # aggregate by the score and weight the score by the number of bases
    fire_scores = (
        pileup.melt(
            value_vars=["score", "score_shuffled"],
            id_vars=["#chrom", "start", "end"],
            variable_name="type",
            value_name="score",
        )
        .with_columns(
            pl.when(pl.col("type") == "score")
            .then(True)
            .otherwise(False)
            .alias("is_real"),
            pl.col("end").sub(pl.col("start")).alias("bp"),
        )
        .group_by(["score", "is_real"])
        .agg(pl.sum("bp").alias("bp"))
        .sort("score", descending=True)
    )

    # count bases in each category
    sums = fire_scores.group_by("is_real").agg(pl.sum("bp").alias("Mbp") / 1_000_000)
    logging.info(f"Number of Mbp in each category:\n{sums}")

    logging.info(f"Done aggregating pileup file:\n{fire_scores}")
    fdr_table = fdr_table_from_scores(fire_scores)
    fdr_table.to_csv(outfile, sep="\t", index=False)
    # raise an error if no threshold below 0.05 is found
    if fdr_table["FDR"].min() > max_fdr:
        raise ValueError(
            f"No FIRE score threshold has an FDR < {max_fdr}. Check the input Fiber-seq data with the QC pipeline and make sure you are using WGS Fiber-seq data."
        )
    return fdr_table


def read_fdr_table(infile):
    fdr_table = my_read_csv(infile, separator="\t").to_pandas()
    logging.info(f"Read FDR table:\n{fdr_table}")
    return fdr_table


def apply_fdr_table(infile, outfile, fdr_table, nrows):
    pileup = read_pileup_file(infile, nrows)
    # there is no input data
    if pileup is None:
        Path(outfile).touch()
        return

    logging.info(f"Applying FDR table to pileup file:\n{pileup}")
    # add a new column that reports the largest score in a centered window of with ROLLING_FIRE_SCORE_WINDOW_SIZE number of bases
    rolling_max_score = (
        pileup.rolling(
            index_column="start",
            period=f"{ROLLING_FIRE_SCORE_WINDOW_SIZE}i",
            offset=f"-{ROLLING_FIRE_SCORE_WINDOW_SIZE // 2}i",
            closed="both",
            group_by="#chrom",
        )
        .agg(
            pl.max("score").alias("max_window_score").fill_null(-1.0),
        )
        .drop("#chrom")
    )

    # add the max window score to the pileup
    pileup = (
        pileup.with_columns(rolling_max_score)
        .with_columns(
            pl.when(
                (pl.col("score") == pl.col("max_window_score"))
                & (pl.col("score") > 0.0)
            )
            .then(True)
            .otherwise(False)
            .alias("is_local_max"),
        )
        .with_columns(
            pl.col("is_local_max").rle_id().alias("local_max_group"),
        )
    )

    # group by local_max_group and find the midpoint (bp) of the group
    middle = (
        pileup.group_by(["#chrom", "local_max_group"])
        .agg(
            pl.col("end").sub(pl.col("start")).sum().alias("width"),
            pl.min("start").alias("first_start"),
        )
        .with_columns(
            # find the middle of the group
            pl.col("first_start").add(pl.col("width") // 2).alias("middle"),
        )
        .drop("width", "first_start")
    )

    # find the middle row of local max groups
    pileup = (
        pileup.join(middle, on=["#chrom", "local_max_group"], how="left")
        .with_columns(
            (
                (pl.col("start") <= pl.col("middle"))
                & (pl.col("end") > pl.col("middle"))
            ).alias("is_middle"),
        )
        .with_columns(
            # only keep the local maxes that are in the middle of the group
            (pl.col("is_local_max") & pl.col("is_middle")).alias("is_local_max"),
        )
        .drop("middle", "is_middle", "local_max_group", "max_window_score")
        .drop(cs.ends_with("_shuffled"))
    )

    # find the FDRs for the thresholds
    fdr_idx = find_nearest(fdr_table.threshold.values, pileup["score"].to_numpy())
    FDRs = fdr_table.FDR.values[fdr_idx]
    pileup = (
        pileup.with_columns(
            FDR=FDRs,
        )
        .with_columns(
            pl.when(pl.col("score") >= 0.0)
            .then(pl.col("FDR"))
            .otherwise(1.0)
            .alias("FDR"),
        )
        .with_columns(
            (-10 * np.log10(pl.col("FDR"))).alias("log_FDR").replace(float("inf"), 100),
        )
    )

    logging.info(f"Done calculating max window score:\n{pileup}")
    # write the pileup to a file
    pileup.write_csv(outfile, separator="\t")


def main(
    infile: Path,
    outfile: Path,
    *,
    fdr_table: Path = None,
    nrows: Optional[int] = None,
    max_cov: Optional[int] = None,
    min_cov: Optional[int] = None,
    max_fdr: float = 0.05,
    verbose: int = 0,
):
    """
    Author Mitchell R. Vollger

    :param infile: ft pileup track with shuffled fiber data
    :param outfile: FIRE score to FDR table
    :param verbose: Set the logging level of the function
    """
    logger = logging.getLogger()
    log_format = "[%(levelname)s][Time elapsed (ms) %(relativeCreated)d]: %(message)s"
    log_level = 10 * (3 - verbose)
    logging.basicConfig(format=log_format)
    logger.setLevel(log_level)

    if fdr_table is not None:
        fdr_table = read_fdr_table(fdr_table)
        apply_fdr_table(infile, outfile, fdr_table, nrows)
    else:
        fdr_table = make_fdr_table(
            infile, outfile, nrows, min_cov=min_cov, max_cov=max_cov, max_fdr=max_fdr
        )
    return 0


if __name__ == "__main__":
    defopt.run(main, show_types=True, version="0.0.1")
