#!/usr/bin/env python
import defopt
import sys
import logging
from pathlib import Path
from typing import Optional
import pandas as pd
import polars as pl
import numpy as np
from numba import njit

FIRE_COLUMNS = [
    "chrom",
    "start",
    "end",
    "fiber",
    "score",
    "strand",
    "thick_start",
    "thick_end",
    "item_rgb",
    "fdr",
    "hap",
]


@njit
def fire_scores_per_chrom(
    starts, ends, q_values, chrom_length, coverage_array, max_add=20
):
    fire_scores = np.zeros(chrom_length, dtype=np.float64)
    q_values_t = -10.0 * np.log10(q_values)
    for start, end, q in zip(starts, ends, q_values_t):
        if end >= chrom_length:
            continue
        fire_scores[start:end] += min(q, max_add)
    # correct for regions with high coverage
    fire_scores = fire_scores / (coverage_array + 1)
    return fire_scores


@njit
def fdr_from_fire_scores(fire_scores, null_fire_scores, thresholds):
    Rs = np.zeros(thresholds.shape[0], dtype=np.int64)
    Vs = np.zeros(thresholds.shape[0], dtype=np.int64)
    for score in fire_scores:
        Rs += score > thresholds
    for score in null_fire_scores:
        Vs += score > thresholds
    fdrs = Vs / Rs
    return (fdrs, Vs, Rs)


@njit
def get_coverage_array(starts, ends, coverages, chrom_length):
    coverage_array = np.zeros(chrom_length, dtype=np.float64)
    for start, end, cov in zip(starts, ends, coverages):
        coverage_array[start:end] = cov
    return coverage_array


@njit
def make_coverage_array(starts, ends, chrom_length):
    coverage_array = np.zeros(chrom_length, dtype=np.float64)
    for start, end in zip(starts, ends):
        coverage_array[start:end] += 1
    return coverage_array


def fire_tracks(fire):
    null_s = []
    fire_s = []
    for chrom, g in fire.groupby("chrom", maintain_order=True):
        # fibers for this chromosome
        fibers = (
            g[
                [
                    "chrom",
                    "fiber",
                    "fiber_start",
                    "fiber_end",
                    "null_fiber_start",
                    "null_fiber_end",
                ]
            ]
            .unique()
            .to_pandas()
        )
        # convert to pandas for easier manipulation
        g = g.to_pandas()

        # get coverage for this chromosome and the shuffled fibers
        chrom_length = g.length[0]
        coverage_array = make_coverage_array(
            fibers.fiber_start.values, fibers.fiber_end.values, chrom_length
        )
        null_coverage_array = make_coverage_array(
            fibers.null_fiber_start.values, fibers.null_fiber_end.values, chrom_length
        )

        # find offset to use based on the shuffled fiber
        g["offset"] = g.null_fiber_start - g.fiber_start
        g["null_start"] = g.start + g.offset
        g["null_end"] = g.end + g.offset

        #
        fire_scores = fire_scores_per_chrom(
            g.start.values,
            g.end.values,
            g.fdr.values,
            g.length.max(),
            coverage_array,
        )
        null_fire_scores = fire_scores_per_chrom(
            g.null_start.values,
            g.null_end.values,
            g.fdr.values,
            g.length.max(),
            null_coverage_array,
        )

        logging.info(
            f"{chrom}:{fire_scores.shape[0]:,}\t"
            f"Max real FIRE score: {fire_scores.max():,.8}\t"
            f"Max null FIRE score: {null_fire_scores.max():,.8}"
        )
        null_s.append(null_fire_scores)
        fire_s.append(fire_scores)

    # all data
    fire_scores = np.concatenate(fire_s)
    null_fire_scores = np.concatenate(null_s)

    # Calculate FDR thresholds
    qs = np.arange(0.95, 1.0, 0.001)
    thresholds = np.quantile(fire_scores, qs)
    FDRs, Vs, Rs = fdr_from_fire_scores(fire_scores, null_fire_scores, thresholds)
    results = pd.DataFrame(
        {"threshold": thresholds, "FDR": FDRs, "shuffled_peaks": Vs, "peaks": Rs}
    )
    logging.info(f"\n{results}")
    logging.info(
        f"all:{fire_scores.shape[0]:,}\t{fire_scores.max():,.8}\t{null_fire_scores.max():,.8}"
    )


def analysis(fire, fibers):
    logging.info("Starting analysis")
    fire = fire.join(fibers, on=["chrom", "fiber"])
    logging.debug(f"Joined fibers\n{fire}")
    fire_tracks(fire)
    return 0


def main(
    infile: Path,
    fiber_locations_file: Path,
    shuffled_locations_file: Path,
    genome_file: Path,
    outfile: Optional[Path] = None,
    *,
    n_rows: Optional[int] = None,
    verbose: int = 0,
):
    """
    Author Mitchell R. Vollger

    :param infile: Input file, stdin by default
    :param outfile: Output file, stdout by default
    :param count: Number of times to display the greeting
    :param verbose: Set the logging level of the function
    """
    if infile is None:
        infile = sys.stdin
    if outfile is None:
        outfile = sys.stdout

    logger = logging.getLogger()
    log_format = "[%(levelname)s][Time elapsed (ms) %(relativeCreated)d]: %(message)s"
    log_level = 10 * (3 - verbose)
    logging.basicConfig(format=log_format)
    logger.setLevel(log_level)

    fai = pl.read_csv(
        genome_file,
        separator="\t",
        has_header=False,
        columns=[0, 1],
        new_columns=["chrom", "length"],
    )
    fiber_locations = pl.read_csv(
        fiber_locations_file,
        separator="\t",
        has_header=False,
        columns=[0, 1, 2, 3],
        new_columns=["chrom", "fiber_start", "fiber_end", "fiber"],
    )
    shuffled_locations = pl.read_csv(
        shuffled_locations_file,
        separator="\t",
        has_header=False,
        columns=[0, 1, 2, 3],
        new_columns=["chrom", "null_fiber_start", "null_fiber_end", "fiber"],
    )
    fibers = fiber_locations.join(shuffled_locations, on=["chrom", "fiber"]).join(
        fai, on="chrom"
    )
    fire = pl.read_csv(
        infile,
        separator="\t",
        has_header=False,
        new_columns=FIRE_COLUMNS,
        comment_char="#",
        batch_size=100_000,
        n_rows=n_rows,
    )
    analysis(fire, fibers)
    return 0


if __name__ == "__main__":
    defopt.run(main, show_types=True, version="0.0.1")
