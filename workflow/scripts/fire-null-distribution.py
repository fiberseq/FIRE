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
FIBER_COLUMNS = [
    "chrom",
    "fiber",
    "fiber_start",
    "fiber_end",
    "null_fiber_start",
    "null_fiber_end",
]


def rle(inarray):
    """run length encoding. Partial credit to R rle function.
    Multi datatype arrays catered for including non Numpy
    returns: tuple (runlengths, startpositions, values)"""
    ia = np.asarray(inarray)  # force numpy
    if ia.size == 0:
        return (ia, ia, ia)
    else:
        n = ia.shape[0]
        y = ia[1:] != ia[:-1]  # pairwise unequal (string safe)
        i = np.append(np.where(y), n - 1)  # must include last element posi
        z = np.diff(np.append(-1, i))  # run lengths
        p = np.cumsum(np.append(0, z))[:-1]  # positions
        return (z, p, ia[i])


def bed_rle(inarray):
    run_lengths, starts, scores = rle(inarray)
    starts = starts.astype(int)
    ends = starts + run_lengths.astype(int)
    return np.array([starts, ends, scores]).transpose()


@njit
def fire_scores_per_chrom(
    starts, ends, q_values, chrom_length, coverage_array, min_allowed_q=0.01
):
    fire_scores = np.zeros(chrom_length, dtype=np.float64)

    multi = -50.0  # a multi of -50 and a min_allowed_q of 0.01 gives a max score of 100
    max_add = multi * np.log10(min_allowed_q)
    q_values_t = multi * np.log10(q_values)
    for start, end, q in zip(starts, ends, q_values_t):
        if end >= chrom_length:
            continue
        fire_scores[start:end] += min(q, max_add)

    # allow division with coverage to always work
    coverage_array[coverage_array == 0] = 1
    # correct for coverage
    fire_scores = fire_scores / coverage_array
    return fire_scores


@njit
def fdr_from_fire_scores(fire_scores):
    Vs = []
    Rs = []
    Ts = []
    cur_R = 0.0
    cur_V = 0.0
    pre_score = 0.0
    for start, end, score, is_real in fire_scores:
        # save the counts and thresholds as long as we have counts
        if score != pre_score and cur_R > 0:
            Rs.append(cur_R)
            Vs.append(cur_V)
            Ts.append(score)
        # update the counts
        counts = end - start
        if is_real:
            cur_R += counts
        else:
            cur_V += counts
        # prepare for next iteration
        pre_score = score
    # set up return values
    Vs = np.array(Vs)
    Rs = np.array(Rs)
    Ts = np.array(Ts)
    FDRs = Vs / Rs
    FDRs[FDRs > 1] = 1.0
    return (Ts, FDRs, Vs, Rs)


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


def fire_tracks(fire, outfile):
    null_s = []
    fire_s = []
    for chrom, g in fire.group_by("chrom", maintain_order=True):
        logging.info(f"Processing {chrom}")
        # fibers for this chromosome
        fibers = g[FIBER_COLUMNS].unique().to_pandas()
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
        expected_median_coverage = np.median(null_coverage_array)

        # find offset to use based on the shuffled fiber
        g["offset"] = g.null_fiber_start - g.fiber_start
        g["null_start"] = g.start + g.offset
        g["null_end"] = g.end + g.offset

        #
        rle_fire_scores = bed_rle(
            fire_scores_per_chrom(
                g.start.values,
                g.end.values,
                g.fdr.values,
                g.length.max(),
                coverage_array,
            )
        )
        rle_null_scores = bed_rle(
            fire_scores_per_chrom(
                g.null_start.values,
                g.null_end.values,
                g.fdr.values,
                g.length.max(),
                null_coverage_array,
            )
        )

        logging.info(
            f"{chrom}: {rle_fire_scores.shape[0]:,}\t"
            f"Max real FIRE score: {rle_fire_scores[:,2].max():,.8}\t"
            f"Max null FIRE score: {rle_null_scores[:,2].max():,.8}\t"
            f"Expected median coverage: {expected_median_coverage}"
        )
        fire_s.append(rle_fire_scores)
        null_s.append(rle_null_scores)

    # all data
    fire_scores = np.concatenate(fire_s)
    null_fire_scores = np.concatenate(null_s)
    logging.debug(f"rle fire score shape: {fire_scores.shape}")
    logging.info(
        f"all: {fire_scores.shape[0]:,}\t"
        f"Max real FIRE score: {fire_scores[:,2].max():,.8}\t"
        f"Max null FIRE score: {null_fire_scores[:,2].max():,.8}"
    )

    # convert to pandas for easier manipulation
    fire_scores = pd.DataFrame(fire_scores, columns=["start", "end", "score"])
    fire_scores["is_real"] = 1.0
    null_fire_scores = pd.DataFrame(null_fire_scores, columns=["start", "end", "score"])
    null_fire_scores["is_real"] = 0.0
    fire_scores = (
        pd.concat([fire_scores, null_fire_scores])
        .sort_values("score", ascending=False)
        .to_numpy()
    )
    logging.debug(f"Fire scores\n{fire_scores}")

    # Calculate FDR thresholds
    Ts, FDRs, Vs, Rs = fdr_from_fire_scores(fire_scores)
    results = pd.DataFrame(
        {
            "threshold": Ts,
            "FDR": FDRs,
            "shuffled_peaks": Vs,
            "peaks": Rs,
        }
    )
    # simplify the results a little, don't want 100,000s of thresholds
    results = results.groupby("FDR", sort=False).tail(1).reset_index(drop=True)
    results["threshold"] = results["threshold"].round(2)
    results = results.groupby("threshold", sort=False).tail(1).reset_index(drop=True)
    logging.info(f"FDR results\n{results}")
    results.to_csv(outfile, sep="\t", index=False)


def make_fdr_table(fire, fibers, outfile):
    logging.info("Starting analysis")
    fire = fire.join(fibers, on=["chrom", "fiber"])
    logging.debug(f"Joined fibers\n{fire}")
    fire_tracks(fire, outfile)
    return 0


def write_bed(chrom, rle_scores, out, first=True):
    df = pd.DataFrame(
        {
            "chrom": chrom,
            "st": rle_scores[:, 0],
            "en": rle_scores[:, 1],
            "score": rle_scores[:, 2],
        }
    )
    if first:
        mode = "w"
    else:
        mode = "a"
    df.to_csv(out, mode=mode, header=False, index=False, sep="\t")


def write_scores(fire, fibers, outfile):
    fire = fire.join(fibers, on=["chrom", "fiber"])
    first = True
    for chrom, g in fire.group_by("chrom", maintain_order=True):
        logging.info(f"Processing {chrom}")
        # fibers for this chromosome
        fibers = g[FIBER_COLUMNS].unique().to_pandas()
        # convert to pandas for easier manipulation
        g = g.to_pandas()

        # get coverage for this chromosome and the shuffled fibers
        chrom_length = g.length[0]
        coverage_array = make_coverage_array(
            fibers.fiber_start.values, fibers.fiber_end.values, chrom_length
        )
        #
        rle_fire_scores = bed_rle(
            fire_scores_per_chrom(
                g.start.values,
                g.end.values,
                g.fdr.values,
                g.length.max(),
                coverage_array,
            )
        )
        write_bed(chrom, rle_fire_scores, outfile, first=first)
        first = False


def main(
    infile: Path,
    fiber_locations_file: Path,
    genome_file: Path,
    outfile: Path,
    *,
    shuffled_locations_file: Optional[Path] = None,
    fdr_table_file: Optional[Path] = None,
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

    logging.info(f"Reading FIRE file: {infile}")
    fire = pl.read_csv(
        infile,
        separator="\t",
        has_header=False,
        new_columns=FIRE_COLUMNS,
        comment_char="#",
        batch_size=100_000,
        n_rows=n_rows,
    )
    logging.info(f"Reading genome file: {genome_file}")
    fai = pl.read_csv(
        genome_file,
        separator="\t",
        has_header=False,
        columns=[0, 1],
        new_columns=["chrom", "length"],
    )
    logging.info(f"Reading fiber locations file: {fiber_locations_file}")
    fiber_locations = pl.read_csv(
        fiber_locations_file,
        separator="\t",
        has_header=False,
        columns=[0, 1, 2, 3],
        new_columns=["chrom", "fiber_start", "fiber_end", "fiber"],
    ).join(fai, on="chrom")

    if shuffled_locations_file is not None:
        logging.info(
            f"Reading shuffled fiber locations file: {shuffled_locations_file}"
        )
        shuffled_locations = pl.read_csv(
            shuffled_locations_file,
            separator="\t",
            has_header=False,
            columns=[0, 1, 2, 3],
            new_columns=["chrom", "null_fiber_start", "null_fiber_end", "fiber"],
        )
        fibers = fiber_locations.join(shuffled_locations, on=["chrom", "fiber"])
        make_fdr_table(fire, fibers, outfile)
    else:
        fdr_table = pl.read_csv(fdr_table_file, sep="\t")
        write_scores(fire, fiber_locations, outfile)
    return 0


if __name__ == "__main__":
    defopt.run(main, show_types=True, version="0.0.1")
