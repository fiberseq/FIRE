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
from scipy.signal import argrelextrema

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
HAPS = ["H1", "H2"]


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
    starts,
    ends,
    q_values,
    chrom_length,
    coverage_array,
    min_allowed_q=0.01,
    min_coverage=4,
):
    fire_scores = np.zeros(chrom_length, dtype=np.float64)

    multi = -50.0  # a multi of -50 and a min_allowed_q of 0.01 gives a max score of 100
    max_add = multi * np.log10(min_allowed_q)
    q_values_t = multi * np.log10(q_values)
    for start, end, q in zip(starts, ends, q_values_t):
        if end >= chrom_length:
            continue
        fire_scores[start:end] += min(q, max_add)

    # correct for coverage
    fire_scores = fire_scores / coverage_array
    # correct divide by zeros
    fire_scores[np.isnan(fire_scores)] = 0.0
    # drop the scores that have no coverage
    fire_scores[coverage_array < min_coverage] = -1.0
    return fire_scores


@njit
def fdr_from_fire_scores(fire_scores):
    Vs = []
    Rs = []
    Ts = []
    cur_R = 0.0
    cur_V = 0.0
    pre_score = -1.0
    first = True
    for start, end, score, is_real in fire_scores:
        # save the counts and thresholds as long as we have counts
        if score != pre_score and cur_R > 0 and not first:
            Rs.append(cur_R)
            Vs.append(cur_V)
            Ts.append(pre_score)
        # don't add negative scores to the fdr data, since they have no coverage.
        if score < 0.0:
            break
        # update the counts
        counts = end - start
        if is_real:
            cur_R += counts
        else:
            cur_V += counts
        # prepare for next iteration
        pre_score = score
        first = False
    # set up return values
    Vs = np.array(Vs)
    Rs = np.array(Rs)
    Ts = np.array(Ts)
    FDRs = Vs / Rs
    FDRs[FDRs > 1] = 1.0
    return (Ts, FDRs, Vs, Rs)


@njit
def get_coverage_from_array(starts, ends, coverage_array, stat="median"):
    out_coverage = np.zeros(starts.shape[0], dtype=np.float64)
    idx = 0
    for start, end in zip(starts, ends):
        if stat == "median":
            val = np.median(coverage_array[start:end])
        elif stat == "max":
            val = np.max(coverage_array[start:end])
        else:
            val = np.mean(coverage_array[start:end])
        out_coverage[idx] = val
        idx += 1
    return out_coverage


@njit
def make_coverage_array(starts, ends, chrom_length):
    coverage_array = np.zeros(chrom_length, dtype=np.float64)
    for start, end in zip(starts, ends):
        coverage_array[start:end] += 1
    return coverage_array


def fire_tracks(fire, outfile, min_coverage=4):
    null_s = []
    fire_s = []
    logging.info(f"Fire data\n{fire}")
    for chrom, g in fire.groupby("chrom", maintain_order=True):
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
        expected_median_coverage = np.median(
            null_coverage_array[null_coverage_array > 0]
        )

        # find offset to use based on the shuffled fiber
        g["offset"] = g.null_fiber_start - g.fiber_start
        g["null_start"] = g.start + g.offset
        g["null_end"] = g.end + g.offset

        logging.info(
            f"real bp: {(g.end-g.start).sum():,}\t"
            f"null bp: {(g.null_end-g.null_start).sum():,}"
        )

        #
        rle_fire_scores = bed_rle(
            fire_scores_per_chrom(
                g.start.values,
                g.end.values,
                g.fdr.values,
                g.length.max(),
                coverage_array,
                min_coverage=min_coverage,
            )
        )
        rle_null_scores = bed_rle(
            fire_scores_per_chrom(
                g.null_start.values,
                g.null_end.values,
                g.fdr.values,
                g.length.max(),
                null_coverage_array,
                min_coverage=min_coverage,
            )
        )
        n_bp_considered = (coverage_array >= min_coverage).sum()
        logging.info(
            f"{chrom}: {n_bp_considered:,} of {chrom_length:,}\t"
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
    results = results[results.threshold > 0.0]
    results = results.groupby("FDR", sort=False).tail(1).reset_index(drop=True)
    results = (
        results.groupby("shuffled_peaks", sort=False).tail(1).reset_index(drop=True)
    )
    results = results.groupby("peaks", sort=False).tail(1).reset_index(drop=True)
    # limit the number of thresholds that can be in the table
    results["threshold"] = results["threshold"].round(2)
    results = results.groupby("threshold", sort=False).tail(1).reset_index(drop=True)
    logging.info(f"FDR results\n{results}")
    results.to_csv(outfile, sep="\t", index=False)


def make_fdr_table(fire, fibers, outfile, min_coverage=4):
    logging.info("Starting analysis")
    fire = fire.join(fibers, on=["chrom", "fiber", "hap"])
    logging.debug(f"Joined fibers\n{fire}")
    fire_tracks(fire, outfile, min_coverage=min_coverage)
    return 0


def find_nearest(array, value):
    idx = np.searchsorted(array, value, side="left")
    idx[idx < 0] = 0
    idx[idx >= len(array)] = len(array) - 1
    return idx


def write_bed(output_dict, out, first=True):
    # make df
    df = pd.DataFrame(output_dict)
    if first:
        header = True
        mode = "w"
    else:
        header = False
        mode = "a"
    df.to_csv(out, mode=mode, header=header, index=False, sep="\t")


def extra_output_columns(fire, fibers, fdr_table, min_coverage=4):
    return_data = {}
    # get the inital data
    chrom_length = fire.length[0]
    # get coverage for this chromosome and the shuffled fibers
    coverage_array = make_coverage_array(
        fibers.fiber_start.values, fibers.fiber_end.values, chrom_length
    )
    # get the FIRE scores in bed format
    fire_scores = fire_scores_per_chrom(
        fire.start.values,
        fire.end.values,
        fire.fdr.values,
        fire.length.max(),
        coverage_array,
        min_coverage=min_coverage,
    )
    rle_fire_scores = bed_rle(fire_scores)
    # ranges to make calculations on
    starts, ends = (
        rle_fire_scores[:, 0],
        rle_fire_scores[:, 1],
    )
    return_data["#chrom"] = fire.chrom[0]
    return_data["starts"] = starts
    return_data["ends"] = ends

    # get fire info per haplotype
    for hap in [""] + HAPS:
        # select data we are working with
        if hap == "":
            tag = ""
            cur_fire = fire
            cur_fibers = fibers
            cur_rle_fire_scores = rle_fire_scores
            cur_coverage_array = coverage_array
            cur_fire_scores = fire_scores
        else:
            logging.info(f"Processing {hap}")
            tag = f"_{hap}"
            cur_fire = fire[fire.hap == hap]
            if cur_fire.shape[0] == 0:
                for x in [
                    "fire_coverage",
                    "coverage",
                    "score",
                    "FDR",
                    "log_FDR",
                ]:
                    return_data[f"{x}{tag}"] = -1
                continue
            cur_fibers = fibers[fibers.hap == hap]
            cur_coverage_array = make_coverage_array(
                cur_fibers.fiber_start.values, cur_fibers.fiber_end.values, chrom_length
            )
            # get the FIRE scores in bed format
            cur_fire_scores = fire_scores_per_chrom(
                cur_fire.start.values,
                cur_fire.end.values,
                cur_fire.fdr.values,
                cur_fire.length.max(),
                cur_coverage_array,
                min_coverage=min_coverage,
            )
            cur_rle_fire_scores = bed_rle(cur_fire_scores)
        #
        # calculate a bunch of different stats per haplotype
        #
        # fire coverage
        fire_coverage = get_coverage_from_array(
            starts,
            ends,
            make_coverage_array(
                cur_fire.start.values, cur_fire.end.values, chrom_length
            ),
        )
        return_data[f"fire_coverage{tag}"] = fire_coverage

        # total coverage
        coverage = get_coverage_from_array(starts, ends, cur_coverage_array)
        return_data[f"coverage{tag}"] = coverage

        # save the scores
        cur_scores = get_coverage_from_array(starts, ends, cur_fire_scores, stat="max")
        return_data[f"score{tag}"] = cur_scores

        # find the FDRs for the thresholds
        fdr_idx = find_nearest(fdr_table.threshold.values, cur_scores)
        FDRs = fdr_table.FDR.values[fdr_idx]
        return_data[f"FDR{tag}"] = FDRs

        # log the FDRs
        tmp_FDR = FDRs.copy()
        tmp_FDR[tmp_FDR <= 0] = tmp_FDR[tmp_FDR > 0].min()
        log_FDRs = -10 * np.log10(tmp_FDR)
        return_data[f"log_FDR{tag}"] = log_FDRs

        # find local maxima
        if hap == "":
            local_max = argrelextrema(cur_scores, np.greater)
            is_local_max = np.zeros(FDRs.shape[0], dtype=int)
            is_local_max[local_max] = True
            return_data[f"is_local_max{tag}"] = is_local_max

    for key, data in return_data.items():
        if isinstance(data, np.ndarray):
            assert (
                data.shape[0] == starts.shape[0]
            ), f"{key} is not the expected size: {data.shape} instead of {starts.shape}."

    logging.info(f"Finished making data, starting to write")
    return return_data


def write_scores(fire, fibers, fdr_table, outfile, min_coverage=4):
    fire = fire.join(fibers, on=["chrom", "fiber", "hap"])
    first = True
    for chrom, g in fire.groupby("chrom", maintain_order=True):
        logging.info(f"Processing {chrom}")
        # fibers for this chromosome
        fibers = (
            g[["chrom", "fiber", "fiber_start", "fiber_end", "hap"]]
            .unique()
            .to_pandas()
        )
        # convert to pandas for easier manipulation
        g = g.to_pandas()

        # get a bunch of extra columns + per haplotype
        output_dict = extra_output_columns(
            g, fibers, fdr_table, min_coverage=min_coverage
        )

        # write data
        write_bed(output_dict, outfile, first=first)
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
    min_coverage: Optional[int] = 4,
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
        columns=[0, 1, 2, 3, 9, 10],
        new_columns=["chrom", "start", "end", "fiber", "fdr", "hap"],
        comment_char="#",
        n_rows=n_rows,
    )
    logging.debug(f"FIRE peaks {fire}")
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
        columns=[0, 1, 2, 3, 5],
        new_columns=["chrom", "fiber_start", "fiber_end", "fiber", "hap"],
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
        make_fdr_table(fire, fibers, outfile, min_coverage=min_coverage)
    else:
        fdr_table = (
            pl.read_csv(fdr_table_file, separator="\t")
            .to_pandas()
            .sort_values("threshold")
        )
        write_scores(
            fire, fiber_locations, fdr_table, outfile, min_coverage=min_coverage
        )
    return 0


if __name__ == "__main__":
    defopt.run(main, show_types=True, version="0.0.1")
