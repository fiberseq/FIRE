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

ROLLING_FIRE_SCORE_WINDOW_SIZE = 200

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
def is_local_max(array):
    output = []
    for idx in range(array.shape[0]):
        if idx - 1 < 0 or idx + 1 >= array.shape[0]:
            output.append(False)
            continue
        cur_res = False
        pre = array[idx - 1]
        cur = array[idx]
        next = array[idx + 1]
        if cur >= pre and cur >= next:
            cur_res = True
        output.append(cur_res)

    return output


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
    fire_scores = np.zeros(int(chrom_length), dtype=np.float64)

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
    coverage_array = np.zeros(int(chrom_length), dtype=np.float64)
    for start, end in zip(starts, ends):
        coverage_array[start:end] += 1
    return coverage_array


def fire_tracks(fire, outfile, min_coverage=4):
    null_s = []
    fire_s = []
    logging.info(f"Fire data\n{fire}")
    # number of elements where start and fiber_start are null
    null_count = fire.filter(
        pl.col("start").is_null() & pl.col("fiber_start").is_null()
    ).shape[0]
    if null_count > 0:
        logging.warn(f"Null count: {null_count}")

    for chrom, g in fire.group_by("chrom", maintain_order=True):
        logging.info(f"Processing {chrom}")
        # fibers for this chromosome
        fibers = (
            g[FIBER_COLUMNS]
            .filter(~pl.col("fiber_start").is_null())
            .unique()
            .to_pandas()
        )
        # convert to pandas for easier manipulation
        g = (
            g.filter(~pl.col("start").is_null())
            .filter(~pl.col("fiber_start").is_null())
            .to_pandas()
        )
        logging.debug(f"Grouped fire data\n{g}\n{g.dtypes}")

        if g.shape[0] == 0:
            logging.warning(f"No data for {chrom}")
            continue

        # get coverage for this chromosome and the shuffled fibers
        chrom_length = g.length[0].astype(int)
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


def make_fdr_table(fire, outfile, min_coverage=4):
    logging.info("Starting analysis")
    logging.debug(f"Joined fibers\n{fire}")
    fire_tracks(fire, outfile, min_coverage=min_coverage)
    return 0


def find_nearest(array, value):
    idx = np.searchsorted(array, value, side="left")
    idx[idx < 0] = 0
    idx[idx >= len(array)] = len(array) - 1
    return idx


def write_bed(chrom, output_dict, out, first=True):
    chrom_length = output_dict["coverage"].shape[0]
    # make df
    if first:
        header = True
        mode = "w"
    else:
        header = False
        mode = "a"
    logging.info("Making data frame")
    df = pl.DataFrame(output_dict)
    del output_dict
    gc.collect()

    # finding maxes within windows
    df = df.with_columns(
        max_window_score=pl.col("score").rolling_max(
            window_size=ROLLING_FIRE_SCORE_WINDOW_SIZE,
            center=True,
        )
    )
    original_columns = df.columns
    # find and clear the duplicates
    logging.info("Finding duplicates")
    # float array that says if a row is different from the previous row
    diff = (((df != df.shift(periods=1)).sum(axis=1)) > 0) * 1.0
    # turn the diff array into a group number
    logging.info("Merging duplicates")
    df = (
        df.with_columns(
            diff.cumsum().alias("group"),
        )
        .with_row_count(name="end", offset=1)
        .unique(keep="last", subset=original_columns + ["group"], maintain_order=True)
        .with_columns(
            pl.col("end").shift_and_fill(periods=1, fill_value=0).alias("start"),
            pl.lit(chrom).alias("#chrom"),
        )
        .select(["#chrom", "start", "end"] + original_columns)
        .sort(["#chrom", "start", "end"])
    ).to_pandas()

    # can only find local maxes after de duplicating
    # window_score = (
    #    df.rolling(ROLLING_FIRE_SCORE_WINDOW_SIZE, on="start", center=True)
    #    .score.max()
    #    .values
    # )
    # df["is_local_max"] = df["score"] == window_score  # df["score"].values
    df["is_local_max"] = (df["score"] == df["max_window_score"]) & (df["score"] > 0.0)

    logging.info(f"Found {df.is_local_max.sum():,} local maximums.")

    # checks
    final_end = df.end.max()
    assert final_end == chrom_length, f"{final_end} != {chrom_length}"

    logging.info(f"Writing {chrom}")
    df.to_csv(out, mode=mode, header=header, index=False, sep="\t")
    logging.info(f"Done writing {chrom}")
    return


def extra_output_columns(fire, fibers, fdr_table, min_coverage=4):
    return_data = {}
    # get the inital data
    chrom_length = fire.length[0]

    # get fire info per haplotype
    for hap in [""] + HAPS:
        # select data we are working with
        if hap == "":
            logging.info("Processing all haplotypes")
            tag = ""
            cur_fire = fire
            cur_fibers = fibers
        else:
            logging.info(f"Processing {hap}")
            tag = f"_{hap}"
            cur_fibers = fibers[fibers.hap == hap]
            cur_fire = fire[fire.hap == hap]
            # if no data in the hap write empty values
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
        # fire coverage
        fire_coverage = make_coverage_array(
            cur_fire.start.values, cur_fire.end.values, chrom_length
        )
        return_data[f"fire_coverage{tag}"] = fire_coverage

        # total coverage
        return_data[f"coverage{tag}"] = cur_coverage_array

        # save the scores
        return_data[f"score{tag}"] = cur_fire_scores

        # find the FDRs for the thresholds
        fdr_idx = find_nearest(fdr_table.threshold.values, cur_fire_scores)
        FDRs = fdr_table.FDR.values[fdr_idx]
        return_data[f"FDR{tag}"] = FDRs

        # log the FDRs
        tmp_FDR = FDRs.copy()
        tmp_FDR[tmp_FDR <= 0] = tmp_FDR[tmp_FDR > 0].min()
        log_FDRs = -10 * np.log10(tmp_FDR)
        return_data[f"log_FDR{tag}"] = log_FDRs

    for key, data in return_data.items():
        if isinstance(data, np.ndarray):
            assert (
                data.shape[0] == chrom_length
            ), f"{key} is not the expected size: {data.shape} instead of {chrom_length.shape}."
    return return_data


def write_scores(fire, fdr_table, outfile, min_coverage=4):
    first = True
    for chrom, g in fire.group_by("chrom", maintain_order=True):
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
        write_bed(chrom, output_dict, outfile, first=first)
        first = False


def main(
    infile: Path,
    fiber_locations_file: Path,
    genome_file: Path,
    *,
    outfile: Optional[Path] = None,
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
        fiber_locations = fiber_locations.join(
            shuffled_locations, on=["chrom", "fiber"]
        )

    logging.info("Joining FIRE elements and fibers and then sorting")
    fire = fire.join(fiber_locations, on=["chrom", "fiber", "hap"], how="outer").sort(
        ["chrom", "start", "end"]
    )

    if shuffled_locations_file is not None:
        make_fdr_table(fire, outfile, min_coverage=min_coverage)
    else:
        fdr_table = (
            pl.read_csv(fdr_table_file, separator="\t")
            .to_pandas()
            .sort_values("threshold")
        )
        write_scores(fire, fdr_table, outfile, min_coverage=min_coverage)
    return 0


if __name__ == "__main__":
    defopt.run(main, show_types=True, version="0.0.1")
