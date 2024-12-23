#!/usr/bin/env python
import defopt
import logging
import polars as pl
import io
import sys


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
        overlap_bp = min(pre_en, cur_en) - max(pre_st, cur_st)
        reciprocal_overlap = min(
            overlap_bp / (pre_en - pre_st), overlap_bp / (cur_en - cur_st)
        )
        if (
            overlap_frac >= min_frac_overlap
            or reciprocal_overlap >= min_reciprocal_overlap
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
        df.sort(["#chrom", "start"])
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
            (~pl.col("shares_FIREs")).cum_sum().alias("group"),
        )
        .sort("group")
        .with_columns(
            pl.col("score").max().over("group").name.suffix("_max"),
            peak_start=pl.col("peak_start").min().over("group").cast(pl.UInt32),
            peak_end=pl.col("peak_end").max().over("group").cast(pl.UInt32),
            local_max_count=pl.col("group").len().over("group"),
            # FIRE_IDs=pl.col("FIRE_IDs").flatten().over("group"),
        )
        .filter(pl.col("score") == pl.col("score_max"))
        # filter multiple maxes
        .group_by("group")
        .agg(pl.all().head(1))
        .explode(pl.all().exclude("group"))
        # add the peak length
        .with_columns(
            peak_length=pl.col("peak_end") - pl.col("peak_start"),
        )
    )
    return df


def iterative_merge(
    df, min_frac_overlap=0.5, min_reciprocal_overlap=0.75, max_grouping_iterations=10
):
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
        if min_frac_overlap >= 1.0:
            logging.info(
                f"{min_reciprocal_overlap:.0%} reciprocal overlap merging round {i} is done: {df.shape[0]:,}"
            )
        else:
            logging.info(
                f"Merging when {min_frac_overlap:.0%} of FIRE elements are shared. Round {i} is done. {n_row:,}"
            )
    return df


def main(
    *,
    max_score_every: int = None,
    min_frac_overlap: float = 0.5,
    min_reciprocal_overlap: float = 0.90,
    min_frac_accessible: float = 0.0,
    max_grouping_iterations: int = 10,
    min_cov: int = 0,
    max_cov: int = 100_000_000_000,
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
    df = pl.read_csv(inf, separator="\t", null_values=".")
    if df.shape[0] == 0:
        logging.info("No peaks to merge")
        return 0
    df = df.with_columns(
        FIRE_IDs=pl.col("FIRE_IDs").str.split(",").cast(pl.List(pl.UInt32)),
    )

    if max_score_every is not None:
        df = (
            # group into sliding X bp windows and only keep the highest score
            df.with_columns(roll_start=pl.col("start"))
            .sort(["#chrom", "roll_start"])
            .groupby_rolling("roll_start", period=f"{max_score_every}i", by="#chrom")
            .agg([pl.exclude("roll_start").sort_by("score").last()])
            .drop("roll_start")
            # remove any peaks that are the highest score for multiple X bp windows
            .unique(subset=["#chrom", "peak_start", "peak_end", "start", "end"])
            .sort(["#chrom", "peak_start"])
        )
    logging.info(
        f"Dynamic window merging over {max_score_every} bp is done: {df.shape[0]:,}"
    )
    # reciprocal overlap merging
    df = iterative_merge(
        df,
        min_frac_overlap=2.0,
        min_reciprocal_overlap=min_reciprocal_overlap,
        max_grouping_iterations=max_grouping_iterations,
    )
    #  fire overlap merging
    df = iterative_merge(
        df,
        min_frac_overlap=min_frac_overlap,
        min_reciprocal_overlap=2.0,
        max_grouping_iterations=max_grouping_iterations,
    )
    # one final round of reciprocal overlap merging
    df = iterative_merge(
        df,
        min_frac_overlap=1.0,
        min_reciprocal_overlap=min_reciprocal_overlap,
        max_grouping_iterations=2,
    )
    # add a column indicating if the peak passes coverage filters
    df = df.with_columns(
        pass_coverage=(pl.col("coverage") >= min_cov) & (pl.col("coverage") <= max_cov),
    ).filter(pl.col("fire_coverage") / pl.col("coverage") >= min_frac_accessible)

    # write to stdout
    (
        df.sort(["#chrom", "peak_start", "peak_end"])
        .drop("score_max", "group", "FIRE_IDs", "shares_FIREs", "is_local_max")
        .write_csv("/dev/stdout", separator="\t")
    )
    return 0


if __name__ == "__main__":
    defopt.run(main, show_types=True, version="0.0.1")
