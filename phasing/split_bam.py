#!/usr/bin/env python3
import argparse
import pysam 
import pandas as pd
import tqdm
import logging

def parse():
    """Console script for fibertools."""
    parser = argparse.ArgumentParser(
        description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("bam", help="aligned bam file from actc")
    parser.add_argument("--reads", nargs="+", help="read lists")
    parser.add_argument("-t", "--threads", help="n threads to use", type=int, default=16)
    parser.add_argument(
        "-v", "--verbose", help="increase logging verbosity", action="store_true"
    )
    args = parser.parse_args()
    log_format = "[%(levelname)s][Time elapsed (ms) %(relativeCreated)d]: %(message)s"
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(format=log_format, level=log_level)
    return args


def main():
    args = parse()
    bam = pysam.AlignmentFile(args.bam, threads=args.threads, check_sq=False)
    outs = []
    for readset in args.reads:
        x = {line.strip() for line in open(readset).readlines()}
        out = pysam.AlignmentFile(f"{readset}.bam", "wb", template=bam, threads=args.threads)
        outs.append((out, x))
        logging.info(f"{readset} had {len(x)} reads.")

    for rec in tqdm.tqdm(bam.fetch(until_eof=True)):
        for out, readset in outs:
            if rec.query_name in readset:
                readset.remove(rec.query_name)
                out.write(rec)
                break

    for out, readset in outs:
        logging.info(f"{len(readset)} reads from {out} were unplaced.")
        out.close()
    bam.close()

if __name__ == "__main__":
    main()



