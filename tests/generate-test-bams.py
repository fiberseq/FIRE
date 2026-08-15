"""Make the generated test BAMs from the test CRAM.

Two BAMs, each skipped when already present:

- generated/test-chr20.bam: the header lists only chr20. This sample
  tests a bam with fewer contigs than its fasta.
- generated/test-rev.bam: the header lists chr21 before chr20. This
  sample tests header-order preservation, because its header order
  differs from lexicographic order.

Run with the working directory set to fire-test-data (see the test-multi
pixi task).
"""

from pathlib import Path

import pysam

CRAM = "test.cram"
REF = "test.fa.gz"
OUT_DIR = Path("generated")


def remap(read: pysam.AlignedSegment, out: pysam.AlignmentFile) -> pysam.AlignedSegment:
    """Point the read's tids at `out`'s header by contig name."""
    if read.reference_name is not None:
        read.reference_id = out.get_tid(read.reference_name)
    if read.next_reference_id >= 0 and read.next_reference_name is not None:
        read.next_reference_id = out.get_tid(read.next_reference_name)
    return read


def make_chr20() -> None:
    out = OUT_DIR / "test-chr20.bam"
    if out.exists():
        print(f"{out} already present, skipping generation")
        return
    with pysam.AlignmentFile(CRAM, "rc", reference_filename=REF) as cram:
        header = cram.header.to_dict()
        header["SQ"] = [sq for sq in header["SQ"] if sq["SN"] == "chr20"]
        with pysam.AlignmentFile(out, "wb", header=header) as bam:
            for read in cram.fetch("chr20"):
                bam.write(remap(read, bam))
    pysam.index(str(out))
    print(f"wrote {out}")


def make_reversed() -> None:
    out = OUT_DIR / "test-rev.bam"
    if out.exists():
        print(f"{out} already present, skipping generation")
        return
    unsorted = OUT_DIR / "test-rev.unsorted.bam"
    with pysam.AlignmentFile(CRAM, "rc", reference_filename=REF) as cram:
        header = cram.header.to_dict()
        header["SQ"] = list(reversed(header["SQ"]))
        with pysam.AlignmentFile(unsorted, "wb", header=header) as bam:
            for read in cram.fetch(until_eof=True):
                bam.write(remap(read, bam))
    pysam.sort("-o", str(out), str(unsorted))
    unsorted.unlink()
    pysam.index(str(out))
    print(f"wrote {out}")


if __name__ == "__main__":
    OUT_DIR.mkdir(exist_ok=True)
    make_chr20()
    make_reversed()
