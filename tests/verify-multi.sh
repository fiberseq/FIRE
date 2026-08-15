#!/usr/bin/env bash
# Checks for the multi-sample run (four samples, per-sample references).
# Run with the working directory set to fire-test-data after `pixi run test-multi`.
set -euo pipefail

V="v$(echo "$PIXI_PROJECT_VERSION" | cut -d. -f1-2)"
FAILURES=0

assert_genome() {
    local sm=$1 name=$2
    local genomes="results/$sm/trackHub-$V/genomes.txt"
    if ! grep -qx "genome $name" "$genomes"; then
        echo "FAIL: $genomes does not contain 'genome $name'"
        FAILURES=$((FAILURES + 1))
    else
        echo "ok: $sm uses genome $name"
    fi
}

assert_genome test hg38
assert_genome test2 testref
assert_genome test20 testref2
assert_genome test-rev hg38

# the chr20-only sample must have no chr21 peaks; grep -c reads the whole
# stream, so no SIGPIPE can mask a hit under pipefail
peaks20="results/test20/test20-fire-$V-peaks.bed.gz"
chr21_count=$(gunzip -c "$peaks20" | cut -f 1 | { grep -cx "chr21" || true; })
if [ "$chr21_count" -gt 0 ]; then
    echo "FAIL: $peaks20 contains $chr21_count chr21 records"
    FAILURES=$((FAILURES + 1))
else
    echo "ok: test20 output is chr20 only"
fi

# the reversed-header sample must keep bam header order (chr21 first) in
# concatenated per-chromosome outputs
pileup_rev="results/test-rev/test-rev-fire-$V-pileup.bed.gz"
# || true absorbs the SIGPIPE that head sends up the pipeline
first_chrom=$(gunzip -c "$pileup_rev" | grep -v "^#" | head -n 1 | cut -f 1 || true)
if [ "$first_chrom" != "chr21" ]; then
    echo "FAIL: $pileup_rev starts with $first_chrom, expected chr21 (header order)"
    FAILURES=$((FAILURES + 1))
else
    echo "ok: test-rev output preserves bam header order"
fi

# test-rev holds the same reads as test with a reversed header, so output
# CONTENT must be identical modulo row order — the only assertion that can
# catch a silent (non-crashing) ordering bug
assert_same_content() {
    local a=$1 b=$2 label=$3
    if diff <(gunzip -c "$a" | grep -v "^#" | LC_ALL=C sort) \
        <(gunzip -c "$b" | grep -v "^#" | LC_ALL=C sort) >/dev/null; then
        echo "ok: $label content identical between test and test-rev"
    else
        echo "FAIL: $label content differs between test and test-rev"
        FAILURES=$((FAILURES + 1))
    fi
}
assert_same_content \
    "results/test/test-fire-$V-peaks.bed.gz" \
    "results/test-rev/test-rev-fire-$V-peaks.bed.gz" "peaks"
assert_same_content \
    "results/test/test-fire-$V-pileup.bed.gz" \
    "results/test-rev/test-rev-fire-$V-pileup.bed.gz" "pileup"

if [ "$FAILURES" -gt 0 ]; then
    echo "verify-multi: $FAILURES failure(s)"
    exit 1
fi
echo "verify-multi: all checks passed"
