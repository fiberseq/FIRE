#!/usr/bin/env bash
# Golden-count checks for the single-sample regression test.
# Run with the working directory set to fire-test-data after `pixi run test`.
set -euo pipefail

V="v$(echo "$PIXI_PROJECT_VERSION" | cut -d. -f1-2)"
PEAKS="results/test/test-fire-$V-peaks.bed.gz"
ELEMENTS="results/test/additional-outputs-$V/fire-peaks/test-$V-fire-elements.bed.gz"
GENOMES="results/test/trackHub-$V/genomes.txt"

peaks=$(gunzip -c "$PEAKS" | wc -l | tr -d ' ')
elements=$(gunzip -c "$ELEMENTS" | wc -l | tr -d ' ')
echo "peaks file lines: $peaks (expected 1181)"
echo "fire elements: $elements (expected 46223)"

FAILURES=0
[ "$peaks" -eq 1181 ] || FAILURES=$((FAILURES + 1))
[ "$elements" -eq 46223 ] || FAILURES=$((FAILURES + 1))
if ! grep -qx "genome hg38" "$GENOMES"; then
    echo "FAIL: $GENOMES does not contain 'genome hg38'"
    FAILURES=$((FAILURES + 1))
fi

if [ "$FAILURES" -gt 0 ]; then
    echo "verify-results: $FAILURES failure(s)"
    exit 1
fi
echo "verify-results: all checks passed"
