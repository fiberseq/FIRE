#!/usr/bin/env bash
# Dry-run validation matrix for the FIRE manifest and reference handling.
# Run with the working directory set to fire-test-data (see the test-dry pixi task).
set -uo pipefail

SNAKEFILE="$PIXI_PROJECT_ROOT/workflow/Snakefile"
CFG="$PIXI_PROJECT_ROOT/tests/config"
FAILURES=0

snk() {
    snakemake -s "$SNAKEFILE" -n -q rules --configfile "$@"
}

expect_pass() {
    local config=$1
    if ! out=$(snk "$config" 2>&1); then
        echo "FAIL (expected pass): $config"
        echo "$out" | tail -5
        FAILURES=$((FAILURES + 1))
    else
        echo "ok (pass): $config"
    fi
}

expect_fail() {
    local config=$1 substring=$2
    if out=$(snk "$config" 2>&1); then
        echo "FAIL (expected failure): $config"
        FAILURES=$((FAILURES + 1))
    elif ! grep -qF "$substring" <<<"$out"; then
        echo "FAIL (wrong message): $config, wanted: $substring"
        echo "$out" | tail -5
        FAILURES=$((FAILURES + 1))
    else
        echo "ok (fail): $config"
    fi
}

expect_stderr() {
    local config=$1 substring=$2
    if ! out=$(snk "$config" 2>&1); then
        echo "FAIL (expected pass): $config"
        FAILURES=$((FAILURES + 1))
    elif ! grep -qF "$substring" <<<"$out"; then
        echo "FAIL (missing stderr line): $config, wanted: $substring"
        FAILURES=$((FAILURES + 1))
    else
        echo "ok (stderr): $config"
    fi
}

# positive cases
expect_pass test.yaml
expect_pass "$CFG/single-4col.yaml"
expect_pass "$CFG/two-sample.yaml"
expect_pass "$CFG/sentinel.yaml"
expect_stderr "$CFG/override-info.yaml" "manifest ref/ref_name columns override config-level values"

# sample-attributed rendered-shell checks for exclude_from_shuffle;
# -R forces the rule so cached results do not hide the rendering. The
# rendered shell contains the sample-scoped output path, so paragraphs
# (awk RS='') can be attributed to one sample.
shell_block() {
    # config, sample -> the rendered exclude_from_shuffle shell for sample
    snakemake -s "$SNAKEFILE" -n -p -R exclude_from_shuffle --configfile "$1" 2>&1 \
        | awk -v RS='' -v sm="results/$2/" '/bedtools genomecov/ && index($0, sm)'
}

check_block() {
    local label=$1 block=$2 must=$3 must_not=$4
    if [ -z "$block" ]; then
        echo "FAIL: no rendered exclude_from_shuffle shell for $label"
        FAILURES=$((FAILURES + 1))
        return
    fi
    if [ -n "$must" ] && ! grep -qF "$must" <<<"$block"; then
        echo "FAIL: $label rendered shell is missing: $must"
        FAILURES=$((FAILURES + 1))
        return
    fi
    if [ -n "$must_not" ] && grep -qF "$must_not" <<<"$block"; then
        echo "FAIL: $label rendered shell must not contain: $must_not"
        FAILURES=$((FAILURES + 1))
        return
    fi
    echo "ok (rendered shell): $label"
}

# with config excludes set: test gets hg38 blacklists + the config exclude,
# test2 gets only the config exclude (a leak of hg38 paths into test2 fails)
block_test=$(shell_block "$CFG/two-sample.yaml" test)
block_test2=$(shell_block "$CFG/two-sample.yaml" test2)
check_block "test (hg38 excludes)" "$block_test" "hg38.blacklist" ""
check_block "test (config exclude)" "$block_test" "extra-exclude.bed" ""
check_block "test2 (no hg38 leak)" "$block_test2" "extra-exclude.bed" "hg38.blacklist"
# the shell template references {params.exclude} twice (the [ -n ] guard
# and gunzip), so each exclude path renders exactly twice; more means a
# get_excludes mutation bug duplicated the list
gap_count=$(grep -oF "hg38.gap.bed.gz" <<<"$block_test" | wc -l | tr -d ' ')
if [ "$gap_count" -ne 2 ]; then
    echo "FAIL: hg38.gap.bed.gz appears $gap_count times for test, expected 2"
    FAILURES=$((FAILURES + 1))
else
    echo "ok (rendered shell): no exclude duplication for test"
fi

# without config excludes: test2 exercises the truly-empty excludes branch
# (the [ -n ] guard must render with an empty parameter)
block_empty=$(shell_block "$CFG/no-excludes.yaml" test2)
check_block "test2 (empty excludes)" "$block_empty" "" "gunzip -cf ."
if grep -qF 'if [ -n "" ]' <<<"$block_empty"; then
    echo "ok (rendered shell): empty excludes guard renders for test2"
else
    echo "FAIL: empty excludes guard did not render empty for test2"
    FAILURES=$((FAILURES + 1))
fi


# error cases, one per validator branch; substrings are specific enough
# that a wrong error cannot satisfy the assertion
expect_fail "$CFG/err-no-ref-anywhere.yaml" "no reference specified"
expect_fail "$CFG/err-ref-col-only.yaml" "manifest columns 'ref' and 'ref_name' must be provided together"
expect_fail "$CFG/err-config-ref-only.yaml" "config options 'ref' and 'ref_name' must be provided together"
expect_fail "$CFG/err-nan-cell.yaml" "missing or malformed manifest fields"
expect_fail "$CFG/err-extra-column.yaml" "cannot parse manifest"
expect_fail "$CFG/err-sentinel-no-config.yaml" "is not set in config.yaml"
expect_fail "$CFG/err-dup-sample.yaml" "duplicate sample names"
expect_fail "$CFG/err-missing-refpath.yaml" "reference file"
expect_fail "$CFG/err-missing-fai.yaml" "reference index file"
expect_fail "$CFG/err-missing-bam.yaml" "cannot read input bam"
expect_fail "$CFG/err-keepchrs.yaml" "no chromosomes left"

if [ "$FAILURES" -gt 0 ]; then
    echo "dry-run-checks: $FAILURES failure(s)"
    exit 1
fi
echo "dry-run-checks: all checks passed"
