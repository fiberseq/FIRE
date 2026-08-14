# FIRE test fixtures and scripts

Everything in this directory is test material. Nothing here is an example
configuration for real use — see `config/` for that. All configs reference
the downloaded `fire-test-data/` directory and run with that directory as
the working directory.

## Entry points (pixi tasks)

| Task | What it does |
|---|---|
| `pixi run test-dry` | Runs `dry-run-checks.sh`: the full validation matrix below as fast `snakemake -n` checks. Runs first in CI. |
| `pixi run test` | Full single-sample regression run (unchanged, uses `fire-test-data/test.yaml`). |
| `pixi run test-verify` | Runs `verify-results.sh`: golden output counts for the regression run. Runs in CI after `test`. |
| `pixi run test-multi` | Local only. Generates a chr20-only bam and a reversed-header bam (`generate-test-bams.py`), runs the four-sample `config/multi.yaml`, then checks it with `verify-multi.sh`. |

## Fixtures in `config/`

Positive cases (the dry run must pass):

| Fixture | Exercises |
|---|---|
| `single-4col.yaml` + `.tbl` | Four-column manifest, no config-level reference. |
| `two-sample.yaml` + `.tbl` | Two samples, one with `ref_name: hg38` (auto-excludes branch) and one without (empty-excludes branch). Sets `excludes:` to guard the list-copy bug in `get_excludes`. |
| `sentinel.yaml` + `.tbl` | A `.` cell in `ref`/`ref_name` falls back to the config value. |
| `override-info.yaml` | Manifest columns and config values both set; asserts the override INFO line. |
| `multi.yaml` + `.tbl` | Four samples for `test-multi`: `generated/test-chr20.bam` (bam contigs are a subset of the fasta), `generated/test-rev.bam` (header order chr21,chr20 — checks that outputs keep bam header order), and a sentinel row. |
| `no-excludes.yaml` | Two samples with no `excludes:` key; the rendered shell must show the truly-empty excludes branch. |
| `extra-exclude.bed` | Small exclude file used by `two-sample.yaml` and `multi.yaml`. |

Error cases (`err-*`, the dry run must fail with a specific message):

| Fixture | Validator branch |
|---|---|
| `err-no-ref-anywhere.yaml` | No reference in the manifest or the config. |
| `err-ref-col-only.yaml` + `.tbl` | Manifest has `ref` without `ref_name`. |
| `err-config-ref-only.yaml` | Config has `ref` without `ref_name`. |
| `err-nan-cell.yaml` + `.tbl` | Short manifest row (missing cells). |
| `err-sentinel-no-config.yaml` + `.tbl` | `.` cell with no config value to fall back to. |
| `err-dup-sample.yaml` + `.tbl` | Duplicate sample names. |
| `err-missing-refpath.yaml` + `.tbl` | Reference fasta does not exist. |
| `err-missing-fai.yaml` + `.tbl` + `no-fai.fa` | Fasta exists but has no `.fai`. |
| `err-missing-bam.yaml` + `.tbl` | Input bam does not exist. |
| `err-keepchrs.yaml` | `keep_chromosomes` filters out every chromosome. |
