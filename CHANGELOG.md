# Change Log

All notable changes to this project will be documented in this file.

## v0.2.0

- `min_frac_accessible` now defaults to 0.10: FIRE peaks must have at least 10%
  of reads accessible in addition to passing FDR peak calling (set to 0.0 to
  restore the old behavior)
- Update fibertools-rs to 0.13.x (requires `ft pileup --rgn`; also updates samtools,
  htslib, and bedtools pins in `workflow/envs/env.yaml`)
- Polars fixes: explicit schemas for chrom/score columns (#48, #52, #53, #56)
- `exclude_from_shuffle` now works when no `excludes` are configured (non-hg38 refs)
- Test data now downloaded via rclone; test runs start from a clean state
  (`test-clean` task added)

## v0.1.2

- fix #34
- fix #33
- Update to pull test data from a new s3 bucket

## v0.1.1

Added more informative error messages if an FDR distribution cannot be made or there is not enough coverage.

## v0.1.0

First major release of the FIRE pipeline. This release includes a refactor to reduce the computation by increased use of ft, changes to the output file names to include the fire version among other things, and finally a new launching method for the pipeline that uses pixi. Results are very similar to v0.0.7 of the pipeline; however, there are minor differences in the peak calls and the output names.
