# Change Log

All notable changes to this project will be documented in this file.
## v0.2.0
- Update fibertools-rs version (also samtools, htslib, and bedtools) in `workflow/envs/env.yaml`
- Fix Polars issue in #52


## v0.1.2

- fix #34
- fix #33
- Update to pull test data from a new s3 bucket

## v0.1.1

Added more informative error messages if an FDR distribution cannot be made or there is not enough coverage.

## v0.1.0

First major release of the FIRE pipeline. This release includes a refactor to reduce the computation by increased use of ft, changes to the output file names to include the fire version among other things, and finally a new launching method for the pipeline that uses pixi. Results are very similar to v0.0.7 of the pipeline; however, there are minor differences in the peak calls and the output names.
