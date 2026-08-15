# Change Log

All notable changes to this project will be documented in this file.

## [0.3.0](https://github.com/fiberseq/FIRE/compare/v0.2.0...v0.3.0) (2026-08-15)


### Features

* per-sample reference support via manifest ref/ref_name columns ([#71](https://github.com/fiberseq/FIRE/issues/71)) ([6e1c24e](https://github.com/fiberseq/FIRE/commit/6e1c24e2b6fe0d81b2df6abdb1ca8b03b690fbbe))


### Refactoring

* remove the unused leviosam2 chain-lift rules ([#70](https://github.com/fiberseq/FIRE/issues/70)) ([150d68e](https://github.com/fiberseq/FIRE/commit/150d68e963f50c441f319a90e8239d52bf45ee70))

## [0.2.0](https://github.com/fiberseq/FIRE/compare/v0.1.3...v0.2.0) (2026-08-14)


### Features

* default min_frac_accessible to 0.10 ([#66](https://github.com/fiberseq/FIRE/issues/66)) ([2f52149](https://github.com/fiberseq/FIRE/commit/2f52149e086fab3a7bd1c26eb74e4f0f487c4591))
* update fibertools-rs to 0.13 (the MA update) and move mac support to Apple Silicon ([#62](https://github.com/fiberseq/FIRE/issues/62)) ([6e8ca80](https://github.com/fiberseq/FIRE/commit/6e8ca80c5d58963cf84beb21177b3915764f0322))


### Bug Fixes

* exclude_from_shuffle when no excludes are configured ([#65](https://github.com/fiberseq/FIRE/issues/65)) ([56966ff](https://github.com/fiberseq/FIRE/commit/56966ffb7a487fe07c2bc49d5425bd6eaae87e83))
* explicitly set schema for reading bed via polars csv ([#53](https://github.com/fiberseq/FIRE/issues/53)) ([fd019e5](https://github.com/fiberseq/FIRE/commit/fd019e5bf9295c04a141057f4e74ca61400eee01)), closes [#52](https://github.com/fiberseq/FIRE/issues/52)
* force chrom as utf8 in pl.read_csv callsites ([#56](https://github.com/fiberseq/FIRE/issues/56)) ([b66739e](https://github.com/fiberseq/FIRE/commit/b66739ea3755068b72baa807540937ea398a8055))
* guard against bioawk phantom empty record in fire_coverage/coverage division ([#48](https://github.com/fiberseq/FIRE/issues/48)) ([835868d](https://github.com/fiberseq/FIRE/commit/835868dee2a56a72986b6ff4559f63385f82a359))

## v0.1.3

- Workaround for a bioawk edge-case failure producing phantom empty records
- Misc fixes: `min_coverage` handling in the Snakefile, env and
  decorated-reads updates (#37)

## v0.1.2

- fix #34
- fix #33
- Update to pull test data from a new s3 bucket

## v0.1.1

Added more informative error messages if an FDR distribution cannot be made or there is not enough coverage.

## v0.1.0

First major release of the FIRE pipeline. This release includes a refactor to reduce the computation by increased use of ft, changes to the output file names to include the fire version among other things, and finally a new launching method for the pipeline that uses pixi. Results are very similar to v0.0.7 of the pipeline; however, there are minor differences in the peak calls and the output names.
