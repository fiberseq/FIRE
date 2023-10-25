# 🔥 **FIRE**: <ins>F</ins>iber-seq <ins>I</ins>nferred <ins>R</ins>egulatory <ins>E</ins>lements

[![DOI](https://zenodo.org/badge/561430995.svg)](https://zenodo.org/doi/10.5281/zenodo.10023811)

A Snakemake workflow for calling Fiber-seq Inferred Regulatory Elements (FIREs) on single molecules.

## Configuring

See `config/config.yaml` and `config/config.tbl` for configuration options.

## Install

You will need **snakemake** and all the **UCSC Kent utilities** and the **latest version** of them (v455).

Add a snakemake conda prefix to your `bashrc`, e.g. in the Stergachis lab add:

```bash
export SNAKEMAKE_CONDA_PREFIX=/mmfs1/gscratch/stergachislab/snakemake-conda-envs
```

Then snakemake installs all the additional requirements as conda envs in that directory.

## Run

We have a run script that executes the FIRE snakemake called `fire`, and any extra parameters are passed directly to snakemake. For example:

```bash
./fire --configfile config/config.yaml
```

If you want to do a dry run:

```bash
./fire --configfile config/config.yaml -n
```

If you want to execute across a cluster (modify `profiles/compute` as needed for distributed execution):

```bash
./fire --configfile config/config.yaml --profile profiles/compute
```

You can also run snakemake directly, e.g.:

```bash
snakemake \
  --configfile config/config.yaml \
  --profile profiles/compute \
  --local-cores 8 -k
```

# Outputs

| Directory     | File                                            | Description                                                             |
| ------------- | ----------------------------------------------- | ----------------------------------------------------------------------- |
| coverage      |                                                 |                                                                         |
|               | {all,hap1,hap2}\_element_coverages.bed.gz       | Coverage tracks for FIREs, nucleosomes, and linkers                     |
|               | exclude-from-shuffles.bed.gz                    | Regions to exclude when making null distributions                       |
|               | {sample}.{bed.gz,d4}                            | bedgraph of coverages                                                   |
|               | {sample}.{median, maximum,minimum}.coverage.txt | Allowed coverage range for analysis                                     |
| fiber-calls   |                                                 |                                                                         |
|               | FIRE.bed.gz                                     | Every FIRE element in every FIRE                                        |
|               | model.results.bed.gz                            | Every FIRE, linker, and nucleosome in every fiber                       |
|               | fire-fibers.bed.gz                              | bed12 start and end of every fiber                                      |
|               | fire-fiber-decorators.bed.gz                    | decorator file that adds annotations of FIRE elements to the bed12 file |
| FDR-peaks     |                                                 |                                                                         |
|               | FDR-FIRE-peaks.bed.gz                           | Fiber-seq peak calls                                                    |
|               | FDR-wide-peaks.bed.gz                           | Fiber-seq wide peak calls                                               |
|               | FDR.track.bed.gz                                | Track of FDR significance of accessibility                              |
|               | {sm}.peaks-vs-percent.pdf                       | Number of peaks vs % accessible                                         |
| all/hap1/hap2 |                                                 |                                                                         |
|               | percent.accessible.bed.gz                       | % of (haplotype) fibers that are accessible                             |
| hap1-vs-hap2  |                                                 |                                                                         |
|               | FIRE.hap.differences.bed                        | Large table of FIREs that are different between hap1 and hap2           |
|               | hap1-vs-hap2-volcano.pdf                        | Volcano plot of FIREs that are different between hap1 and hap2          |
|               | hap1-vs-hap2.pdf                                | Scatter plot of FIREs and their percent accessibility for each hap      |
| trackHub      |                                                 |                                                                         |
|               | \*                                              | a trackHub directory that can be loaded into the UCSC browser           |

---

## Test data

You can find input data to test against at [this url](https://s3-us-west-2.amazonaws.com/stergachis-public1/index.html?prefix=Projects/Phased-GM12878/fire-test/).
