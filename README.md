# 🔥 **FIRE**: <ins>F</ins>iber-seq <ins>I</ins>nferred <ins>R</ins>egulatory <ins>E</ins>lements

[![DOI](https://zenodo.org/badge/561430995.svg)](https://zenodo.org/doi/10.5281/zenodo.10023811)

A Snakemake workflow for calling Fiber-seq Inferred Regulatory Elements (FIREs) on single molecules. For a more detailed description and methods see the [docs](/docs/README.md).

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

## Test data

You can find input data to test against at [this url](https://s3-us-west-2.amazonaws.com/stergachis-public1/index.html?prefix=Projects/Phased-GM12878/fire-test/).
