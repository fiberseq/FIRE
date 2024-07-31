# 🔥 **FIRE**: <ins>F</ins>iber-seq <ins>I</ins>nferred <ins>R</ins>egulatory <ins>E</ins>lements

[![DOI](https://zenodo.org/badge/561430995.svg)](https://zenodo.org/doi/10.5281/zenodo.10023811)

A Snakemake workflow for calling Fiber-seq Inferred Regulatory Elements (FIREs) on single molecules. For a more detailed description and methods see the [docs](/docs/README.md), or [watch](https://youtu.be/RiZrMltAiWM?si=sSo64goaNQxgyfcc) my lab meeting on FIRE.

## Install

Please install `snakemake` and all the UCSC Kent utilities. For detailed instructions see the [installation README](/INSTALL.md).

## Configuring

See the [configuration README](/config/README.md), the example [configuration file](/config/config.yaml), and the example [manifest file](/config/config.tbl) for configuration options.


## Run

We have a run script that executes the FIRE snakemake called `fire`, and any extra parameters are passed directly to snakemake. For example:

```bash
./fire --configfile config/config.yaml
```

If you want to do a dry run:

```bash
./fire --configfile config/config.yaml -n
```

If you want to execute across a cluster (modify `profiles/slurm-executor` as needed for distributed execution):

```bash
./fire --configfile config/config.yaml --profile profiles/slurm-executor
```

You can also run snakemake directly, e.g.:

```bash
snakemake \
  --configfile config/config.yaml \
  --profile profiles/slurm-executor \
  --local-cores 8 -k
```

## Test data

You can find input data to test against at [this url](https://s3-us-west-2.amazonaws.com/stergachis-public1/index.html?prefix=FIRE/test-data/).
