# 🔥 **FIRE**: <ins>F</ins>iber-seq <ins>I</ins>nferred <ins>R</ins>egulatory <ins>E</ins>lements

[![DOI](https://zenodo.org/badge/561430995.svg)](https://zenodo.org/doi/10.5281/zenodo.10023811)

A Snakemake workflow for calling Fiber-seq Inferred Regulatory Elements (FIREs) on single molecules. For a more detailed description and methods see the [docs](/docs/README.md), or [watch](https://youtu.be/RiZrMltAiWM?si=sSo64goaNQxgyfcc) my lab meeting on FIRE.

## Install

Please start by installing [pixi](https://pixi.sh/latest/) which handles the environment of the FIRE workflow.

Then install FIRE using `git` and `pixi`:

```bash
git clone https://github.com/fiberseq/FIRE.git
pixi install
```

We then recommend quickly testing your installation by running the test suite:

```bash
pixi run test
```

Further installation instructions can be found in the [INSTALL.md](/INSTALL.md) file.

## Configuring

See the [configuration README](/config/README.md), the example [configuration file](/config/config.yaml), and the example [manifest file](/config/config.tbl) for configuration options.

## Run

The `FIRE` workflow can be executed using the `pixi run fire` command. Under the hood this runs a `snakemake` workflow and any extra parameters are passed directly to snakemake. For example:

```bash
pixi run fire --configfile config/config.yaml
```

If you want to do a dry run:

```bash
pixi run fire --configfile config/config.yaml -n
```

If you want to execute across a cluster (modify `profiles/slurm-executor` as needed for distributed execution):

```bash
pixi run fire --configfile config/config.yaml --profile profiles/slurm-executor
```

You can also run snakemake directly, e.g.:

```bash
pixi shell
snakemake \
  --configfile config/config.yaml \
  --profile profiles/slurm-executor \
  --local-cores 8 -k
```
