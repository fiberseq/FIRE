# 🔥 **FIRE**: <ins>F</ins>iber-seq <ins>I</ins>nferred <ins>R</ins>egulatory <ins>E</ins>lements

[![DOI](https://zenodo.org/badge/561430995.svg)](https://zenodo.org/doi/10.5281/zenodo.10023811)

A Snakemake workflow for calling Fiber-seq Inferred Regulatory Elements (FIREs) on single molecules. For a more detailed description and methods see the [docs](/docs/README.md), or [watch](https://youtu.be/RiZrMltAiWM?si=sSo64goaNQxgyfcc) my lab meeting on FIRE.


## Install

You will need **snakemake** and all the **UCSC Kent utilities** and the **latest version** of them (v455).

You can install snakemake using conda/mamba, e.g.:
```
mamba create -c conda-forge -c bioconda -n snakemake 'snakemake>=8.4'
```

You can find the UCSC kent utilities at [this url](http://hgdownload.soe.ucsc.edu/admin/exe/). You will need to add the directory containing the utilities to your `PATH` environment variable.

Finally, if you wish to distribute jobs across a cluster you will need to install the appropriate [snakemake executor plugin](https://snakemake.github.io/snakemake-plugin-catalog/). For example, to use SLURM you can install the `snakemake-executor-slurm` plugin using pip:
```  
pip install snakemake-executor-plugin-slurm
```

We recommend adding a snakemake conda prefix to your `bashrc`, e.g. in the Stergachis lab add:
```bash
export SNAKEMAKE_CONDA_PREFIX=/mmfs1/gscratch/stergachislab/snakemake-conda-envs
```
Then snakemake installs all the additional requirements as conda envs in that directory.

## Configuring

See `config/config.yaml` and `config/config.tbl` for configuration options.


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

You can find input data to test against at [this url](https://s3-us-west-2.amazonaws.com/stergachis-public1/index.html?prefix=Projects/Phased-GM12878/fire-test/).
