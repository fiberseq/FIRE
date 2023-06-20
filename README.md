# 🔥 **FIRE**: <ins>F</ins>iber-seq <ins>I</ins>mplicated <ins>R</ins>egulatory <ins>E</ins>lements
A pipeline for calling Fiber-seq Implicated Regulatory Elements (FIREs) on single molecules.

## Configuring
See `config/config.yaml` and `config/config.tbl` for configuration options.

## Install
Add a snakemake conda prefix to your `bashrc`, e.g. in the Stergachis lab add:
```
export SNAKEMAKE_CONDA_PREFIX=/mmfs1/gscratch/stergachislab/snakemake-conda-envs
```
Then snakemake should install all the additional requirements as a conda env in that directory. 

## Run
```bash
snakemake \
  --configfile config/config.yaml \
  --local-cores $(nproc) \
  --cores $(nproc) \
  --rerun-triggers mtime \
  --rerun-incomplete \
  --scheduler greedy \
  --resources load=1000 \
  --show-failed-logs \
  -p 
```
And modify as needed for distributed execution. 


### Updates (ignore)
Then activate the env and update `ft` to the latest version of `fibertools-rs`:
```bash
mamba install -c conda-forge -c bioconda 'fibertools-rs>=0.2.0' && ft --version
```
Then update to the latest python fibertools:
```bash
yes | pip uninstall fibertools && pip install git+ssh://git@github.com/mrvollger/fibertools.git; fibertools -h
```

