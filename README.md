# 🔥 **FIRE**: <ins>F</ins>iber-seq <ins>I</ins>mplicated <ins>R</ins>egulatory <ins>E</ins>lements
A pipeline for calling Fiber-seq Implicated Regulatory Elements (FIREs) on single molecules.

## Configuring
See `config/config.yaml` and `config/config.tbl` for configuration options.

## Install
Go to the [fiberseq-smk repo](https://github.com/fiberseq/fiberseq-smk) and follow the instructions to install fiberseq-smk.
Then activate the env and update `ft` to the latest version of `fibertools-rs`:
```bash
mamba install -c conda-forge -c bioconda fibertools-rs>=0.2.0 && ft --version
```
Then update to the latest python fibertools:
```bash
yes | pip uninstall fibertools && pip install git+https://github.com/fiberseq/fibertools && fibertools -h 
```

## Model
Unless directed otherwise it would be best to use this model for your data:
```
yes | pip uninstall fibertools && pip install git+https://github.com/fiberseq/fibertools && fibertools -h 
```

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