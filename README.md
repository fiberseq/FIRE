# 🔥 **FIRE**: <ins>F</ins>iber-seq <ins>I</ins>mplicated <ins>R</ins>egulatory <ins>E</ins>lements
A pipeline for calling Fiber-seq Implicated Regulatory Elements (FIREs) on single molecules.

## Install
```bash
# make fibertools env for running snakemake
mamba env update -n fibertools --file env.yaml

# activate env
conda activate fibertools

# partially test install
fibertools -h
```

To get the latest fibertools without rebuilding the whole env you can do:
```bash
yes | pip uninstall fibertools && pip install git+ssh://git@github.com/fiberseq/fibertools.git && fibertools model -h
```
