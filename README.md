# fiberseq-fdr
A pipeline for making fiberseq FDR calls for MSPs

## install
```bash
# make fibertools env for running snakemake
mamba env update -n fibertools --file env.yaml

# activate env
conda activate fibertools

# partially test install
fibertools -h
```

To get the latest fibertools without rebuiling the whole env you can do:
```bash
yes | pip uninstall fibertoolds && pip install git+ssh://git@github.com/mrvollger/fibertools.git; fibertools -h
```
