# Install

You will need **snakemake** and all the **UCSC Kent utilities** (version >= 455).

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
export APPTAINER_CACHEDIR=/mmfs1/gscratch/stergachislab/snakemake-conda-envs/apptainer-cache
```
Then snakemake installs all the additional requirements as conda envs in that directory.
