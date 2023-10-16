# 🔥 **FIRE**: <ins>F</ins>iber-seq <ins>I</ins>nferred <ins>R</ins>egulatory <ins>E</ins>lements

A Snakemake workflow for calling Fiber-seq Inferred Regulatory Elements (FIREs) on single molecules.

## Configuring

See `config/config.yaml` and `config/config.tbl` for configuration options.

## Install

Add a snakemake conda prefix to your `bashrc`, e.g. in the Stergachis lab add:

```bash
export SNAKEMAKE_CONDA_PREFIX=/mmfs1/gscratch/stergachislab/snakemake-conda-envs
```

Then snakemake should install all the additional requirements as a conda env in that directory.

## Model

Unless directed otherwise it would be best to use this model for your data:

```bash
models/model.dat
```

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

```bash
FIRE.bed.gz # Every FIRE element on every fiber
FIRE.peaks.with.coverage.bed # Every FIRE peak in the genome, filtered for coverage
percent-in-clusters.txt # Percent of FIREs that are in clusters (QC metric)
{sample name}.peaks-vs-percent.pdf # Plot of of the above
clustering-vs-null.bed.gz # A null distribution of FIRE.bed.gz for measuring clustering
coverage/ # Whole genome coverage tracks
___ {sample name}.bed.gz # bed graph of coverage
___ {sample name}.d4 # bed graph of coverage in d4 format
___ {sample name}.median.chromosome.coverage.bed # median coverage per chromosome
___ {sample name}.median.coverage.txt # median coverage across the genome ignoring regions with zero coverage
all/ # results for all fibers
___ acc.model.results.bed.gz # All nucleosomes and MSPs with their FIRE scores for all fibers
___ acc.model.results.bed.gz.tbi # index for above
___ peak.calls.bed # FIRE peaks without coverage or significance filtering
___ fire.peaks.and.coverages.bed.gz # FIRE peaks with coverage but no significance filtering
___ fire.peaks.and.coverages.bed.gz.tbi # index for above
hap1/ # results for hap1 fibers, same structure as the all directory
hap2/ # results for hap2 fibers, same structure as the all directory
unk/ # results for fibers from an unknown haplotype, same structure as the all directory
hap1-vs-hap2 # results for comparing hap1 and hap2 fibers
___ FIRE.hap.differences.bed # Large table of FIREs that are different between hap1 and hap2
___ FIRE.hap.differences.bed9 # Same as above but in bed9 format
___ hap1-vs-hap2-volcano.pdf # Volcano plot of FIREs that are different between hap1 and hap2
___ hap1-vs-hap2.pdf # Scatter plot of FIREs and their percent accessibility for each hap
trackHub/ #
```

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
|               | asdf                                            |                                                                         |
|               | asdf                                            |                                                                         |
| all/hap1/hap2 |                                                 |                                                                         |
|               | percent.accessible.bed.gz                       | % of fibers that are accessible on the haplotype                        |
| hap1-vs-hap2  |                                                 |                                                                         |
|               | asdf                                            |                                                                         |
|               | asdf                                            |                                                                         |
| trackHub      |                                                 |                                                                         |
|               | \*                                              | a trackHub directory that can be loaded into the UCSC browser           |
| .             |                                                 |                                                                         |
|               | GM12878_FDR.peaks-vs-percent.pdf                | Number of peaks vs % accessible                                         |

---

## Test data

You can find input data to test against at [this url](https://s3-us-west-2.amazonaws.com/stergachis-public1/index.html?prefix=Projects/Phased-GM12878/fire-test/).
