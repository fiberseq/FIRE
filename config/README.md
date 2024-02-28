# Inputs and options for the configuration yaml file
see `config.yaml` for an example.

## Required input options
Reference name, this is the name that will be used as the genome name in the UCSC track hub, so be sure to use a valid UCSC genome name when possible. 
```
ref_name: hg38
```
Reference `fasta` file:
```
ref: /path/to/hg38.fa
```
Manifest of input sample(s), must have two white-space separated columns: sample name (`sample`) and input bam file path (`bam`). See `config.tbl` for an example.
```
manifest: config/config.tbl
```


## Optional input options
Max number of threads to use in very distributed steps:
```
max_t: 4
```

### Coverage options
Filter peaks that are more than `x` standard deviations from the mean coverage:
```
coverage_within_n_sd: 5
```

### FDR calling regions
Exclude the following regions when creating a null distribution of FIRE scores. This is not needed if your reference is `hg38` as the defaults will automatically load!
```
excludes:
    - annotations/gaps.bed
    - annotations/centromeres.bed
    - annotations/cnvs.bed
```

# Developer options, not for general use
Min coverage of FIRE elements for calling a FIRE peak:
```
min_coverage: 4
```
Forgo the use of FDR peak calling and instead call peaks for regions with at least this fraction of reads actuated.
```
min_per_acc_peak = 0.25
```
Apply a percent actuation filter on top of the FDR peak calling. Default is `0.0` for no filter.
```
min_frac_accessible: 0.0
```
Process only chromosomes matching this regular expression:
```
keep_chromosomes: "chr[0-9XY]+$"
```
The false discovery rate (FDR) for calling FIRE peaks:
```
max_peak_fdr: 0.05
```
The allowed false discovery rate for calling FIRE elements. 
```
min_fire_fdr: 0.10
```
