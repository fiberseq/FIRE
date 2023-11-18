# FIRE: descriptions, methods, and outputs 

For details on running the FIRE pipeline see the [README.md](/README.md).

## Fiber-seq inferred regulatory elements (FIREs)

FIREs are MTase sensitive patches (MSPs) that are inferred to be regulatory elements on single chromatin fibers. To do this we used semi-supervised machine learning to identify MSPs that are likely to be regulatory elements using the `Mokapot` framework and `XGBoost`. Every individual FIRE element is associated with a precision value, which indicates the probability that the FIRE element is a true regulatory element. The precision of FIREs elements are estimated using `Mokapot` and validation data not used in training. We train our model targeting FIRE elements with at least 90% precision, MSPs with less than 90% precision are considered to have average level of accessibility expected between two nucleosomes, and are referred to as linker regions.

Semi-superivized machine learning with `Mokapot` requires a mixed-positive training set and a clean negative training set. To create mixed positive training data we selected MSPs that overlapped DNase hypersensitive sites (DHSs) and CTCF ChIP-seq peaks. And to create a clean negative training set we selected MSPs that did not overlap DHSs or CTCF ChIP-seq peaks. 

The following features were used as features with `Mokapot` in FIRE element classification:

| Feature | Description |
| ------- | ----------- |
| msp_len | Length of the MSP |
| msp_len_times_m6a_fc | Length of the MSP times the fold-change of m6A in the MSP |
| ccs_passes | Number of CCS passes in the MSP |
| fiber_m6a_count | Number of m6A sites in the Fiber-seq read |
| fiber_AT_count | Number of AT sites in the Fiber-seq read |
| fiber_m6a_frac | Fraction of m6A sites in the Fiber-seq read |
| msp_m6a | Number of m6A sites in the MSP |
| msp_AT | Number of AT sites in the MSP |
| msp_m6a_frac | Fraction of m6A sites in the MSP |
| msp_fc | Fold-change of m6A fraction in the MSP relative to the whole reads m6A fraction |
| m6a_count_X | Number of m6A sites in the Xth 40 bp window |
| AT_count_X | Number of AT sites in the Xth 40 bp window |
| m6a_frac_X | Fraction of m6A sites in the Xth 40 bp window |
| m6a_fc_X | Fold-change of m6A fraction in the Xth 40 bp window relative to the whole reads m6A fraction |

There are nine 40 bp windows for each MSP (X = 1, 2, ..., 9), with the 5th window centered on the MSP. 

## The FIRE score track

The FIRE score track is a genome-wide track that reflects the aggregation of FIRE elements for all the Fiber-seq reads overlapping a genomic position. The FIRE score track for Fiber-seq is comparable to DNase-seq and ATAC-seq accessibility tracks. Details on the FIRE score calculation are provided in the next section.

### FIRE score calculation

The FIRE score ($S$) for a position in the genome ($g$) is calculated using the following formula:

$$ S_g = -\frac{50}{C_g} \sum_{i=1}^{C_g} \log_{10}(p_i) $$

Where $C_g$ is the Fiber-seq coverage overlapping position $g$, and $p_i$ is $1$ minus the precision of the FIRE element, i.e., the chance that the FIRE element is a false positive. Note, the highest confidence a FIRE element can have is is $0.99$, i.e., $s_i$ is at least $0.01$.  Therefore, the FIRE score ($S_g$) can take on values between $0$ (none of the overlapping Fiber-seq reads have FIRE elements) and $100$ (all of the overlapping Fiber-seq reads have FIRE elements with a precision of $0.99$).

Regions covered by less than $3$ FIRE elements are not scored and given a value of $-1$.

## FDR track calculation
FDR calculation begins by shuffling the locations of all the fibers across the genome and recalculating the FIRE score for each position in the genome. The FDR is then defined as the number of bases that have shuffled FIRE scores above a threshold divided by the number of bases in the un-shuffled data. The shuffle and FDR calculations are performed per sample and the results are recorded in the file `FDR-peaks/FIRE.score.to.FDR.tbl` with the following columns:

| Column | description |
| ------ | ----------- |
| threshold | FIRE score threshold |
| FDR | FDR at threshold |
| shuffled_peaks | Number of base-pairs at threshold in shuffled data |
| peaks | Number of base-pairs at threshold in un-shuffled data |

## Peak calling

Peaks are called by identifying FIRE score local-maxima that have FDR values below a threshold. By default the pipeline reports results for both 1% and 5% FDR thresholds. Once a local-maxima is identified, the start and end positions of the peak are determined by the median start and end positions of the underlying FIRE elements. We also calculate and report wide peaks by taking the union of the FIRE peaks and all regions below the FDR threshold and then merging resulting regions that are within one nucleosome (147 bp) of one another.

The peak results are reported in the file `FDR-peaks/FDR-FIRE-peaks.bed.gz` with the following columns:

| Column           | Description                                                       |
| ---------------- | ----------------------------------------------------------------- |
| chrom            | Chromosome                                                        |
| peak_start       | Start of the peak                                                 |
| peak_end         | End of the peak                                                   |
| start            | Start of the local-maxima of the peak                             |
| end              | End of the local-maxima the peak                                  |
| fire_coverage    | # of fire elements overlapping the peak                           |
| coverage         | Coverage of the peak                                              |
| score            | Max FIRE-score across the peak                                    |
| FDR              | FDR of the peak                                                   |
| log_FDR          | $-10 \log_{10}(FDR)$ of the peak                                  |
| fire_coverage_H1 | # of fire elements overlapping the peak in H1                     |
| coverage_H1      | Coverage of the peak in H1                                        |
| score_H1         | Max FIRE-score across the peak in H1                              |
| FDR_H1           | FDR of the peak in H1                                             |
| log_FDR_H1       | $-10 \log_{10}(FDR)$  of the peak in H1                           |
| fire_coverage_H2 | # of fire elements overlapping the peak in H2                     |
| coverage_H2      | Coverage of the peak in H2                                        |
| score_H2         | Max FIRE-score across the peak in H2                              |
| FDR_H2           | FDR of the peak in H2                                             |
| log_FDR_H2       | $-10 \log_{10}(FDR)$  of the peak in H2                           |
| max_window_score | Max FIRE-score across the peak                                    |
| FIRE_size_mean   | Mean size of FIRE elements in the peak                            |
| FIRE_size_ssd    | Standard deviation of the size of FIRE elements in the peak       |
| FIRE_start_ssd   | Standard deviation of start position of FIRE elements in the peak |
| FIRE_end_ssd     | Standard deviation of end position of FIRE elements in the peak   |
| local_max_count  | Number of local maxima in the peak                                |
| peak_length      | Length of the peak                                                |

## General description of the outputs of the FIRE pipeline

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
|               | FDR-FIRE-peaks.bed.gz                           | Fiber-seq peak calls                                                    |
|               | FDR-wide-peaks.bed.gz                           | Fiber-seq wide peak calls                                               |
|               | FDR.track.bed.gz                                | Track of FDR significance of accessibility                              |
|               | {sm}.peaks-vs-percent.pdf                       | Number of peaks vs % accessible                                         |
| all/hap1/hap2 |                                                 |                                                                         |
|               | percent.accessible.bed.gz                       | % of (haplotype) fibers that are accessible                             |
| hap1-vs-hap2  |                                                 |                                                                         |
|               | FIRE.hap.differences.bed                        | Large table of FIREs that are different between hap1 and hap2           |
|               | hap1-vs-hap2-volcano.pdf                        | Volcano plot of FIREs that are different between hap1 and hap2          |
|               | hap1-vs-hap2.pdf                                | Scatter plot of FIREs and their percent accessibility for each hap      |
| trackHub      |                                                 |                                                                         |
|               | \*                                              | a trackHub directory that can be loaded into the UCSC browser           |

---