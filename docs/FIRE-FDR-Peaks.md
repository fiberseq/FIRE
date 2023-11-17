# Description of the FIRE peaks

| Column           | Description                                                       |
| ---------------- | ----------------------------------------------------------------- |
| chrom            | Chromosome                                                        |
| peak_start       | Start of the peak                                                 |
| peak_end         | End of the peak                                                   |
| start            | Start of the peak                                                 |
| end              | End of the peak                                                   |
| fire_coverage    | # of fire elements overlapping the peak                           |
| coverage         | Coverage of the peak                                              |
| score            | Max FIRE-score across the peak                                    |
| FDR              | FDR of the peak                                                   |
| log_FDR          | -10\*log10(FDR) of the peak                                       |
| fire_coverage_H1 | # of fire elements overlapping the peak in H1                     |
| coverage_H1      | Coverage of the peak in H1                                        |
| score_H1         | Max FIRE-score across the peak in H1                              |
| FDR_H1           | FDR of the peak in H1                                             |
| log_FDR_H1       | -10\*log10(FDR) of the peak in H1                                 |
| fire_coverage_H2 | # of fire elements overlapping the peak in H2                     |
| coverage_H2      | Coverage of the peak in H2                                        |
| score_H2         | Max FIRE-score across the peak in H2                              |
| FDR_H2           | FDR of the peak in H2                                             |
| log_FDR_H2       | -10\*log10(FDR) of the peak in H2                                 |
| max_window_score | Max FIRE-score across the peak                                    |
| FIRE_size_mean   | Mean size of FIRE elements in the peak                            |
| FIRE_size_ssd    | Standard deviation of the size of FIRE elements in the peak       |
| FIRE_start_ssd   | Standard deviation of start position of FIRE elements in the peak |
| FIRE_end_ssd     | Standard deviation of end position of FIRE elements in the peak   |
| local_max_count  | Number of local maxima in the peak                                |
| peak_length      | Length of the peak                                                |


## FIRE score calculation

The FIRE score ($S$) for a position in the genome ($g$) is calculated using the following formula:

$$ S_g = -\frac{50}{C_g} \sum_{i=1}^{C_g} log_{10}(p_i) $$

Where $C_g$ is the Fiber-seq coverage overlapping position $g$, and $p_i$ is $1$ minus the precision of the FIRE element, i.e., the chance that the FIRE element is a false positive. Note, the highest confidence a FIRE element can have is is $0.99$, i.e., $s_i$ is at least $0.01$.  Therefore, the FIRE score ($S_g$) can take on values between $0$ (none of the overlapping Fiber-seq reads have FIRE elements) and $100$ (all of the overlapping Fiber-seq reads have FIRE elements with a precision of $0.99$).

Regions covered by less than $3$ FIRE elements are not scored and given a value of $-1$.

## FDR calculation
FDR is calculated by shuffling the locations of all the fibers across the genome and recalculating the FIRE score for each position. The FDR is then calculated as the fraction of bases that have shuffled FIRE scores that are higher than the FIRE score of the un-shuffled data at various FIRE score thresholds. The shuffle and FDR calculations are performed per sample and the results are recorded in the file `FDR-peaks/FIRE.score.to.FDR.tbl` with the following columns:

| Column | description |
| ------ | ----------- |
| threshold | FIRE score threshold |
| FDR | FDR at threshold |
| shuffled_peaks | Number of base-pairs at threshold in shuffled data |
| peaks | Number of base-pairs at threshold in un-shuffled data |

## Peak calling
TODO 

## Determination of peak start and end
TODO