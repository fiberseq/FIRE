# Description of the FIRE peaks

TODO in progress

#chrom peak_start peak_end start end fire_coverage coverage score FDR log_FDR fire_coverage_H1 coverage_H1 score_H1 FDR_H1 log_FDR_H1 fire_coverage_H2 coverage_H2 score_H2 FDR_H2 log_FDR_H2 max_window_score FIRE_size_mean FIRE_size_ssd FIRE_start_ssd FIRE_end_ssd local_max_count peak_length

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
