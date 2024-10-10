#!/bin/bash
# author : sjn
# date : July.2023

set -euox pipefail

if [[ $# != 3 ]]; then
    printf "Expect $0 <clustering-vs-null.bed.gz> <fire.peaks.and.coverages.bed.gz> <outfile>\n"
    exit 1
fi

oe=$(zcat $1 |
    awk '{ \
             coverage[$4][$5]+=$3-$2; \
           } END { \
             real_bp = 0; \
             over_expected = 0; \
             for ( cov in coverage ) { \
               real_bp += cov * coverage[cov]["Real"]; \
               val = cov * (coverage[cov]["Real"] - coverage[cov]["Null"]); \
               if ( val > 0 ) { over_expected += val; } \
             } \
             print 100*over_expected/real_bp; \
          }')
n_tests=$(zcat $2 | wc -l)
min_fdr=$(echo "-10.0*(l(0.01/${n_tests})/l(10.0))" | bc -lq)
# n_peaks=$(zcat $2 | awk -v m="${min_fdr}" '$4 >= m' | wc -l)

printf "percent-of-MSPs-preferentially-clustered-along-the-genome\tmin_fdr\n" >$3
printf "%s%%\t%s\n" ${oe} ${min_fdr} >>$3

exit 0
