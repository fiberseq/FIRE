b=results/GM12878_FDR/all/percent.accessible.bed.gz

zcat data/ENCFF762CRQ_DNase_peaks.bed.gz | cut -f 1-3,7 \
  | bedmap --ec --echo --max --delim '\t' - <(zcat $b | bioawk -t '{print $1,$2,$3,".",$4}') \
  | bgzip \
  > ./Tables/DNase-peaks-vs-percent-fire.bed.gz
  



