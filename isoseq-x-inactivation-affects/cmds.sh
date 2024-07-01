bedtools makewindows -w 50000 -g ~/assemblies/hg38.analysisSet.fa.fai  > windows.bed


bedtools intersect -sorted -c \
  -a windows.bed \
  -b H1.bam \
  | bedtools intersect -sorted -c \
  -a - \
  -b H2.bam \
  > windows.with.counts.bed

