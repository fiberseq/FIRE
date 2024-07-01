whatshap haplotag \
  --ignore-read-groups \
  ../../phased-bams/GM12878_WGS/GM12878_WGS.haplotagged.vcf.gz \
  ../../isoseq/GM12878.isoseq.bam \
  -r ~/assemblies/hg38.analysisSet.fa \
  | samtools view \
    -b -@ 16 \
    -o gm12878.haplotagged.bam --write-index

