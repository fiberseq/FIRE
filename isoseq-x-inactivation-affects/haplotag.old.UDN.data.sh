whatshap haplotag \
  --ignore-read-groups \
  ../../phased-bams/UDN318336_WGS/UDN318336_WGS.haplotagged.vcf.gz \
  ./UDN318336.isoseq.bam \
  -r ~/assemblies/hg38.analysisSet.fa \
  | samtools view \
    -b -@ 16 \
    -o mat.pat.old.data.tagged.bam --write-index

