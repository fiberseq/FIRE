bedtools closest \
  -a <(zcat ./transcripts.old.data.bed.gz | rg -v "^chrom" | bedtools sort ) \
  -b ~/projects/fire-figures/data/gencode.v42.annotation_TSS.gff3.gz \
  -f 0.1 \
  | grep -vw "\-1" \
  | sed 's/\tchr.*gene_name=/\t/g' | sed 's/;.*//g' \
  | bgzip -@ 16 \
  > ./transcripts.old.data.genes.bed.gz
