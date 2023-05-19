#chr1    28867   29204   chr1    28966   28977   126     7.000000        1       3       chr1    29152   29154   37      2.000000        0       1       UDN318336

OUT=./Tables/raw-hap1-vs-hap2-accessability.tsv
printf "ct\tst\ten\t" > $OUT
printf "hap1_ct\thap1_st\thap1_en\t" >> $OUT
printf "hap1_fdr\thap1_acc\thap1_link\thap1_nuc\t" >> $OUT
printf "hap2_ct\thap2_st\thap2_en\t" >> $OUT
printf "hap2_fdr\thap2_acc\thap2_link\thap2_nuc\t" >> $OUT
printf "imprinted\tsample\tcov\n" >> $OUT


for SM in UDN318336 GM12878_pacbiome HG002_pacbiome GM12878_130X; do
  COV=$(cat results/$SM/average.coverage.txt )
  paste \
    <(bedmap --delim '\t' --echo --max-element  <(cut -f 1-3 results/$SM/*peaks*.bed) <(zcat results/$SM/hap1/fdr.peaks.and.coverages.bed.gz)) \
    <(bedmap --max-element  <(cut -f 1-3 results/$SM/*peaks*.bed) <(zcat results/$SM/hap2/fdr.peaks.and.coverages.bed.gz)) \
    <(bedmap --indicator  <(cut -f 1-3 results/$SM/*peaks*.bed) <(sort-bed data/lcl_dmr_coordinates_Akbari.bed)) \
    | sed "s/$/\t$SM\t$COV/g" \
    >> $OUT

#  bedmap --range 50000 --indicator <(grep -v '^ct' Tables/raw-hap1-vs-hap2-accessability.tsv | sort-bed - ) <(sort-bed data/lcl_dmr_coordinates_Akbari.bed) > Tables/imprinted.tbl
done

head $OUT | column -t
tail $OUT | column -t

