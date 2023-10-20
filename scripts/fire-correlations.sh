a=results/GM12878_FDR/FDR-peaks/FDR-FIRE-peaks.bed.gz
b=results/PS00356_COLO829BL_2/FDR-peaks/FDR-FIRE-peaks.bed.gz
c=results/PS00388_COLO829BL_1/FDR-peaks/FDR-FIRE-peaks.bed.gz

printf "#chrom\tstart\tend\tGM12878\tGM12878_30X\tPS00388\tPS00356\tDNase\n" >"./Tables/fire-correlations.tbl"

bedops --ec -m <(zcat $a) <(zcat $b) <(zcat $c) |
    bedmap --ec --echo --max --delim '\t' - <(hck -z -f 1-4 -F score results/GM12878_FDR/FDR-peaks/FDR.track.bed.gz) |
    bedmap --ec --echo --max --delim '\t' - <(hck -z -f 1-4 -F score results/GM12878_FDR_30X/FDR-peaks/FDR.track.bed.gz) |
    bedmap --ec --echo --max --delim '\t' - <(hck -z -f 1-4 -F score results/PS00388_COLO829BL_1/FDR-peaks/FDR.track.bed.gz) |
    bedmap --ec --echo --max --delim '\t' - <(hck -z -f 1-4 -F score results/PS00356_COLO829BL_2/FDR-peaks/FDR.track.bed.gz) |
    bedmap --ec --echo --max --delim '\t' - <(bioawk -t '{print $1,$2,$3,".",$4}' data/bedgraph_annotations/ENCFF960FMM_dnase_signal.bed) \
        >>"./Tables/fire-correlations.tbl"
head ./Tables/fire-correlations.tbl
bgzip -f Tables/fire-correlations.tbl
