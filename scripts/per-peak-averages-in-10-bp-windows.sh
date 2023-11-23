PEAKS="results/GM12878_FDR/FDR-peaks/FDR-FIRE-peaks.bed.gz"
DATA1="results/GM12878_FDR/all/percent.accessible.bed.gz"
DATA2="data/bedgraph_annotations/ENCFF960FMM_dnase_signal.bed"
DATA3="data/ATAC/10X_GM12878_aggr_scATAC.bg.gz"
FAI="/mmfs1/home/mvollger/assemblies/hg38.analysisSet.fa.fai"
S=2000
W=25

bioawk -tc hdr '{print $1,$2,$3,$1"_"$2"_"$3,$score,$FDR}' $PEAKS \
  | bedtools slop -b $S -g $FAI -i - \
  | bedtools makewindows -w $W -b - -i src \
  | sort-bed - \
  | bedmap --ec --echo --max --delim '\t' - <(zcat -f $DATA1 | bioawk -t '{print $0,$4}' ) \
  | bedmap --ec --echo --max --delim '\t' - <(bioawk -t '{print $0,$4}' $DATA2) \
  | bedmap --ec --echo --max --delim '\t' - <(zcat $DATA3 | bioawk -t '{print $0,$4}') \
  | bgzip -@ 8 \
> Tables/per_site_windows_dnase_and_percent_fire.tbl.gz



