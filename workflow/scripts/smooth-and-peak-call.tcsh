#!/bin/tcsh -ef
# author : sjn

if ( $#argv != 5 ) then
  printf "%s\n" $0
  printf "  <Chromosome>\n"
  printf "  <BED>\n"
  printf "  <FDR-Threshold>\n"
  printf "  <Output>\n"
  printf "  <details>\n"
  exit -1
endif

set echo # now echo all lines

set chrom = $1
set bed_sorted = $2 # a BED file sorted per sort-bed of BEDOPS
set fdr_thold = $3
set output = $4
set details = $5

# settings determined by Andrew
set waveletlvl=7
set filter_type=Haar
set boundary_type=reflected

set tmpd = $TMPDIR/`whoami`/$$
if ( -d $tmpd ) then
  rm -rf $tmpd
endif
mkdir -p $tmpd

# create values per base
tabix $bed_sorted $chrom \
  | awk 'BEGIN {s=$2} ; { while(s<$3) { print $4; s+=1 } }' \
  | tee $tmpd/$output:t.$chrom.perbase \
  | modwt --operation smooth --level $waveletlvl --to-stdout --boundary $boundary_type --filter $filter_type - \
 >! $tmpd/$output:t.$chrom.waves

# $0:h/per-chrom-peakcall-details.tcsh 
$details \
  $tmpd \
  $chrom \
  $tmpd/$output:t.$chrom.perbase \
  $tmpd/$output:t.$chrom.waves \
  $fdr_thold \
  $output

rm -rf $tmpd

exit 0
