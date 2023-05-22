#!/bin/tcsh -ef
# author : sjn

if ( $#argv != 6 ) then
  printf "%s\n" $0
  printf "  <working-dir>\n"
  printf "  <Chromosome>\n"
  printf "  <File of Chromosome's per-base FDR values>\n"
  printf "  <File of Chromosome's smoothed per-base FDR values>\n"
  printf "  <FDR-thold>\n"
  printf "  <Output>\n"
  exit -1
endif

set echo # now echo all lines

set tmpd = $1     # cleaned up by the calling script
set chrom = $2
set perbase = $3  # raw values
set smoothed = $4 # smoothed values
set min_fdr = $5
set output = $6

# Andrew's settings and directions:
#   to merge or not-to-merge
#    Take 2 neighboring peaks, if intervening trough's height (raw) is X% ($merge_trough_check) of height (or less)
#      keep the peaks separate.  Otherwise, merge them.
set merge_trough_check = 0.5 # (times peak's raw height)
set peak_border_height_fraction = 0.25 # (times peak's raw height)
set buffer = 500

set lraw = $tmpd/.pb.$chrom

mkfifo $lraw.1
mkfifo $lraw.2
mkfifo $lraw.3

(awk -v chr=$chrom 'BEGIN {OFS="\t"} ; { print chr, NR-1, NR, "ID-"NR, $1 }' $perbase \
 >! $lraw.1) &

(awk -v chr=$chrom 'BEGIN {OFS="\t"} ; { print chr, NR-1, NR, "ID-"NR, $1 }' $perbase \
 >! $lraw.2) &

(awk -v chr=$chrom 'BEGIN {OFS="\t"} ; { print chr, NR-1, NR, "ID-"NR, $1 }' $perbase \
 >! $lraw.3) &

paste $smoothed $perbase \
  | awk \
    'BEGIN {OFS="\t"} ; { \
      if ( $1 < 0 ) { $1=0; } \
      if ( NR == 1 ) { \
        prev=$1; \
      } else { \
        print $1-prev, $1, $2; \
        prev = $1; \
      } \
    }' \
  | awk -v chr=$chrom \
    'BEGIN {OFS="\t";start=0;} ; { \
      if ( 0 == start ) { \
        if ( $1 != 0 ) { \
          cu = ($1 > 0); \
          #print chr, NR-1, NR, "v", smooths, raws; \
          smooths = $2; \
          raws = $3; \
          start = 1; \
        } \
      } else if ( $1 != 0 ) { \
        ncu = ($1 > 0); \
        if ( cu != ncu ) { \
          if ( 0 != ncu ) { \
            #print chr, NR-1, NR, "v", smooths, raws; \
          } else { \
            if ( raws > 0 ) { \
              print chr, NR-1, NR, "p", smooths, raws; \
            } \
          } \
          cu = ncu; \
        } else if ( $2 == 0 ) { \
          #print chr, NR-1, NR, "v", smooths, raws; \
          cu = 0; \
        } \
        smooths = $2; \
        raws = $3; \
      } \
    }' \
  | bedmap --sweep-all --chrom $chrom --prec 0 --range $buffer --echo --echo-map-score - $lraw.1 \
  | awk -v thold=$peak_border_height_fraction -v merge=$merge_trough_check -v buff=$buffer \
    'BEGIN {FS="|"; OFS="\t"; if(merge<thold) {exit -1}} ; { \
      split($1, a, "\t"); \
      goal=thold*a[6]; \
      split_goal=merge*a[6]; \
      n=split($2, scores, ";"); \
      summit_loc=a[2]; \
      start_loc=summit_loc-buff; \
      if ( start_loc < 0 ) { next; } \
      idx=buff; \
      cntr = 1; \
      split_left = -1; \
      found1 = 0; \
      found2 = 0; \
      while ( idx > 0 ) { \
        val=scores[idx]; \
        if ( val <= split_goal ) { \
          if ( split_left < 0 ) {found2=1; split_left=summit_loc-cntr+1;} \
        } \
        if ( val <= goal ) { \
          start_loc=summit_loc-cntr+1; \
          found1 = 1; \
          break; \
        } \
        idx -= 1; \
        cntr += 1; \
      } \
      if ( found2 == 0 ) { split_left=summit_loc-cntr+1;} \
      if ( found1 == 0 ) { start_loc=split_left;} \
      stop_loc=a[3]; \
      idx=buff+2; \
      end_idx=buff+2+buff; \
      start_idx = idx; \
      split_right = -1; \
      found1 = 0; \
      found2 = 0; \
      while ( idx < end_idx ) { \
        val=scores[idx]; \
        if ( val <= split_goal ) { \
          if ( split_right < 0 ) {found2=1; split_right=summit_loc+(idx-start_idx);} \
        } \
        if ( val <= goal ) { \
          stop_loc=summit_loc+(idx-start_idx); \
          found1 = 1; \
          break; \
        } \
        ++idx; \
      } \
      if ( found2 == 0 ) { split_right=summit_loc+(idx-start_idx); } \
      if ( found1 == 0 ) { stop_loc=split_right; } \
      print a[1], start_loc, stop_loc+1, "ID-"NR, a[5], a[6], split_left, split_right+1, a[2]; \
    }' \
  | awk \
    'BEGIN {OFS="\t"; min_start=0} ; { \
      if ( NR > 1 ) { \
        if ( $3 < min_start ) { next } \
        if ( $2 < min_start ) { $2 = min_start; } \
        if ( $7 < min_start ) { $7 = min_start; } \
        split(prev_peak, pp, "\t"); \
        if ( pp[3] > $2 ) { \
          if ( pp[8] < $9 ) { # split \
            pp[3] = pp[8]; \
            $2 = pp[8]; \
            min_start = pp[3]; \
            check = 1; \
            if ( $2 < $7 ) { $2 = $7; check=0; } \
            print pp[1], pp[2], pp[3], pp[4], pp[5], pp[6], pp[7], pp[8], pp[9]; \
            if ( check > 0 ) { $4 = "check-left"; } \
            prev_peak = $0; \
          } else if ( $7 > pp[9] ) { # split \
            pp[3] = $7-1; \
            $2 = $7; \
            min_start = pp[3]; \
            check = 1; \
            if ( pp[8] < pp[3] ) { pp[3] = pp[8]; check=0; } \
            print pp[1], pp[2], pp[3], pp[4], pp[5], pp[6], pp[7], pp[8], pp[9]; \
            if ( check > 0 ) { $4 = "check-right"; } \
            prev_peak = $0; \
          } else { # merge \
            chr = pp[1]; \
            start_coord = pp[2]; \
            end_coord = pp[3]; \
            id = pp[4]; \
            val1 = pp[5]; \
            val2 = pp[6]; \
            split_left = pp[7]; \
            split_right = pp[8]; \
            summit = pp[9]; \
            if ( $2 > start_coord ) { start_coord = $2; } \
            if ( $3 < end_coord ) { end_coord = $3; } \
            if ( $6 > val2 ) { val1 = $5; val2 = $6; summit = $9 } \
            if ( $7 > split_left ) { split_left = $7; } \
            if ( $8 < split_right ) { split_right = $8; } \
            if ( start_coord < min_start ) { start_coord = min_start } \
            prev_peak = chr"\t"start_coord"\t"end_coord"\t"id"\t"val1"\t"val2"\t"split_left"\t"split_right"\t"summit; \
          } \
        } else { \
          print prev_peak; \
          prev_peak = $0; \
        } \
      } else { \
        prev_peak = $0; \
      } \
    } END { print prev_peak }' \
  | bedmap --sweep-all --chrom $chrom --prec 0 --delim "\t" --echo --echo-map-score - $lraw.2 \
  | awk 'BEGIN {OFS="\t"} ; { \
         if ( $4 ~ /check/ ) { \
           n=split($NF, scores, ";"); \
           start=$2; \
           s_val=scores[1]; \
           end=$3; \
           e_val=scores[n]; \
           i=2; \
           if ( $4 == "check-left" ) { \
             for ( ; i<=$9-$2; ++i ) { \
               if ( scores[i] == s_val ) { start=$2+i; } \
               else { break; } \
             } \
           } else { \
             for ( i=n-1; i>$9-$2+1; --i ) { \
               if ( scores[i] == e_val ) { end=end-1; } \
               else { break; } \
             } \
           } \
           $2=start; $3=end; \
           print $0; \
         } else { print; } \
       }' \
  | bedmap --sweep-all --chrom $chrom --delim "\t" --echo --mean --sum --max - $lraw.3 \
  | awk 'BEGIN {OFS="\t"} ; { printf "%s\t%s\t%s\tID-%s\t%.0f\t%.0f\t%.0f\t%s\t%s\t%.0f\n", $1, $2, $3, NR, $11, $12, $13, $9, $6, $5 }' \
  | awk -v thold=$min_fdr '$7 >= thold' \
 >! $output

rm -f $lraw.1 $lraw.2 $lraw.3

exit 0
