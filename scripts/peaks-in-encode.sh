for f in results/*{GM20129_PS00447,HG02630_PS00445,UDN318336,_pacbiome}*/FIRE.peaks.*bed ; do 
  echo $f
  a=$(bedtools intersect -u -a $f -b ./data/ENCODE3_consensus_DHS_ENCFF503GCK.tsv.gz -sorted | wc -l)
  b=$(grep -v '^#' $f | wc -l)
  echo "$a   $b"
  echo "$a   $b" | awk '{print 100*$1/$2}'
  echo ""
done

exit 
