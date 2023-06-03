#ls results/*{UDN318336,_pacbiome}/FIRE.bed.gz | parallel -n 1 $'bgzip -@ 16 -cd {} | awk \'$5 <= 5\' | bgzip -@ 16 > {.}.only_acc.bed.gz'
#ls results/*{UDN318336,_pacbiome}/all/acc.model.results.bed.gz | parallel -n 1 $'printf "{}\t"; bgzip -@ 16 -cd {} | wc -l'
ls results/*{UDN318336,_pacbiome}/FIRE.bed.gz | parallel -n 1 $'printf "{.}\t"; bedtools intersect -u -a {} -b ./data/ENCODE3_consensus_DHS_ENCFF503GCK.tsv -sorted | wc -l'
ls results/*{UDN318336,_pacbiome}/FIRE.bed.gz | parallel -n 1 $'printf "{.}\t"; bgzip -cd -@16 {} | wc -l'
exit 
ls results/*{UDN318336,_pacbiome}/all/acc*.only_acc.bed.gz | parallel -n 1 $'bedtools merge -i {} -o count -c 1 | bgzip -@16 > {.}.merged.count.bed.gz'
ls results/*{UDN318336,_pacbiome}/all/acc*.only_acc.bed.gz | parallel -n 1 $'bgzip -cd -@16 {} | awk \'$5<=5\' | bedtools merge -i - -o count -c 1 | bgzip -@16 > {.}.merged.count.bed.gz'

