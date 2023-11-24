# shell
SHELL=/bin/bash -euo pipefail

# Variables
BAM := ../test/chr20.hifi.bam
BAM := ../../k-mer-variant-phasing/results/GM12878/GM12878.haplotagged.bam
SBAM := small.bam
ALL_POSITIVES := ./GM12878-all-accessible.bed.gz
MIXED_POSITIVES := ./GM12878-mixed-positive-accessible.bed.gz
ALL_NEGATIVES := ./GM12878-not-accessible.bed.gz
SHUF_NEGATIVES := ./GM12878-shuf-not-accessible.bed.gz
HG38 := ~/assemblies/hg38.analysisSet.fa.fai
N_SITES := 10000
TRAININ_DATA := ./training-data.bed.gz
FEATS := ./tmp.feats.gz
EXCLUDE := tmp.exclude.bed.gz
# Targets and their dependencies
all: $(FEATS) $(TRAININ_DATA)


$(EXCLUDE): ../workflow/annotations/hg38.gap.bed.gz ../workflow/annotations/SDs.merged.hg38.bed.gz 
	less $^ | cut -f 1-3 | bedtools sort | bedtools merge | bgzip > $@

$(MIXED_POSITIVES): ./mixed-positives/ATAC.bed.gz ./mixed-positives/GM12878_DHS.bed.gz $(EXCLUDE)
	bedtools intersect -a $(word 1,$^) -b $(word 2,$^) \
		| cut -f 1-3 \
		| (shuf | head -n $(N_SITES) || true) \
		| bedtools sort -g $(HG38) \
		| rg -v _ \
		| bedtools merge -d 147 \
		| bedtools subtract -a - -b $(EXCLUDE) \
		| bgzip > $@ 

$(ALL_POSITIVES): ./mixed-positives/*.bed*
	less $^ | cut -f 1-3 | bedtools sort -g $(HG38) | bedtools merge -d 147 | rg -v _ | bgzip > $@

$(ALL_NEGATIVES): $(ALL_POSITIVES) $(EXCLUDE)
	bedtools complement -i $(ALL_POSITIVES) -g $(HG38) \
		| rg -v _ \
		| bedtools subtract -a - -b $(EXCLUDE) \
		| bgzip > $(ALL_NEGATIVES)

$(SHUF_NEGATIVES): $(ALL_POSITIVES) $(ALL_NEGATIVES) $(MIXED_POSITIVES)
	bedtools shuffle -excl $(ALL_POSITIVES) -incl $(ALL_NEGATIVES) -i $(MIXED_POSITIVES) -chrom -seed 42 -g $(HG38) | sort -k1,1 -k2,2n > $@

$(SBAM): $(BAM) $(SHUF_NEGATIVES) $(MIXED_POSITIVES) 
	samtools view -@ `nproc` -b -M \
		-L <(less $(SHUF_NEGATIVES) $(MIXED_POSITIVES)  ) \
		$(BAM) -o $@

$(FEATS): $(SBAM)
	ft fire -t `nproc` -m 0 -f $(SBAM) \
		| awk 'NR == 1 || ($$3>$$2 && $$3-$$2 < 10000)' \
		| bgzip -@ 8 > $(FEATS)
		
		#| bedtools sort -header \

$(TRAININ_DATA): $(FEATS) $(MIXED_POSITIVES) $(ALL_NEGATIVES) $(SHUF_NEGATIVES) $(ALL_POSITIVES)
	(zcat $(FEATS) | head -n 1 || true) | sed 's/$$/\tLabel/' > training-data.bed
	bedtools intersect -f 0.25 -u -a $(FEATS) -b $(MIXED_POSITIVES) | sed 's/$$/\t1/' >> training-data.bed
	bedtools intersect -f 1.0 -u -a $(FEATS) -b $(SHUF_NEGATIVES) \
		| bedtools intersect -v -a - -b $(ALL_POSITIVES) \
		| sed 's/$$/\t-1/' >> training-data.bed
	hck -F Label training-data.bed | uniq -c 
	bgzip -f -@ 16 training-data.bed

FIRE.FDR.pdf: $(TRAININ_DATA)
	python train-fire-model.py training-data.bed.gz --threads `nproc` --train-fdr 0.01


.PHONY: clean
clean:
	rm -f tmp.pos.bed tmp.neg.bed tmp.bed tmp.bam training-data.bed.gz $(EXCLUDE)  $(ALL_POSITIVES) $(MIXED_POSITIVES) 


