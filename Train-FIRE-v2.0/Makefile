# shell
SHELL=/bin/bash -euo pipefail

# Variables
BAM := ../test/chr20.hifi.bam
BAM := ../../k-mer-variant-phasing/results/GM12878/GM12878.haplotagged.bam
SBAM := small.bam
REGIONS_BAM := regions.bam
ALL_POSITIVES := ./GM12878-all-accessible.bed.gz
MIXED_POSITIVES := ./GM12878-mixed-positive-accessible.bed.gz
ALL_NEGATIVES := ./GM12878-not-accessible.bed.gz
SHUF_NEGATIVES := ./GM12878-shuf-not-accessible.bed.gz
HG38 := ~/assemblies/hg38.analysisSet.fa.fai
N_SITES := 100000
SAMPLE_RATE := 0.10
TRAININ_DATA := ./training-data.bed.gz
FEATS := ./tmp.feats.gz
EXCLUDE := tmp.exclude.bed.gz
FT := ./fibertools-rs/target/release/ft
# Targets and their dependencies
all: $(FEATS) $(TRAININ_DATA)


$(EXCLUDE): ../workflow/annotations/hg38.gap.bed.gz ../workflow/annotations/SDs.merged.hg38.bed.gz 
	( \
		less $^; \
		cat ~/assemblies/hg38.analysisSet.fa.fai  | rg "_|chrEBV|chrMT|chrX|chrY|chrM" | awk -v OFS='\t' '{print $$1,0,$$2}' \
	) \
	| cut -f 1-3 | bedtools sort | bedtools merge | bgzip > $@

$(MIXED_POSITIVES): ./mixed-positives/ATAC.bed.gz ./mixed-positives/GM12878_DHS.bed.gz $(EXCLUDE)
	( \
		zcat ./mixed-positives/ATAC.bed.gz | awk '$$5 >= 1000'; \
		zcat ./mixed-positives/peaks_ENCFF073ORT.bed.gz; \
		zcat ./mixed-positives/peaks_ENCFF598KWZ.bed.gz; \
		zcat ./mixed-positives/peaks_CTCF_ENCFF951PEM.bed.gz \
	) \
		| cut -f 1-3 \
		| bedtools sort -g $(HG38) \
		| rg -v _ \
		| bedtools merge -d 147 \
		| (shuf | head -n $(N_SITES) || true) \
		| bedtools sort -g $(HG38) \
		| bedtools subtract -a - -b $(EXCLUDE) \
		| bgzip > $@ 
	printf "Number of mixed positives: "
	rb bl -r $@ 
	echo ""

$(ALL_POSITIVES): ./mixed-positives/*.bed*
	less $^ | cut -f 1-3 | bedtools sort -g $(HG38) | bedtools merge -d 147 | rg -v _ | bgzip > $@
	printf "Number of all positives: "
	rb bl -r $@ 
	echo ""

$(ALL_NEGATIVES): $(ALL_POSITIVES) $(EXCLUDE)
	bedtools complement -i $(ALL_POSITIVES) -g $(HG38) \
		| rg -v _ \
		| bedtools subtract -a - -b $(EXCLUDE) \
		| bgzip > $(ALL_NEGATIVES)
	printf "Number of all negatives: "
	rb bl -r $@ 
	echo ""


$(SHUF_NEGATIVES): $(ALL_POSITIVES) $(ALL_NEGATIVES) $(MIXED_POSITIVES)
	bedtools shuffle -excl $(ALL_POSITIVES) -incl $(ALL_NEGATIVES) -i $(MIXED_POSITIVES) -chrom -seed 42 -g $(HG38) | sort -k1,1 -k2,2n | bgzip > $@
	printf "Number of shuf negatives: "
	rb bl -r $@
	echo ""


$(SBAM): $(BAM)
	samtools view -@ `nproc` -b -s $(SAMPLE_RATE) $(BAM) -o $@ --write-index 

$(REGIONS_BAM): $(SBAM) $(SHUF_NEGATIVES) $(MIXED_POSITIVES) 
	samtools view -@ `nproc` -b -M \
		-L <(less $(SHUF_NEGATIVES) $(MIXED_POSITIVES)  ) \
		$(SBAM) -o $@

$(FEATS): $(REGIONS_BAM)
	pushd fibertools-rs/ && git pull && cargo build --release && popd
	$(FT) fire -t `nproc` -m 85 -f $^ \
		| awk 'NR == 1 || ($$3>$$2 && $$3-$$2 < 10000)' \
		| bgzip -@ 8 > $(FEATS)
		
		#| bedtools sort -header \

$(TRAININ_DATA): $(FEATS) $(MIXED_POSITIVES) $(ALL_NEGATIVES) $(SHUF_NEGATIVES) $(ALL_POSITIVES)
	(zcat $(FEATS) | head -n 1 || true) | sed 's/$$/\tLabel/' > training-data.bed
	bedtools intersect -f 0.25 -u -a $(FEATS) -b $(MIXED_POSITIVES) \
		| sed 's/$$/\t1/' \
		| rg '^chr' \
		>> training-data.bed
	bedtools intersect -f 0.25 -u -a $(FEATS) -b $(SHUF_NEGATIVES) \
		| bedtools intersect -v -a - -b $(ALL_POSITIVES) \
		| sed 's/$$/\t-1/' \
		| rg '^chr' \
		>> training-data.bed
	hck -F Label training-data.bed | uniq -c 
	bgzip -f -@ 16 training-data.bed

FIRE.FDR.pdf: $(TRAININ_DATA)
	python train-fire-model.py training-data.bed.gz --threads `nproc` --train-fdr 0.05 --direction msp_len_times_m6a_fc


.PHONY: clean
clean:
	rm -f tmp.pos.bed tmp.neg.bed tmp.bed tmp.bam training-data.bed.gz $(EXCLUDE)  $(ALL_POSITIVES) $(MIXED_POSITIVES) $(REGIONS_BAM)



