#
# Index the chain file for leviosam2.
#
# Input here is a chain file that defines the alignment between the DSA and the reference genome
# at a contig level (>100 kbp of alignment).
#
# The output is a special leviosam2 index file that is used to lift over the alignments from the DSA to the reference genome.
#
rule leviosam2_index:
    input:
        chain=DSA_CHAIN,
        fai=FAI,
    output:
        index=temp("temp/{sm}/leviosam2-index/index.clft"),
    conda:
        DEFAULT_ENV
    threads: 1
    resources:
        mem_mb=64 * 1024,
        runtime=16 * 60,
    shell:
        """
        {LEVIO_EXE} index \
            -p results/leviosam2-index/index \
            -c {input.chain} \
            -F {input.fai}
        """


#
# Lift over the alignments from the DSA to the reference genome using the chain file / leviosam2 index.
#
# This is not a realignment, but a lift over of the reads from the DSA to the reference genome.
#
rule leviosam2:
    input:
        bam=rules.fire.output.bam,
        levio_index=rules.leviosam2_index.output.index,
        ref=REF,
    output:
        lifted=temp("temp/{sm}/leviosam2/{sm}-{chrom}-committed.bam"),
        deferred=temp("temp/{sm}/leviosam2/{sm}-{chrom}-deferred.bam"),
        unliftable=temp("temp/{sm}/leviosam2/{sm}-{chrom}-unliftable.bam"),
    threads: MAX_THREADS
    resources:
        mem_mb=MAX_THREADS * 4 * 1024,
        runtime=16 * 60,
    conda:
        DEFAULT_ENV
    params:
        # maximum number of CIGAR opts to change, also the max gap size that can be spanned
        G=config.get("levio_G", 100_000),
        # Using -S clipped_frac 0.05 means when a read has >5% clipped bases, it is deferred. A lower value is more stringent (by deferring more reads).
        # aln_score is the minumum score before the alignment is lifted over
        S=config.get(
            "levio_S",
            f"-S mapq:0 -S hdist:{100_000} -S isize:{100_000} -S clipped_frac:0.95 -S aln_score:100",
        ),
        # number of reads per thread
        T=config.get("levio_T", 4 * 256),
    shell:
        """
        PRE="temp/{wildcards.sm}/leviosam2/{wildcards.sm}-{wildcards.chrom}"
        {LEVIO_EXE} lift -t {threads} -a {input.cram} \
            -T {params.T} -G {params.G} {params.S} \
            -C {input.levio_index} -p $PRE -f {input.ref} -m -O bam
        """
        # ^ bam is the only option, no CRAM.
        #samtools view -@ {threads} -u {input.cram} \


#
# This step sorted the leviosam2 output and fixes some tags in the CRAM file.
#
# Specifically, the MAPQ is reset to 60 for all reads that were previously aligned to the DSA.
# And the XS tag is set to zero for all reads that were aligned to the DSA.
# This is a hueristic that we may need to return to in the future.
#
# Other tags and fields like CIGAR, bitflags, and MD are correctly updated by
# leviosam2 during liftover.
#
rule leviosam2_sorted:
    input:
        lifted=rules.leviosam2.output.lifted,
        ref=REF,
    output:
        bam=temp("temp/{sm}/leviosam2/{sm}-{chrom}-sorted.bam"),
    threads: SORT_THREADS
    resources:
        mem_mb=SORT_THREADS * 4 * 1024,
        runtime=16 * 60,
    conda:
        DEFAULT_ENV
    shell:
        """
        samtools sort {input.lifted} \
            -@ {threads} -m 3G \
            -o {output.bam}
        """


# params:
# reset_mapq=workflow.source_path("../scripts/reset-mapq.py"),
# python {params.reset_mapq} -t {threads} {input.lifted} \
