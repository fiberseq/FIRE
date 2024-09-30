#
# Applying the model
#
rule fire:
    input:
        bam=ancient(get_input_bam),
    output:
        bam=temp("temp/{sm}/fire/{chrom}.fire.bam"),
    threads: 8
    resources:
        mem_mb=8 * 1024,
    params:
        min_msp=config.get("min_msp", 10),
        min_ave_msp_size=config.get("min_ave_msp_size", 10),
        use_ont=USE_ONT,
        flag=FILTER_FLAG,
    conda:
        DEFAULT_ENV
    shell:
        """
        samtools view -F {params.flag} -u -@ {threads} {input.bam} {wildcards.chrom} \
            | {FT_EXE} fire -t {threads} \
                {params.use_ont} \
                --min-msp {params.min_msp} \
                --min-ave-msp-size {params.min_ave_msp_size} \
                --skip-no-m6a \
                - {output.bam}
        """


rule merged_fire_bam:
    input:
        ref=ancient(REF),
        fai=ancient(FAI),
        bams=expand(rules.fire.output.bam, chrom=get_chroms(), allow_missing=True),
    output:
        cram="results/{sm}/{sm}-fire.cram",
        crai="results/{sm}/{sm}-fire.cram.crai",
    threads: 16
    resources:
        mem_mb=16 * 1024,
        runtime=300,
    conda:
        DEFAULT_ENV
    benchmark:
        "results/{sm}/additional-outputs/benchmarks/{sm}-merged-fire-bam.txt"
    shell:
        """
        samtools merge -@ {threads} -u {input.bams} -o - \
            | samtools view -C -@ {threads} -T {input.ref} \
                --output-fmt-option embed_ref=1 \
            | samtools view -C -@ {threads} -T {input.ref} \
                --output-fmt-option embed_ref=1 \
                --input-fmt-option required_fields=0x1bff \
                --write-index -o {output.cram}
        # the second samtools view of CRAM file is needed to drop the quality scores
        # this halves the size of the CRAM file
        """


rule fire_sites_chrom:
    input:
        cram=rules.merged_fire_bam.output.cram,
    output:
        bed=temp("temp/{sm}/chrom/{chrom}.sorted.bed.gz"),
    threads: 4
    conda:
        DEFAULT_ENV
    resources:
        mem_mb=16 * 1024,
    params:
        min_fdr=MIN_FIRE_FDR,
    shell:
        """
        samtools view -@ {threads} -u {input.cram} {wildcards.chrom} \
            | {FT_EXE} fire -t {threads} --extract - \
                | LC_ALL=C sort --parallel={threads} \
                    -k1,1 -k2,2n -k3,3n -k4,4 \
                | bioawk -tc hdr '$10<={params.min_fdr}' \
                | (grep '\\S' || true) \
                | (grep -v '^#' || true) \
                | bgzip -@ {threads} \
            > {output.bed}
        """


rule fire_sites:
    input:
        beds=expand(
            rules.fire_sites_chrom.output.bed, chrom=get_chroms(), allow_missing=True
        ),
    output:
        bed="results/{sm}/additional-outputs/fire-peaks/{sm}-fire-elements.bed.gz",
    threads: 8
    conda:
        DEFAULT_ENV
    shell:
        """
        cat {input.beds} > {output.bed}
        """


rule fire_sites_index:
    input:
        bed=rules.fire_sites.output.bed,
    output:
        tbi=rules.fire_sites.output.bed + ".tbi",
    threads: 1
    conda:
        DEFAULT_ENV
    shell:
        """
        tabix -p bed {input.bed}
        """
