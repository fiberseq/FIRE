#
# Applying the model
#
rule fire:
    input:
        bam=ancient(get_input_bam),
        ref=ancient(REF),
    output:
        cram="results/{sm}/{sm}-fire-{v}-filtered.cram",
        crai="results/{sm}/{sm}-fire-{v}-filtered.cram.crai",
    threads: 32
    resources:
        mem_mb=32 * 1024,
        runtime=600,
    params:
        min_msp=config.get("min_msp", 10),
        min_ave_msp_size=config.get("min_ave_msp_size", 10),
        use_ont=USE_ONT,
        flag=FILTER_FLAG,
    benchmark:
        "results/{sm}/additional-outputs-{v}/benchmarks/{sm}-fire-bam.txt"
    conda:
        DEFAULT_ENV
    shell:
        """
        samtools view -@ {threads} -u -F {params.flag} {input.bam} \
            | {FT_EXE} fire -F {params.flag} -t {threads} \
                {params.use_ont} \
                --min-msp {params.min_msp} \
                --min-ave-msp-size {params.min_ave_msp_size} \
                --skip-no-m6a \
                - - \
            | samtools view -C -@ {threads} -T {input.ref} \
                --output-fmt-option embed_ref=1 \
            | samtools view -C -@ {threads} -T {input.ref} \
                --output-fmt-option embed_ref=1 \
                --input-fmt-option required_fields=0x1bff \
                --write-index -o {output.cram}
        """


rule fire_sites_chrom:
    input:
        cram=rules.fire.output.cram,
    output:
        bed=temp("temp/{sm}/chrom/{v}-{chrom}.sorted.bed.gz"),
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
        bed="results/{sm}/additional-outputs-{v}/fire-peaks/{sm}-{v}-fire-elements.bed.gz",
    threads: 1
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
