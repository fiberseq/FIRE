rule decorate_fibers_chromosome:
    input:
        cram=rules.merged_fire_bam.output.cram,
        crai=rules.merged_fire_bam.output.crai,
    output:
        bed=temp("temp/{sm}/decorate/{chrom}.bed.gz"),
        decorated=temp("temp/{sm}/decorate/{chrom}.dec.bed.gz"),
    threads: 8
    resources:
        mem_mb=get_large_mem_mb,
    conda:
        DEFAULT_ENV
    shell:
        """
        samtools view -@ {threads} -u {input.cram} {wildcards.chrom} \
            | ft track-decorators -t {threads} --bed12 {output.bed} \
            | sort -k1,1 -k2,2n -k3,3n -k4,4 \
            | bgzip -@ {threads} \
        > {output.decorated}
        """


rule decorate_fibers_1:
    input:
        bed=expand(
            rules.decorate_fibers_chromosome.output.bed,
            chrom=get_chroms(),
            allow_missing=True,
        ),
        fai=ancient(FAI),
    output:
        bed="results/{sm}/fiber-calls/fire-fibers.bed.gz",
        bb="results/{sm}/trackHub/bb/fire-fibers.bb",
    benchmark:
        "results/{sm}/benchmarks/decorate_fibers_1/{sm}.txt"
    threads: 1
    resources:
        runtime=240,
    conda:
        DEFAULT_ENV
    params:
        bed_as=workflow.source_path("../templates/bed12_filter.as"),
    shell:
        """
        cat {input.bed} > {output.bed}
        bedToBigBed -allow1bpOverlap -type=bed12+ -as={params.bed_as} \
            {output.bed} {input.fai} {output.bb}
        """


rule decorate_fibers_2:
    input:
        decorated=expand(
            rules.decorate_fibers_chromosome.output.decorated,
            chrom=get_chroms(),
            allow_missing=True,
        ),
        fai=ancient(FAI),
    output:
        decorated=temp("temp/{sm}/fiber-calls/fire-fiber-decorators.bed.gz"),
        bb="results/{sm}/trackHub/bb/fire-fiber-decorators.bb",
    benchmark:
        "results/{sm}/benchmarks/decorate_fibers_2/{sm}.txt"
    threads: 1
    resources:
        runtime=60 * 16,
    conda:
        DEFAULT_ENV
    params:
        dec_as=workflow.source_path("../templates/decoration.as"),
    shell:
        """
        cat {input.decorated} > {output.decorated}
        bedToBigBed -allow1bpOverlap -type=bed12+ -as={params.dec_as} \
            {output.decorated} {input.fai} {output.bb}
        """
