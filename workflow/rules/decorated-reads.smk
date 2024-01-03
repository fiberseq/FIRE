rule decorate_fibers_chromosome:
    input:
        bam=rules.merged_fire_bam.output.bam,
    output:
        bed=temp("temp/{sm}/decorate/{chrom}.bed"),
        decorated=temp("temp/{sm}/decorate/{chrom}.dec.bed.gz"),
    threads: 8
    resources:
        mem_mb=get_large_mem_mb,
    conda:
        "../envs/env.yaml"
    shell:
        """
        samtools view -@ {threads} -u {input.bam} {wildcards.chrom} \
            | ft track-decorators -t {threads} --bed12 {output.bed} \
            | sort -k1,1 -k2,2n -k3,3n -k4,4 \
            | bgzip -@ {threads} \
        > {output.decorated}
        """


rule decorate_fibers:
    input:
        bed=expand(
            rules.decorate_fibers_chromosome.output.bed,
            chrom=get_chroms(),
            allow_missing=True,
        ),
        decorated=expand(
            rules.decorate_fibers_chromosome.output.decorated,
            chrom=get_chroms(),
            allow_missing=True,
        ),
        fai=f"{ref}.fai",
    output:
        bed="results/{sm}/fiber-calls/fire-fibers.bed.gz",
        bb="results/{sm}/trackHub/bb/fire-fibers.bb",
        decorated="results/{sm}/fiber-calls/fire-fiber-decorators.bed.gz",
        bbd="results/{sm}/trackHub/bb/fire-fiber-decorators.bb",
    threads: 4
    conda:
        conda
    params:
        bed_as=workflow.source_path("../templates/bed12_filter.as"),
        dec_as=workflow.source_path("../templates/decoration.as"),
    shell:
        """
        cat {input.bed} > {output.bed}
        bedToBigBed -allow1bpOverlap -type=bed12+ -as={params.bed_as} {output.bed} {input.fai} {output.bb}
        cat {input.decorated} > {output.decorated}
        bedToBigBed -allow1bpOverlap -type=bed12+ -as={params.dec_as} {output.decorated} {input.fai} {output.bbd}
        """
