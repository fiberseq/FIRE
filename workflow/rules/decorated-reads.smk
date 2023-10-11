rule decorate_reads_chromosome:
    input:
        bed=rules.merge_model_results.output.bed,
    output:
        bed=temp("temp/{sm}/decorate/{chrom}.bed.gz"),
        decorated=temp("temp/{sm}/decorate/{chrom}.dec.bed.gz"),
    threads: 4
    conda:
        conda
    params:
        script=workflow.source_path("../scripts/decorated-bed12.py"),
    shell:
        """
        tabix -h {input.bed} {wildcards.chrom} \
            | python {params.script} - {output.bed} \
            |  sort -k1,1 -k2,2n -k3,3n -k4,4 \
            | bgzip -@ {threads} \
        > {output.decorated}
        """


rule decorate_fibers:
    input:
        bed=expand(
            rules.decorate_reads_chromosome.output.bed,
            chrom=get_chroms(),
            allow_missing=True,
        ),
        decorated=expand(
            rules.decorate_reads_chromosome.output.decorated,
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
        dec_as=workflow.source_path("../templates/decoration.as"),
        bed_as=workflow.source_path("../templates/bed12_filter.as"),
    shell:
        """
        cat {input.bed} > {output.bed}
        bedToBigBed {output.bed} {input.fai} {output.bb}
        cat {input.decorated} > {output.decorated}
        bedToBigBed {output.decorated} {input.fai} {output.bbd}
        """
