rule decorate_fibers_chromosome:
    input:
        bed=rules.merge_model_results.output.bed,
        tbi=rules.index_model_results.output.tbi,
    output:
        bed=temp("temp/{sm}/decorate/{chrom}.bed.gz"),
        decorated=temp("temp/{sm}/decorate/{chrom}.dec.bed.gz"),
    threads: 8
    resources:
        mem_mb=get_large_mem_mb,
    conda:
        "../envs/python.yaml"
    params:
        script=workflow.source_path("../scripts/decorated-bed12.py"),
    shell:
        """
        INBED={resources.tmpdir}/tmp.{wildcards.sm}.{wildcards.chrom}.bed.gz
        OUTBED={resources.tmpdir}/tmp.out.{wildcards.sm}.{wildcards.chrom}.bed
        tabix -h {input.bed} {wildcards.chrom} | bgzip -@ {threads} > $INBED
        python {params.script} -v 1 $INBED $OUTBED \
            | sort -k1,1 -k2,2n -k3,3n -k4,4 \
            | bgzip -@ {threads} \
        > {output.decorated}

        # sort the other file 
        sort -k1,1 -k2,2n -k3,3n -k4,4 $OUTBED \
            | bgzip -@ {threads} \
        > {output.bed}

        rm $INBED $OUTBED
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
