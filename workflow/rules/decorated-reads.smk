# Number of items to bundle in r-tree [default: 256]
BLOCK_SIZE = 256 * 8
# Number of data points bundled at lowest level [default: 1024]
ITEMS_PER_SLOT = 1024 * 8


rule decorate_fibers_chromosome:
    input:
        cram=rules.merged_fire_bam.output.cram,
        crai=rules.merged_fire_bam.output.crai,
    output:
        bed=temp("temp/{sm}/decorate/{v}-{chrom}.bed.gz"),
        decorated=temp("temp/{sm}/decorate/{v}-{chrom}.dec.bed.gz"),
    threads: 4
    resources:
        mem_mb=get_mem_mb,
    conda:
        DEFAULT_ENV
    shell:
        """
        samtools view -@ {threads} -u {input.cram} {wildcards.chrom} \
            | {FT_EXE} track-decorators -t {threads} --bed12 {output.bed} \
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
        #bed=temp("temp/{sm}/fiber-calls/fire-fibers.bed.gz"),
        bb="results/{sm}/trackHub-{v}/bb/fire-fibers.bb",
    benchmark:
        "results/{sm}/additional-outputs-{v}/benchmarks/decorate_fibers_1/{sm}.txt"
    threads: 8
    resources:
        runtime=240,
    conda:
        DEFAULT_ENV
    params:
        bed_as=workflow.source_path("../templates/bed12_filter.as"),
        nzooms=NZOOMS,
        items_per_slot=ITEMS_PER_SLOT,
        block_size=BLOCK_SIZE,
    shell:
        # bigtools version
        """
        cat {input.bed} \
            | bgzip -cd -@ {threads} \
            | bigtools bedtobigbed \
                --inmemory \
                --block-size {params.block_size} --items-per-slot {params.items_per_slot} \
                --nzooms {params.nzooms} \
                -s start -a {params.bed_as} \
                - {input.fai} {output.bb}
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
        bb="results/{sm}/trackHub-{v}/bb/fire-fiber-decorators.bb",
        #bed=temp("temp/{sm}/trackHub-{v}/bb/fire-fiber-decorators.bed.gz"),
    benchmark:
        "results/{sm}/additional-outputs-{v}/benchmarks/decorate_fibers_2/{sm}.txt"
    threads: 8
    resources:
        runtime=60 * 16,
    conda:
        DEFAULT_ENV
    params:
        dec_as=workflow.source_path("../templates/decoration.as"),
        nzooms=NZOOMS,
        items_per_slot=ITEMS_PER_SLOT,
        block_size=BLOCK_SIZE,
    shell:
        # bigtools version
        # for some reason filtering out NUCs removes the display bug for bigtools
        # at least in my test cases
        """
        cat {input.decorated} \
            | bgzip -cd -@ {threads} \
            | rg -v '^#' \
            | rg -vw 'NUC' \
            | bigtools bedtobigbed \
                --inmemory \
                --block-size {params.block_size} --items-per-slot {params.items_per_slot} \
                --nzooms {params.nzooms} \
                -a {params.dec_as} -s start \
                - {input.fai} {output.bb}
        """


if False:
    # UCSC version
    """
        cat {input.decorated} > {output.bed}
        bedToBigBed \
            -allow1bpOverlap -type=bed12+ -as={params.dec_as} \
            {output.bed} {input.fai} {output.bb}
        """
