#
# Coverage calculations
#
rule genome_bedgraph:
    input:
        ref=ancient(REF),
        fai=ancient(FAI),
        cram=rules.merged_fire_bam.output.cram,
        crai=rules.merged_fire_bam.output.crai,
    output:
        bg=temp("temp/{sm}/coverage/{sm}-{v}.bed.gz"),
        tbi=temp("temp/{sm}/coverage/{sm}-{v}.bed.gz.tbi"),
    threads: 16
    shadow:
        "minimal"
    conda:
        DEFAULT_ENV
    shell:
        """ 
        mosdepth -f {input.ref} -t {threads} tmp {input.cram}
        bgzip -cd tmp.per-base.bed.gz \
            | LC_ALL=C sort --parallel={threads} -k1,1 -k2,2n -k3,3n -k4,4  \
            | bgzip -@ {threads} \
        > {output.bg}
        tabix -f -p bed {output.bg}
        """


rule coverage:
    input:
        bg=rules.genome_bedgraph.output.bg,
    output:
        cov="results/{sm}/additional-outputs-{v}/coverage/{sm}-{v}-median-coverage.txt",
        minimum="results/{sm}/additional-outputs-{v}/coverage/{sm}-{v}-minimum-coverage.txt",
        maximum="results/{sm}/additional-outputs-{v}/coverage/{sm}-{v}-maximum-coverage.txt",
    conda:
        "../envs/python.yaml"
    threads: 1
    resources:
        mem_mb=64 * 1024,
    benchmark:
        "results/{sm}/additional-outputs-{v}/benchmarks/coverage/{sm}.txt"
    params:
        coverage_within_n_sd=COVERAGE_WITHIN_N_SD,
        min_coverage=MIN_COVERAGE,
        chroms=get_chroms(),
    script:
        "../scripts/cov.py"


#
# fiber locations and coverages
#
rule fiber_locations_chromosome:
    input:
        cram=rules.merged_fire_bam.output.cram,
        crai=rules.merged_fire_bam.output.crai,
    output:
        bed=temp("temp/{sm}/coverage/{v}-{chrom}.fiber-locations.bed.gz"),
    threads: 8
    conda:
        DEFAULT_ENV
    shell:
        """
        # get fiber locations
        (samtools view -@ {threads} -u {input.cram} {wildcards.chrom} \
            | {FT_EXE} extract -t {threads} -s --all - \
            | hck -F '#ct' -F st -F en -F fiber -F strand -F HP ) \
            | (grep -v "^#" || true) \
            | bgzip -@ {threads} \
        > {output.bed}
        """


rule fiber_locations:
    input:
        fibers=expand(
            rules.fiber_locations_chromosome.output.bed,
            chrom=get_chroms(),
            allow_missing=True,
        ),
        bg=rules.genome_bedgraph.output.bg,
        minimum=rules.coverage.output.minimum,
        maximum=rules.coverage.output.maximum,
    output:
        bed=temp("temp/{sm}/coverage/{v}-fiber-locations.bed.gz"),
        bed_tbi=temp("temp/{sm}/coverage/{v}-fiber-locations.bed.gz.tbi"),
        filtered=temp(
            "temp/{sm}/coverage/filtered-for-coverage/{v}-fiber-locations.bed.gz"
        ),
        filtered_tbi=temp(
            "temp/{sm}/coverage/filtered-for-coverage/{v}-fiber-locations.bed.gz.tbi"
        ),
    threads: 4
    conda:
        DEFAULT_ENV
    params:
        max_frac_overlap=0.2,
    shell:
        """
        cat {input.fibers} > {output.bed}
        tabix -f -p bed {output.bed}
        
        # get filtered fiber locations
        MIN=$(cat {input.minimum})
        MAX=$(cat {input.maximum})
        bedtools intersect -header -sorted -v -f {params.max_frac_overlap} \
            -a {output.bed} \
            -b <(bgzip -cd {input.bg} | awk -v MAX="$MAX" -v MIN="$MIN" '$4 <= MIN || $4 >= MAX') \
        | bgzip -@ {threads} \
        > {output.filtered}
        tabix -f -p bed {output.filtered}
        """


#
#
#
rule exclude_from_shuffle:
    input:
        filtered=rules.fiber_locations.output.filtered,
        fai=ancient(FAI),
    output:
        bed="results/{sm}/additional-outputs-{v}/coverage/exclude-from-shuffles.bed.gz",
    threads: 4
    conda:
        DEFAULT_ENV
    params:
        exclude=EXCLUDES,
    shell:
        """

        ( \
            bedtools genomecov -bga -i {input.filtered} -g {input.fai} | awk '$4 == 0'; \
            less {params.exclude} \
        ) \
            | cut -f 1-3 \
            | bedtools sort \
            | bedtools merge \
            | bgzip -@ {threads} \
        > {output.bed}
        """


rule unreliable_coverage_regions:
    input:
        bg=rules.genome_bedgraph.output.bg,
        minimum=rules.coverage.output.minimum,
        maximum=rules.coverage.output.maximum,
        fai=ancient(FAI),
    output:
        bed="results/{sm}/additional-outputs-{v}/coverage/unreliable-coverage-regions.bed.gz",
        bed_tbi="results/{sm}/additional-outputs-{v}/coverage/unreliable-coverage-regions.bed.gz.tbi",
        bb="results/{sm}/trackHub-{v}/bb/unreliable-coverage-regions.bb",
    threads: 4
    params:
        min_len=MIN_UNRELIABLE_COVERAGE_LEN,
        bed3_as=workflow.source_path("../templates/bed3.as"),
    conda:
        DEFAULT_ENV
    shell:
        """
        MIN=$(cat {input.minimum})
        MAX=$(cat {input.maximum})
        bgzip -cd {input.bg} \
            | awk '$4>0' \
            | awk -v MAX="$MAX" -v MIN="$MIN" '$4 <= MIN || $4 >= MAX' \
            | bedtools merge -i - \
            | awk '$3-$2 >= {params.min_len}' \
            | bgzip -@ {threads} \
        > {output.bed}

        # bigbed
        bgzip -cd {output.bed} -@ {threads} \
            | bigtools bedtobigbed \
                -s start -a {params.bed3_as} \
                - {input.fai} {output.bb}

        # index 
        tabix -f -p bed {output.bed}
        """
