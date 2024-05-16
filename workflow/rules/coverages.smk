#
# Coverage calculations
#
rule genome_bedgraph:
    input:
        cram=rules.merged_fire_bam.output.cram,
        crai=rules.merged_fire_bam.output.crai,
    output:
        bg="results/{sm}/coverage/{sm}.bed.gz",
        tbi="results/{sm}/coverage/{sm}.bed.gz.tbi",
    threads: 16
    shadow:
        "minimal"
    conda:
        DEFAULT_ENV
    benchmark:
        "results/{sm}/benchmarks/genome_bedgraph/{sm}.txt"
    shell:
        """ 
        mosdepth -t {threads} tmp {input.cram}
        zcat tmp.per-base.bed.gz \
            | LC_ALL=C sort --parallel={threads} -k1,1 -k2,2n -k3,3n -k4,4  \
            | bgzip -@ {threads} \
        > {output.bg}
        tabix -f -p bed {output.bg}
        """


rule coverage:
    input:
        bg=rules.genome_bedgraph.output.bg,
    output:
        cov="results/{sm}/coverage/{sm}.median.coverage.txt",
        minimum="results/{sm}/coverage/{sm}.minimum.coverage.txt",
        maximum="results/{sm}/coverage/{sm}.maximum.coverage.txt",
    conda:
        "../envs/python.yaml"
    threads: 1
    resources:
        mem_mb=48 * 1024,
    benchmark:
        "results/{sm}/benchmarks/coverage/{sm}.txt"
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
        bed=temp("temp/{sm}/coverage/{chrom}.fiber-locations.bed.gz"),
    threads: 8
    conda:
        DEFAULT_ENV
    shell:
        """
        # get fiber locations
        (samtools view -@ {threads} -F 2308 -u {input.cram} {wildcards.chrom} \
            | ft extract -t {threads} -s --all - \
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
        bed="results/{sm}/coverage/fiber-locations.bed.gz",
        bed_tbi="results/{sm}/coverage/fiber-locations.bed.gz.tbi",
        filtered="results/{sm}/coverage/filtered-for-coverage/fiber-locations.bed.gz",
        filtered_tbi="results/{sm}/coverage/filtered-for-coverage/fiber-locations.bed.gz.tbi",
    threads: 4
    conda:
        DEFAULT_ENV
    shell:
        """
        cat {input.fibers} > {output.bed}
        tabix -f -p bed {output.bed}
        
        # get filtered fiber locations
        MIN=$(cat {input.minimum})
        MAX=$(cat {input.maximum})
        bedtools intersect -header -sorted -v -f 0.2 \
            -a {output.bed} \
            -b <(zcat {input.bg} | awk -v MAX="$MAX" -v MIN="$MIN" '$4 <= MIN || $4 >= MAX') \
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
        bed="results/{sm}/coverage/exclude-from-shuffles.bed.gz",
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
        bed="results/{sm}/coverage/unreliable-coverage-regions.bed.gz",
        bed_tbi="results/{sm}/coverage/unreliable-coverage-regions.bed.gz.tbi",
        bb="results/{sm}/trackHub/bb/unreliable-coverage-regions.bb",
    threads: 4
    conda:
        DEFAULT_ENV
    shell:
        """
        FILE={output.bed}
        TMP="${{FILE%.*}}"

        MIN=$(cat {input.minimum})
        MAX=$(cat {input.maximum})
        zcat {input.bg} \
            | awk '$4>0' \
            | awk -v MAX="$MAX" -v MIN="$MIN" '$4 <= MIN || $4 >= MAX' \
            | bedtools merge -i - \
        > $TMP
        bedToBigBed $TMP {input.fai} {output.bb}
        # compress 
        bgzip -f -@ {threads} $TMP
        tabix -f -p bed {output.bed}
        """
