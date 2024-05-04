#
# Coverage calculations
#
rule genome_bedgraph:
    input:
        bam=rules.merged_fire_bam.output.bam,
        fai=ancient(f"{ref}.fai"),
    output:
        bg="results/{sm}/coverage/{sm}.bed.gz",
        tbi="results/{sm}/coverage/{sm}.bed.gz.tbi",
    threads: 16
    shadow:
        "minimal"
    conda:
        default_env
    benchmark:
        "results/{sm}/benchmarks/genome_bedgraph/{sm}.txt"
    shell:
        """ 
        mosdepth -t {threads} tmp {input.bam}
        ls * 
        zcat tmp.per-base.bed.gz \
            | LC_ALL=C sort -k1,1 -k2,2n -k3,3n -k4,4  \
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
        coverage_within_n_sd=coverage_within_n_sd,
        min_coverage=min_coverage,
        chroms=get_chroms(),
    script:
        "../scripts/cov.py"


#
# fiber locations and coverages
#
rule fiber_locations_chromosome:
    input:
        bam=rules.merged_fire_bam.output.bam,
    output:
        bed=temp("temp/{sm}/coverage/{chrom}.fiber-locations.bed.gz"),
    threads: 8
    conda:
        default_env
    shell:
        """
        # get fiber locations
        (samtools view -@ {threads} -F 2308 -u {input.bam} {wildcards.chrom} \
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
        default_env
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
        fai=ancient(f"{ref}.fai"),
    output:
        bed="results/{sm}/coverage/exclude-from-shuffles.bed.gz",
    threads: 4
    conda:
        default_env
    params:
        exclude=excludes,
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
        fai=ancient(f"{ref}.fai"),
    output:
        bed="results/{sm}/coverage/unreliable-coverage-regions.bed.gz",
        bed_tbi="results/{sm}/coverage/unreliable-coverage-regions.bed.gz.tbi",
        bb="results/{sm}/trackHub/bb/unreliable-coverage-regions.bb",
    threads: 4
    conda:
        default_env
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
