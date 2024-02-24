#
# Coverage calculations
#
rule genome_bedgraph:
    input:
        bam=rules.merged_fire_bam.output.bam,
        fai=ancient(f"{ref}.fai"),
    output:
        d4="results/{sm}/coverage/{sm}.d4",
        bg="results/{sm}/coverage/{sm}.bed.gz",
        median="results/{sm}/coverage/{sm}.median.chromosome.coverage.bed",
    threads: 16
    conda:
        default_env
    shell:
        """ 
        d4tools create -t {threads} -Azr {input.fai} {input.bam} {output.d4}
        d4tools view {output.d4} | bedtools sort | bgzip -@ {threads} > {output.bg}
        zcat {output.bg} \
            | awk '$4>0' \
            | datamash -g 1 min 2 max 3 median 4 \
        > {output.median}
        """


rule coverage:
    input:
        median=rules.genome_bedgraph.output.median,
    output:
        cov="results/{sm}/coverage/{sm}.median.coverage.txt",
        minimum="results/{sm}/coverage/{sm}.minimum.coverage.txt",
        maximum="results/{sm}/coverage/{sm}.maximum.coverage.txt",
    conda:
        "../envs/python.yaml"
    params:
        chrom=get_chroms(),
        coverage_within_n_sd=coverage_within_n_sd,
        min_coverage=min_coverage,
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
            | awk -v MAX="$MAX" -v MIN="$MIN" '$4 <= MIN || $4 >= MAX' \
            | bedtools merge -i - \
        > $TMP
        bedToBigBed $TMP {input.fai} {output.bb}
        # compress 
        bgzip -f -@ {threads} $TMP
        tabix -f -p bed {output.bed}
        """
