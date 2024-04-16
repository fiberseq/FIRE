#
# Coverage calculations
#
rule genome_bedgraph:
    input:
        bam=rules.merged_fire_bam.output.bam,
        fai=ancient(f"{ref}.fai"),
    output:
        bg="results/{sm}/coverage/{sm}.bed.gz",
        csi="results/{sm}/coverage/{sm}.bed.gz.csi",
    threads: 16
    shadow:
        "shallow"
    conda:
        default_env
    shell:
        """ 
        mosdepth -t {threads} tmp {input.bam}
        mv tmp.per-base.bed.gz {output.bg}
        mv tmp.per-base.bed.gz.csi {output.csi}
        """


rule coverage:
    input:
        bam=rules.merged_fire_bam.output.bam,
    output:
        cov="results/{sm}/coverage/{sm}.median.coverage.txt",
        minimum="results/{sm}/coverage/{sm}.minimum.coverage.txt",
        maximum="results/{sm}/coverage/{sm}.maximum.coverage.txt",
    conda:
        default_env
    threads: 16
    params:
        n_sd=coverage_within_n_sd,
        mincov=min_coverage,
    shell:
        """
        samtools depth -@ {threads} {input.bam} | datamash median 3 > {output.cov}
        MEDIAN=$(cat {output.cov})

        # calculate minimum and maximum coverage        
        echo $MEDIAN \
            | awk '{{print int($0 + {params.n_sd} * sqrt($0) + 0.5) }}' \
            > {output.maximum}

        echo $MEDIAN \
            | awk '{{print int($0 - {params.n_sd} * sqrt($0) + 0.5) }}' \
            | awk '{{if $0 < {params.mincov} {{print {params.mincov}}}; else {{print $0}}}}' \
            > {output.minimum}

        echo "Median coverage: $MEDIAN"
        echo "Minimum coverage: $(cat {output.minimum})"
        echo "Maximum coverage: $(cat {output.maximum})"
        """

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
