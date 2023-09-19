#
# Coverage calculations
#
rule genome_bedgraph:
    input:
        bam=ancient(lambda wc: data.loc[wc.sm, "bam"]),
        fai=ancient(f"{ref}.fai"),
    output:
        d4="results/{sm}/coverage/{sm}.d4",
        bg="results/{sm}/coverage/{sm}.bed.gz",
        median="results/{sm}/coverage/{sm}.median.chromosome.coverage.bed",
    threads: 16
    conda:
        conda
    shell:
        """ 
        d4tools create -t {threads} -Azr {input.fai} {input.bam} {output.d4}
        d4tools view {output.d4} | bgzip -@ {threads} > {output.bg}
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
    run:
        find_median_coverage(
            input["median"],
            outfile=output["cov"],
            min_out=output["minimum"],
            max_out=output["maximum"],
        )


#
# fiber locations and coverages
#
rule fiber_locations_chromosome:
    input:
        bam=lambda wc: data.loc[wc.sm, "bam"],
    output:
        bed=temp("temp/{sm}/coverage/{chrom}.fiber-locations.bed.gz"),
    threads: 8
    conda:
        conda
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
        filtered="results/{sm}/FDR-peaks/filtered-for-fdr/fiber-locations.bed.gz",
        filtered_tbi="results/{sm}/FDR-peaks/filtered-for-fdr/fiber-locations.bed.gz.tbi",
    threads: 4
    conda:
        conda
    shell:
        """
        cat {input.fibers} > {output.bed}
        tabix -f -p bed {output.bed}
        
        # get filtered fiber locations
        MIN=$(cat {input.minimum})
        MAX=$(cat {input.maximum})
        bedtools intersect -header -v -f 0.2 \
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
        #exclude=excludes,
    output:
        bed="results/{sm}/coverage/exclude-from-shuffles.bed.gz",
    threads: 4
    conda:
        conda
    shell:
        """

        ( \
            bedtools genomecov -bga -i {input.filtered} -g {input.fai} | awk '$4 == 0'; \
            less {input.exclude} \
        ) \
            | cut -f 1-3 \
            | bedtools sort \
            | bedtools merge \
            | bgzip -@ {threads} \
        > {output.bed}
        """