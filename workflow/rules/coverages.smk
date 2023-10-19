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
    run:
        if force_coverage is not None:
            coverage = force_coverage
        else:
            df = pd.read_csv(
                input.median,
                sep="\t",
                header=None,
                names=["chr", "start", "end", "coverage"],
            )
            df = df[df.coverage > 0]
            df = df[df["chr"].isin(get_chroms())]
            df = df[~df["chr"].isin(["chrX", "chrY", "chrM", "chrEBV"])]
            df["weight"] = df["end"] - df["start"]
            print(df, file=sys.stderr)
            coverage = weighted_median(df, "coverage", "weight")

        min_coverage = get_min_coverage(coverage)
        max_coverage = get_max_coverage(coverage)
        print(coverage, file=sys.stderr)
        if coverage <= 1:
            raise ValueError(
                f"Median coverage is {coverage}! Did you use the correct reference, or is data missing from most of your genome. If so consider the keep_chromosomes parameter in config.yaml"
            )
        open(output.cov, "w").write(str(round(coverage)) + "\n")
        open(output.minimum, "w").write(str(round(min_coverage)) + "\n")
        open(output.maximum, "w").write(str(round(max_coverage)) + "\n")


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
        filtered="results/{sm}/coverage/filtered-for-coverage/fiber-locations.bed.gz",
        filtered_tbi="results/{sm}/coverage/filtered-for-coverage/fiber-locations.bed.gz.tbi",
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
        conda
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
