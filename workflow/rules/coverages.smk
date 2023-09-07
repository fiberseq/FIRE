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
            out=output["cov"],
            min_out=output["minimum"],
            max_out=output["maximum"],
        )
