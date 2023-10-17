#
# stats in peaks
#
rule clustering_vs_null:
    input:
        bed=rules.fire_sites.output.bed,
        fai=ancient(f"{ref}.fai"),
    output:
        tmp=temp("temp/{sm}/tmp.pre.calls.bed"),
        null=temp("temp/{sm}/null.calls.bed"),
        bed="results/{sm}/clustering-vs-null.bed.gz",
    threads: 8
    conda:
        conda
    shell:
        """
        bgzip -cd -@{threads} {input.bed} | cut -f 1-3 > {output.tmp}
        bedtools shuffle -chrom -i {output.tmp} -g {input.fai} > {output.null}

        ( bedtools genomecov -bg -i {output.tmp} -g {input.fai} | sed 's/$/\\tReal/g' ; \
          bedtools genomecov -bg -i {output.null} -g {input.fai} | sed 's/$/\\tNull/g' ) \
            | bedtools sort \
            | bgzip -@ {threads} \
        > {output.bed}
        """

rule percent_in_clusters:
    input:
        bed=rules.clustering_vs_null.output.bed,
        fire=rules.fdr_track_filtered.output.bed,
    output:
        txt="results/{sm}/percent-in-clusters.txt",
    threads: 8
    conda:
        conda
    params:
        script=workflow.source_path("../scripts/percent-in-clusters.sh"),
    shell:
        """
        bash {params.script} {input.bed} {input.fire} {output.txt}
        """

rule hap_differences:
    input:
        bed=rules.fdr_peaks_by_fire_elements.output.bed,
    output:
        fig1="results/{sm}/hap1-vs-hap2/hap1-vs-hap2.pdf",
        fig2="results/{sm}/hap1-vs-hap2/hap1-vs-hap2-volcano.pdf",
        bed="results/{sm}/hap1-vs-hap2/FIRE.hap.differences.bed",
        bed9="results/{sm}/hap1-vs-hap2/FIRE.hap.differences.bed9",
    threads: 8
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/hap-diffs.R"
