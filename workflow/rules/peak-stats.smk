#
# stats in peaks
#
rule clustering_vs_null:
    input:
        bed=rules.fire_sites.output.bed,
        fai=ancient(FAI),
    output:
        tmp=temp("temp/{sm}/tmp.pre.calls.bed"),
        null=temp("temp/{sm}/null.calls.bed"),
        bed="results/{sm}/clustering-vs-null.bed.gz",
    threads: 8
    conda:
        DEFAULT_ENV
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


rule fires_in_peaks:
    input:
        fire=rules.fire_sites.output.bed,
        exclude=rules.unreliable_coverage_regions.output.bed,
        peaks=rules.fdr_peaks_by_fire_elements.output.bed,
    output:
        tmp=temp("temp/{sm}/tmp.FIREs-in-peaks.bed"),
        txt="results/{sm}/tables/FIREs-in-peaks.txt",
    threads: 8
    conda:
        DEFAULT_ENV
    params:
        script=workflow.source_path("../scripts/percent-in-clusters.sh"),
    shell:
        """
        bedtools intersect -sorted -a {input.fire} -b {input.exclude} -v > {output.tmp}

        echo "Total # of FIREs within normal coverage regions" >> {output.txt}
        wc -l {output.tmp} >> {output.txt}

        echo "# of FIREs within peaks" >> {output.txt}
        bedtools intersect -sorted -u -a {input.fire} -b {input.peaks} | wc -l >> {output.txt} 
        """


rule hap_differences:
    input:
        bed=rules.fdr_peaks_by_fire_elements.output.bed,
    output:
        fig1=report(
            "results/{sm}/hap1-vs-hap2/hap1-vs-hap2.pdf",
            category="Haplotype selectivity",
        ),
        fig2=report(
            "results/{sm}/hap1-vs-hap2/hap1-vs-hap2-volcano.pdf",
            category="Haplotype selectivity",
        ),
        bed="results/{sm}/hap1-vs-hap2/FIRE.hap.differences.bed",
        bed9="results/{sm}/hap1-vs-hap2/FIRE.hap.differences.bed9",
    threads: 8
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/hap-diffs.R"
