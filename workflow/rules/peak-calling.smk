#
# shane peak calling
#
rule peak_calls_per_chromosome:
    input:
        bed=rules.merged_fire_track.output.bed,
        tbi=rules.merged_fire_track.output.tbi,
    output:
        bed="temp/{sm}/{hp}/{chrom}.peak.calls.bed",
    benchmark:
        "benchmarks/{sm}/{hp}/{chrom}.peak.calls.tsv"
    threads: 1
    conda:
        conda
    params:
        script=workflow.source_path("../scripts/smooth-and-peak-call.tcsh"),
        details=workflow.source_path("../scripts/per-chrom-peakcall-details.tcsh"),
        fdr=80,
    shell:
        """
        chmod +x {params.script} {params.details}
        {params.script} {wildcards.chrom} {input.bed} {params.fdr} {output.bed} {params.details}
        """


rule merge_peak_calls:
    input:
        beds=expand(
            rules.peak_calls_per_chromosome.output.bed,
            chrom=get_chroms(),
            allow_missing=True,
        ),
    output:
        bed="results/{sm}/{hp}/peak.calls.bed",
    threads: 1
    conda:
        conda
    shell:
        """
        printf \
            "#chromosome\tstart\tstop\tID\tAvgFDR\tSumFDR\tMaxFDR\tWaveletSummit\tFDR_at_WaveletSummit\tSmooth_at_WaveletSummit\n" \
        > {output.bed}
        cat {input.beds} >> {output.bed}
        """


rule fire_peaks:
    input:
        txt=rules.percent_in_clusters.output.txt,
        bed=expand(rules.merge_peak_calls.output.bed, hp="all", allow_missing=True),
    output:
        bed=temp("temp/{sm}/FIRE.peaks.bed"),
    threads: 8
    conda:
        conda
    shell:
        """
        MIN_FDR=$(hck -F min_fdr {input.txt} | tail -n 1)
        echo $MIN_FDR
        awk -v min_fdr=$MIN_FDR '$7 >= min_fdr' {input.bed} > {output.bed}
        """


rule fire_with_coverage:
    input:
        bed=rules.fire_peaks.output.bed,
        cov=rules.average_coverage.output.cov,
        cov_bed=expand(rules.merged_fire_track.output.bed, hp="all", allow_missing=True),
    output:
        bed="results/{sm}/FIRE.peaks.with.coverage.bed",
    threads: 8
    conda:
        conda
    shell:
        """
        COV=$(cat {input.cov})
        MIN=$(echo "$COV" | awk '{{print $1-5*sqrt($1)}}')
        MAX=$(echo "$COV" | awk '{{print $1+5*sqrt($1)}}')
        printf "#ct\tst\ten\t" > {output.bed}
        printf "peak_ct\tpeak_st\tpeak_en\t" >> {output.bed}
        printf "peak_fdr\tpeak_acc\tpeak_link\tpeak_nuc\t" >> {output.bed}
        printf "sample\tcov\\n" >> {output.bed}
        paste \
            <(bedmap --delim '\\t' --echo --max-element \
                <(cut -f 1-3 {input.bed} | tail -n +2 ) \
                <(zcat {input.cov_bed} | tail -n +2) \
            ) \
            | sed "s/$/\t{wildcards.sm}\t$COV/g" \
            | awk -v cov=$COV -v MIN=$MIN -v MAX=$MAX '$8+$9+$10 > MIN && $8+$9+$10 < MAX' \
        >> {output.bed}
        """


rule fire_bw:
    input:
        bed=rules.fire_with_coverage.output.bed,
        fai=ancient(f"{ref}.fai"),
    output:
        bb="results/{sm}/trackHub/bb/FIRE.bb",
        bed=temp("temp/{sm}/trackHub/bb/FIRE.temp.bed"),
    threads: 8
    conda:
        conda
    shell:
        """
        grep -v AvgFDR {input.bed} | cut -f 1-3 > {output.bed}
        if [[ -s {output.bed} ]]; then
            bedToBigBed {output.bed} {input.fai} {output.bb}
        else
          touch {output.bb}
        fi
        """


rule hap_peaks:
    input:
        bed=rules.fire_with_coverage.output.bed,
        cov=rules.average_coverage.output.cov,
        h1=expand(rules.merged_fire_track.output.bed, hp="hap1", allow_missing=True),
        h2=expand(rules.merged_fire_track.output.bed, hp="hap2", allow_missing=True),
    output:
        bed=temp("temp/{sm}/hap1-vs-hap2/FIRE.hap.peaks.bed"),
    threads: 8
    conda:
        conda
    shell:
        """
        COV=$(cat {input.cov})
        printf "#ct\tst\ten\t" > {output.bed}
        printf "hap1_ct\thap1_st\thap1_en\t" >> {output.bed}
        printf "hap1_fdr\thap1_acc\thap1_link\thap1_nuc\t" >> {output.bed}
        printf "hap2_ct\thap2_st\thap2_en\t" >> {output.bed}
        printf "hap2_fdr\thap2_acc\thap2_link\thap2_nuc\t" >> {output.bed}
        printf "sample\tcov\\n" >> {output.bed}
        paste \
            <(bedmap --delim '\\t' --echo --max-element \
                <(cut -f 1-3 {input.bed} | tail -n +2 ) \
                <(zcat {input.h1} | tail -n +2) \
            ) \
            <(bedmap --max-element  \
                <(cut -f 1-3 {input.bed} | tail -n +2) \
                <(zcat {input.h2} | tail -n +2) \
            ) \
            | sed "s/$/\t{wildcards.sm}\t$COV/g" \
        >> {output.bed}
        """


rule hap_differences:
    input:
        bed=rules.hap_peaks.output.bed,
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


#
# stats in peaks
#
rule clustering_vs_null:
    input:
        bed=rules.fire_sites.output.bed,
        fai=ancient(f"{ref}.fai"),
    output:
        tmp=temp("temp/{sm}/acc.calls.bed"),
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


rule peaks_vs_percent:
    input:
        bed=rules.fire_with_coverage.output.bed,
    output:
        fig1="results/{sm}/{sm}.peaks-vs-percent.pdf",
    threads: 8
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/peaks-vs-percent.R"


rule percent_in_clusters:
    input:
        bed=rules.clustering_vs_null.output.bed,
        fire=expand(rules.merged_fire_track.output.bed, hp="all", allow_missing=True),
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
