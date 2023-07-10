rule make_fdr_d4:
    input:
        fai=ancient(f"{ref}.fai"),
        bed=rules.merge_model_results.output.bed,
        tbi=rules.index_model_results.output.tbi,
    output:
        d4=temp("temp/{sm}/{hp}/{chrom}.fdr.coverages.d4"),
        bed=temp("temp/{sm}/{hp}/{chrom}.fdr.coverages.bed"),
    benchmark:
        "benchmarks/{sm}/{hp}/{chrom}.fdr.d4.tsv"
    threads: 4
    resources:
        mem_mb=get_mem_mb,
    conda:
        conda
    shell:
        """
        printf "{wildcards.chrom}\t0\t1\tfake\t100\t+\t0\t1\t230,230,230\t1.0\n" > {output.bed}
        printf "{wildcards.chrom}\t0\t1\tfake\t100\t+\t0\t1\t147,112,219\t1.0\n" >> {output.bed}
        printf "{wildcards.chrom}\t0\t1\tfake\t3\t+\t0\t1\t255,0,0\t0.03\n" >> {output.bed}
        tabix {input.bed} {wildcards.chrom} >> {output.bed}
        
        fibertools -v bed2d4 \
            --chromosome {wildcards.chrom} \
            -g {input.fai} \
            -c score \
            {output.bed} {output.d4}
        """


rule make_fdr_peaks:
    input:
        fai=ancient(f"{ref}.fai"),
        d4=rules.make_fdr_d4.output.d4,
    output:
        d4="temp/{sm}/{hp}/{chrom}.fdr.peaks.d4",
    benchmark:
        "benchmarks/{sm}/{hp}/{chrom}.fdr.peaks.tsv"
    threads: 4
    resources:
        mem_mb=get_mem_mb,
    conda:
        conda
    shell:
        """
        fibertools -v bed2d4 \
            --chromosome {wildcards.chrom} \
            -g {input.fai} \
            -c score \
            -q {input.d4} {output.d4}
        """


rule fdr_bed:
    input:
        peaks=rules.make_fdr_peaks.output.d4,
    output:
        bed="results/{sm}/{hp}/chromosomes/{chrom}.fdr.peaks.and.coverages.bed.gz",
    threads: 4
    resources:
        mem_mb=get_mem_mb,
    conda:
        conda
    shell:
        """
        d4tools view {input.peaks} {wildcards.chrom} \
          | bgzip -@ {threads} \
        > {output.bed}
        #sort --parallel={threads} -S 2G -k1,1 -k2,2n -k3,3n -k4,4 
        """


rule chromosome_coverage_tracks:
    input:
        bed=rules.fdr_bed.output.bed,
    output:
        bed=temp("temp/{sm}/{hp}/trackHub/bw/{chrom}.{types}.cov.bed"),
    threads: 4
    params:
        col=lambda wc: types_to_col[wc.types],
    conda:
        conda
    resources:
        mem_mb=get_mem_mb,
    shell:
        """
        bgzip -cd -@ {threads} {input.bed} | cut -f 1,2,3,{params.col} > {output.bed}
        """


rule coverage_tracks:
    input:
        beds=expand(
            rules.chromosome_coverage_tracks.output.bed,
            chrom=get_chroms(),
            allow_missing=True,
        ),
        fai=ancient(f"{ref}.fai"),
    output:
        bed=temp("temp/{sm}/{hp}/trackHub/bw/{types}.cov.bed"),
        bw="results/{sm}/trackHub/bw/{hp}.{types}.bw",
    threads: 4
    conda:
        conda
    resources:
        mem_mb=get_mem_mb,
    shell:
        """
        cat {input.beds} > {output.bed}
        bedGraphToBigWig {output.bed} {input.fai} {output.bw}
        """


rule merged_fdr_track:
    input:
        beds=expand(
            rules.fdr_bed.output.bed,
            chrom=get_chroms(),
            allow_missing=True,
        ),
        fai=ancient(f"{ref}.fai"),
    output:
        bed="results/{sm}/{hp}/fdr.peaks.and.coverages.bed.gz",
        tbi="results/{sm}/{hp}/fdr.peaks.and.coverages.bed.gz.tbi",
    threads: 4
    conda:
        conda
    resources:
        mem_mb=get_mem_mb,
    shell:
        """
        cat {input.beds} > {output.bed}
        tabix -p bed {output.bed}
        """


rule chromosome_fdr_tracks:
    input:
        bed=rules.fdr_bed.output.bed,
    output:
        bed=temp("temp/{sm}/{hp}/trackHub/bw/{chrom}.fdr.{fdr}.bed"),
    threads: 4
    conda:
        conda
    resources:
        mem_mb=get_mem_mb,
    shell:
        """
        bgzip -cd -@ {threads} {input.bed} | cut -f 1-4 | awk '$4 > {wildcards.fdr}' > {output.bed}
        """


rule fdr_tracks:
    input:
        beds=expand(
            rules.chromosome_fdr_tracks.output.bed,
            chrom=get_chroms(),
            allow_missing=True,
        ),
        fai=ancient(f"{ref}.fai"),
    output:
        bw="results/{sm}/trackHub/bw/fdr.{hp}.{fdr}.bw",
        bed=temp("temp/{sm}/{hp}/trackHub/bw/fdr.{fdr}.bed"),
    threads: 4
    conda:
        conda
    resources:
        mem_mb=get_mem_mb,
    shell:
        """
        head -n 1 {input.fai} | awk '{{print $1"\t0\t1\t0"}}' > {output.bed}
        cat {input.beds} | awk 'NF > 2'  >> {output.bed}
        bedGraphToBigWig {output.bed} {input.fai} {output.bw}
        """


rule average_coverage:
    input:
        median=rules.genome_bedgraph.output.median,
    output:
        cov="results/{sm}/coverage/{sm}.median.coverage.txt",
    run:
        find_median_coverage(input["median"], outfile=output["cov"])


rule binned_fdr_calls:
    input:
        bed=rules.merge_model_results.output.bed,
        tbi=rules.index_model_results.output.tbi,
    output:
        beds=temp(
            expand(
                "temp/{sm}/{hp}/{chrom}.bin.{bin}.bed", bin=bins, allow_missing=True
            )
        ),
    benchmark:
        "benchmarks/{sm}/{hp}/{chrom}.fdr.d4.tsv"
    threads: 4
    resources:
        mem_mb=get_large_mem_mb,
    conda:
        conda
    shell:
        """
        ((zcat {input.bed} | head -n 1) || true; tabix {input.bed} {wildcards.chrom}) \
            | fibertools -v bin - --outs {output.beds}
        """


rule merge_binned_fdr_calls:
    input:
        beds=expand(
            "temp/{sm}/{hp}/{chrom}.bin.{bin}.bed",
            chrom=get_chroms(),
            allow_missing=True,
        ),
        fai=f"{ref}.fai",
    output:
        bed=temp("temp/{sm}/{hp}/chromosomes/{bin}.bed"),
        bb="results/{sm}/trackHub/bins/{hp}.bin.{bin}.bed.bb",
    threads: 1
    resources:
        mem_mb=get_mem_mb,
    conda:
        conda
    params:
        chrom=get_chroms()[0],
    shell:
        """
        printf "{params.chrom}\t0\t1\tfake\t100\t+\t0\t1\t230,230,230\n" > {output.bed}
        cat {input.beds} | awk 'NF == 9' >> {output.bed}
        bedToBigBed {output.bed} {input.fai} {output.bb}
        """


rule fire_sites:
    input:
        bed=expand(rules.merge_model_results.output.bed, hp="all", allow_missing=True),
    output:
        bed="results/{sm}/FIRE.bed.gz",
    threads: 8
    conda:
        conda
    params:
        min_fdr=min_fire_fdr,
    shell:
        """
        bgzip -cd -@{threads} {input.bed} \
            | awk '$5<={params.min_fdr}' \
            | bgzip -@{threads} \
            > {output.bed}
        """


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


# peak calling
rule peak_calls_per_chromosome:
    input:
        bed=rules.merged_fdr_track.output.bed,
        tbi=rules.merged_fdr_track.output.tbi,
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


rule percent_in_clusters:
    input:
        bed=rules.clustering_vs_null.output.bed,
        fdr=expand(rules.merged_fdr_track.output.bed, hp="all", allow_missing=True),
    output:
        txt="results/{sm}/percent-in-clusters.txt",
    threads: 8
    conda:
        conda
    params:
        script=workflow.source_path("../scripts/percent-in-clusters.py"),
    shell:
        """
        python {params.script} {input.bed} {input.fdr} {output.txt}
        """


rule n_peaks:
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
        bed=rules.n_peaks.output.bed,
        cov=rules.average_coverage.output.cov,
        cov_bed=expand(rules.merged_fdr_track.output.bed, hp="all", allow_missing=True),
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
        bedToBigBed {output.bed} {input.fai} {output.bb}
        """


rule hap_peaks:
    input:
        bed=rules.fire_with_coverage.output.bed,
        cov=rules.average_coverage.output.cov,
        h1=expand(rules.merged_fdr_track.output.bed, hp="hap1", allow_missing=True),
        h2=expand(rules.merged_fdr_track.output.bed, hp="hap2", allow_missing=True),
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


rule hap_differences_track:
    input:
        bed9=rules.hap_differences.output.bed9,
        fai=f"{ref}.fai",
    output:
        bed=temp("temp/{sm}/hap_differences/temp.bed"),
        bb="results/{sm}/trackHub/bb/hap_differences.bb",
    threads: 1
    resources:
        mem_mb=get_mem_mb,
    conda:
        conda
    params:
        chrom=get_chroms()[0],
    shell:
        """
        printf "{params.chrom}\t0\t1\tfake\t100\t+\t0\t1\t230,230,230\\n" > {output.bed}
        bedtools sort -i {input.bed9} >> {output.bed}
        bedToBigBed {output.bed} {input.fai} {output.bb}
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


rule trackhub:
    input:
        fai=ancient(f"{ref}.fai"),
        fire=rules.fire_bw.output.bb,
        cov=rules.average_coverage.output.cov,
        hap_diffs=rules.hap_differences.output.bed,
        hap_diffs2=rules.hap_differences_track.output.bb,
        bed=expand(rules.merge_model_results.output.bed, hp=haps, allow_missing=True),
        bw=expand(rules.fdr_tracks.output.bw, hp=haps, fdr=[100], allow_missing=True),
        fdr=expand(
            rules.coverage_tracks.output.bw, hp=haps, types="fdr", allow_missing=True
        ),
        acc=expand(
            rules.coverage_tracks.output.bw, hp=haps, types="acc", allow_missing=True
        ),
        link=expand(
            rules.coverage_tracks.output.bw, hp=haps, types="link", allow_missing=True
        ),
        nuc=expand(
            rules.coverage_tracks.output.bw, hp=haps, types="nuc", allow_missing=True
        ),
    output:
        hub="results/{sm}/trackHub/hub.txt",
    benchmark:
        "benchmarks/{sm}/trackhub.tsv"
    resources:
        load=get_load,
    threads: 4
    conda:
        conda
    params:
        ref=ref_name,
    shell:
        """
        fibertools -v trackhub \
          -r {params.ref} \
          --sample {wildcards.sm} \
          -t results/{wildcards.sm}/trackHub \
          --average-coverage $(cat {input.cov}) \
          {input.fai} \
          --bw {input.acc} {input.link} {input.nuc} {input.bw} {input.fdr}
        """
