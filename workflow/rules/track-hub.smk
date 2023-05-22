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
            chrom=chroms,
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
            chrom=chroms,
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
            rules.chromosome_fdr_tracks.output.bed, chrom=chroms, allow_missing=True
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
        bam=lambda wc: data.loc[wc.sm, "bam"],
    output:
        cov="results/{sm}/average.coverage.txt",
    threads: 16
    resources:
        mem_mb=get_mem_mb,
    conda:
        conda
    params:
        chrom=chroms[0],
        first_n=1_000_000_000,
    shell:
        """
        samtools depth -@ {threads} {input.bam} -r {params.chrom} \
            | head -n {params.first_n} \
            | cut -f 3 \
            | datamash median 1 \
            > {output.cov}
        """


rule trackhub:
    input:
        fai=ancient(f"{ref}.fai"),
        cov=rules.average_coverage.output.cov,
        bed=expand(rules.merge_model_results.output.bed, hp=haps, allow_missing=True),
        bw=expand(
            rules.fdr_tracks.output.bw, hp=haps, fdr=[90, 100], allow_missing=True
        ),
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
            "temp/{sm}/{hp}/{chrom}.bin.{bin}.bed", chrom=chroms, allow_missing=True
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
        chrom=chroms[0],
    shell:
        """
        printf "{params.chrom}\t0\t1\tfake\t100\t+\t0\t1\t230,230,230\n" > {output.bed}
        cat {input.beds} | awk 'NF == 9' >> {output.bed}
        bedToBigBed {output.bed} {input.fai} {output.bb}
        """




# peak calling 
rule peak_calls_per_chromosome:
    input:
        bed=rules.merge_model_results.output.bed,
        tbi=rules.index_model_results.output.tbi,
    output:
        bed="temp/{sm}/{hp}/{chrom}.peak.calls.bed",
    benchmark:
        "benchmarks/{sm}/{hp}/{chrom}.peak.calls.tsv"
    threads: 1
    conda:
        conda
    params:
        script=workflow.source_path("../scripts/qc/smooth-and-peak-call.tcsh"),
        fdr=100,
    shell:
        """
        chmod +x {params.script}
        {params.script} {wildcards.chrom} {input.bed} {params.fdr} {output.bed}
        """

