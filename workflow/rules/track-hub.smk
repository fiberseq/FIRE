# TODO
rule percent_accessible:
    input:
        bed=rules.fdr_track_filtered.output.bed,
        fai=ancient(f"{ref}.fai"),
    output:
        tmp=temp("temp/{sm}/{hp}/percent.accessible.bed"),
        bw="results/{sm}/trackHub/bw/{hp}.percent.accessible.bw",
        bed="results/{sm}/{hp}/percent.accessible.bed.gz",
        tbi="results/{sm}/{hp}/percent.accessible.bed.gz.tbi",
    threads: 4
    conda:
        conda
    resources:
        mem_mb=get_mem_mb,
    shell:
        """
        bgzip -@{threads} -cd {input.bed} \
            | awk -v OFS='\t' '$5 + $6 + $7 > 0 {{print $1,$2,$3,100*$5/($5+$6+$7)}}' \
        > {output.tmp}
        
        bedGraphToBigWig {output.tmp} {input.fai} {output.bw}
        
        bgzip -@{threads} -c {output.tmp} > {output.bed}
        tabix -p bed {output.bed}
        """


rule fdr_track_to_bw:
    input:
        bed=rules.fdr_track.output.bed,
        fai=ancient(f"{ref}.fai"),
    output:
        bw="results/{sm}/trackHub/bw/{col}.bw",
        tmp=temp("temp/{sm}/trackHub/bw/{col}.tmp.bed"),
    threads: 4
    conda:
        conda
    shell:
        """
        hck -z -f 1-3 -F {wildcards.col} {input.bed} > {output.tmp} 
        bedGraphToBigWig {output.tmp} {input.fai} {output.bw}
        """


rule fdr_peaks_by_fire_elements_to_bb:
    input:
        bed=rules.fdr_peaks_by_fire_elements.output.bed,
        fai=ancient(f"{ref}.fai"),
    output:
        bb="results/{sm}/trackHub/bw/FDR-FIRE-peaks.bb",
        tmp=temp("temp/{sm}/trackHub/bw/FDR-FIRE-peaks.bb.tmp"),
    threads: 4
    conda:
        conda
    shell:
        """
        hck -z -f 1-3 {input.bed} -F FDR > {output.tmp} 
        bedToBigBed -type=bed4 {output.tmp} {input.fai} {output.bb}
        """


#
# make single molecule viz for trackhub
#
rule binned_fire_calls:
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
        "benchmarks/{sm}/{hp}/{chrom}.fire.d4.tsv"
    threads: 4
    resources:
        mem_mb=get_large_mem_mb,
    conda:
        "../envs/fibertools.yaml"
    shell:
        """
        ((zcat {input.bed} | head -n 1) || true; tabix {input.bed} {wildcards.chrom}) \
            | fibertools -v bin - --outs {output.beds}
        """


rule merge_binned_fire_calls:
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


# TODO
rule trackhub:
    input:
        fai=ancient(f"{ref}.fai"),
        fire=rules.fire_bw.output.bb,
        cov=rules.coverage.output.cov,
        hap_diffs=rules.hap_differences.output.bed,
        hap_diffs2=rules.hap_differences_track.output.bb,
        bed=expand(rules.merge_model_results.output.bed, hp=haps, allow_missing=True),
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
        script=workflow.source_path("../scripts/trackhub.py"),
        max_bins=max_bins,
    shell:
        """
        python {params.script} -v \
          --trackhub-dir results/{wildcards.sm}/trackHub \
          --reference {params.ref} \
          --sample {wildcards.sm} \
          --max-bins {max_bins} \
          --average-coverage $(cat {input.cov}) 
        """
