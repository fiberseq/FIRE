rule percent_accessible:
    input:
        bed=rules.fdr_track.output.bed,
        fai=ancient(FAI),
    output:
        tmp=temp("temp/{sm}/{hp}/percent.accessible.bed"),
        bw="results/{sm}/trackHub/bw/{hp}.percent.accessible.bw",
        bed="results/{sm}/{hp}/percent.accessible.bed.gz",
        tbi="results/{sm}/{hp}/percent.accessible.bed.gz.tbi",
    threads: 4
    conda:
        DEFAULT_ENV
    resources:
        mem_mb=get_mem_mb,
    params:
        cols=hap_hck_columns,
    shell:
        """
        zcat {input.bed} \
            | hck -f 1-3 {params.cols} \
            | grep -v "^#" \
            | awk -v OFS='\t' '$5 > 0 {{print $1,$2,$3,$4*100/$5}}' \
        > {output.tmp}

        # skip if the file is empty
        if [[ -s {output.tmp} ]]; then
            bedGraphToBigWig {output.tmp} {input.fai} {output.bw}
        else
            touch {output.bw}
        fi
        
        bgzip -@{threads} -c {output.tmp} > {output.bed}
        tabix -p bed {output.bed}
        """


rule element_coverages_bw:
    input:
        bed=rules.element_coverages.output.bed,
        fai=ancient(FAI),
    output:
        tmp=temp("temp/{sm}/trackHub/bw/{hp}.{el_type}.coverage.bed"),
        bw="results/{sm}/trackHub/bw/{hp}.{el_type}.coverage.bw",
    conda:
        DEFAULT_ENV
    shell:
        """
        zcat {input.bed} | hck -f 1-3 -F {wildcards.el_type} > {output.tmp}
        bedGraphToBigWig {output.tmp} {input.fai} {output.bw}
        """


rule fdr_track_to_bw:
    input:
        bed=rules.fdr_track.output.bed,
        fai=ancient(FAI),
    output:
        bw="results/{sm}/trackHub/bw/{col}.bw",
        tmp=temp("temp/{sm}/trackHub/bw/{col}.tmp.bed"),
    threads: 4
    conda:
        DEFAULT_ENV
    shell:
        """
        hck -z -f 1-3 -F {wildcards.col} {input.bed} > {output.tmp} 
        bedGraphToBigWig {output.tmp} {input.fai} {output.bw}
        """


rule fdr_peaks_by_fire_elements_to_bb:
    input:
        bed=rules.fdr_peaks_by_fire_elements.output.bed,
        fai=ancient(FAI),
    output:
        bb="results/{sm}/trackHub/bb/FDR-FIRE-peaks.bb",
    threads: 4
    conda:
        DEFAULT_ENV
    params:
        bedfmt=workflow.source_path("../templates/fire_peak.as"),
    shell:
        """
        zcat {input.bed} \
            | bioawk -tc hdr '{{print $1,$2,$3,"peak-"NR,int($score*10),".",$score,"-1",$log_FDR,int($start/2+$end/2)-$peak_start}}' \
            | bioawk -tc hdr '$5<=1000' \
            | rg -v '^#' \
            | bigtools bedtobigbed \
                -a {params.bedfmt} -s start \
                - {input.fai} {output.bb}
        """


rule hap_differences_track:
    input:
        bed9=rules.hap_differences.output.bed9,
        fai=ancient(FAI),
    output:
        bed=temp("temp/{sm}/hap_differences/temp.bed"),
        bb="results/{sm}/trackHub/bb/hap_differences.bb",
    threads: 1
    resources:
        mem_mb=get_mem_mb,
    conda:
        DEFAULT_ENV
    params:
        chrom=get_chroms()[0],
    shell:
        """
        printf "{params.chrom}\t0\t1\tfake\t100\t+\t0\t1\t230,230,230\\n" > {output.bed}
        bedtools sort -i {input.bed9} >> {output.bed}

        bigtools bedtobigbed -s start {output.bed} {input.fai} {output.bb}
        """


rule trackhub:
    input:
        fai=ancient(FAI),
        fire=rules.fdr_peaks_by_fire_elements_to_bb.output.bb,
        cov=rules.coverage.output.cov,
        hap_diffs=rules.hap_differences_track.output.bb,
        wide=rules.wide_fdr_peaks.output.bb,
        decorators_1=rules.decorate_fibers_1.output.bb,
        decorators_2=rules.decorate_fibers_2.output.bb,
    output:
        hub="results/{sm}/trackHub/hub.txt",
    resources:
        load=get_load,
    threads: 4
    conda:
        "../envs/python.yaml"
    params:
        ref=REF_NAME,
        script=workflow.source_path("../scripts/trackhub.py"),
    shell:
        """
        python {params.script} -v 2 \
          --trackhub-dir results/{wildcards.sm}/trackHub \
          --reference {params.ref} \
          --sample {wildcards.sm} \
          --average-coverage $(cat {input.cov}) 
        """
