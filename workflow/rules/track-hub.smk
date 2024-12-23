rule percent_accessible:
    input:
        bed=rules.pileup.output.bed,
        fai=ancient(FAI),
    output:
        tmp=temp("temp/{sm}/{hp}/{v}-percent.accessible.bed"),
        bw="results/{sm}/trackHub-{v}/bw/{hp}.percent.accessible.bw",
    threads: 4
    conda:
        DEFAULT_ENV
    resources:
        mem_mb=get_mem_mb,
    params:
        suffix=get_hap_col_suffix,
        nzooms=NZOOMS,
        chrom=get_chroms()[0],
    shell:
        """
        bgzip -cd {input.bed} \
            | bioawk -tc hdr '$coverage{params.suffix}>0' \
            | bioawk -tc hdr \
              'NR>1{{print $1,$2,$3,100*$fire_coverage{params.suffix}/$coverage{params.suffix}}}' \
        > {output.tmp}

        # add fake if file is empty
        if [[ -s {output.tmp} ]]; then
            echo "File is not empty"
        else
            echo "File is empty"
            printf "{params.chrom}\t0\t1\t0\\n" > {output.tmp}
        fi


        bigtools bedgraphtobigwig \
            --nzooms {params.nzooms} -s start \
            {output.tmp} {input.fai} {output.bw}
        """


rule element_coverages_bw:
    input:
        bed=rules.pileup.output.bed,
        fai=ancient(FAI),
    output:
        bw="results/{sm}/trackHub-{v}/bw/{hp}.{el_type}.coverage.bw",
    conda:
        DEFAULT_ENV
    params:
        nzooms=NZOOMS,
        cut_cmd=pileup_cut_cmd,
    shell:
        """
        bgzip -cd {input.bed} \
            | {params.cut_cmd} \
            | grep -v "^#" \
            | bigtools bedgraphtobigwig \
                -s start --nzooms {params.nzooms} \
                - {input.fai} {output.bw}
        """


rule fdr_track_to_bw:
    input:
        bed=rules.pileup.output.bed,
        fai=ancient(FAI),
    output:
        bw="results/{sm}/trackHub-{v}/bw/{col}.bw",
    threads: 4
    conda:
        DEFAULT_ENV
    params:
        nzooms=NZOOMS,
    shell:
        """
        hck -z -f 1-3 -F {wildcards.col} {input.bed} \
            | grep -v "^#" \
            | bigtools bedgraphtobigwig \
                -s start --nzooms {params.nzooms} \
                - {input.fai} {output.bw}
        """


rule fire_peaks_bb:
    input:
        bed=rules.fire_peaks.output.bed,
        fai=ancient(FAI),
    output:
        bb="results/{sm}/trackHub-{v}/bb/fire-peaks.bb",
    threads: 4
    conda:
        DEFAULT_ENV
    params:
        bedfmt=workflow.source_path("../templates/fire_peak.as"),
    shell:
        """
        bgzip -cd {input.bed} \
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
        bb="results/{sm}/trackHub-{v}/bb/hap_differences.bb",
    threads: 4
    resources:
        mem_mb=get_mem_mb,
    conda:
        DEFAULT_ENV
    params:
        chrom=get_chroms()[0],
        bed9_as=workflow.source_path("../templates/bed9.as"),
    shell:
        """
        ( \
            printf "{params.chrom}\t0\t1\tfake\t100\t+\t0\t1\t230,230,230\\n"; \
            bedtools sort -i {input.bed9} \
        ) \
            | bigtools bedtobigbed \
                -s start -a {params.bed9_as} \
                - {input.fai} {output.bb}
        """


rule trackhub:
    input:
        cov=rules.coverage.output.cov,
    output:
        hub="results/{sm}/trackHub-{v}/hub.txt",
        description="results/{sm}/trackHub-{v}/fire-description.html",
    resources:
        load=get_load,
    threads: 4
    conda:
        "../envs/python.yaml"
    params:
        ref=REF_NAME,
        script=workflow.source_path("../scripts/trackhub.py"),
        description=workflow.source_path("../templates/fire-description.html"),
    shell:
        """
        python {params.script} -v 2 \
          --trackhub-dir results/{wildcards.sm}/trackHub-{wildcards.v} \
          --reference {params.ref} \
          --sample {wildcards.sm} \
          --average-coverage $(cat {input.cov}) 
        cp {params.description} {output.description}
        """
