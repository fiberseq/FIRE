rule filtered_and_shuffled_fiber_locations_chromosome:
    input:
        filtered=rules.fiber_locations.output.filtered,
        exclude=rules.exclude_from_shuffle.output.bed,
        fai=ancient(FAI),
    output:
        shuffled=temp("temp/{sm}/coverage/{chrom}.fiber-locations-shuffled.bed.gz"),
    threads: 4
    conda:
        DEFAULT_ENV
    shell:
        """
        tabix -h {input.filtered} {wildcards.chrom} \
            | bedtools shuffle -chrom -seed 42 \
                -excl {input.exclude} \
                -i - \
                -g {input.fai} \
            |  sort -k1,1 -k2,2n -k3,3n -k4,4 \
            | bgzip -@ {threads} \
        > {output.shuffled}
        """


rule filtered_and_shuffled_fiber_locations:
    input:
        shuffled=expand(
            rules.filtered_and_shuffled_fiber_locations_chromosome.output.shuffled,
            chrom=get_chroms(),
            allow_missing=True,
        ),
    output:
        shuffled="results/{sm}/coverage/filtered-for-coverage/fiber-locations-shuffled.bed.gz",
    threads: 1
    conda:
        DEFAULT_ENV
    shell:
        """
        cat {input.shuffled} > {output.shuffled}
        """


#
# FIRE sites and FDR tracks
#
rule fdr_table:
    input:
        fire=rules.fire_sites.output.bed,
        fiber_locations=rules.fiber_locations.output.filtered,
        shuffled=rules.filtered_and_shuffled_fiber_locations.output.shuffled,
        fai=ancient(FAI),
    output:
        tbl="results/{sm}/FDR-peaks/FIRE.score.to.FDR.tbl",
    benchmark:
        "results/{sm}/benchmarks/fdr_table/{sm}.txt"
    threads: 8
    conda:
        "../envs/python.yaml"
    params:
        script=workflow.source_path("../scripts/fire-null-distribution.py"),
    resources:
        mem_mb=get_mem_mb,
    shell:
        """
        python {params.script} -v 1 {input.fire} {input.fiber_locations} {input.fai} -s {input.shuffled} -o {output.tbl}
        """


rule fdr_track_chromosome:
    input:
        fire=rules.fire_sites.output.bed,
        fire_tbi=rules.fire_sites_index.output.tbi,
        fiber_locations=rules.fiber_locations.output.bed,
        fai=ancient(FAI),
        fdr_tbl=rules.fdr_table.output.tbl,
    output:
        fire=temp("temp/{sm}/FDR-peaks/{chrom}-fire.bed"),
        fiber=temp("temp/{sm}/FDR-peaks/{chrom}-fiber.bed"),
        bed=temp("temp/{sm}/FDR-peaks/{chrom}-FDR.track.bed"),
    threads: 8
    conda:
        "../envs/python.yaml"
    params:
        script=workflow.source_path("../scripts/fire-null-distribution.py"),
    resources:
        mem_mb=get_mem_mb_xl,
    shell:
        """
        tabix -h {input.fire} {wildcards.chrom} > {output.fire}
        tabix -h {input.fiber_locations} {wildcards.chrom} > {output.fiber}

        # check if file is empty
        if [ ! -s {output.fire} ]; then
            echo "No FIRE sites for {wildcards.chrom}"
            touch {output}
            exit 0
        fi

        python {params.script} -v 1 \
            {output.fire} {output.fiber} \
            {input.fai} -f {input.fdr_tbl} \
            -o {output.bed}
        """


rule fdr_track:
    input:
        beds=expand(
            rules.fdr_track_chromosome.output.bed,
            chrom=get_chroms(),
            allow_missing=True,
        ),
    output:
        fofn=temp("temp/{sm}/FDR-peaks/FDR.track.fofn"),
        bed="results/{sm}/FDR-peaks/FDR.track.bed.gz",
        tbi="results/{sm}/FDR-peaks/FDR.track.bed.gz.tbi",
    threads: 8
    conda:
        DEFAULT_ENV
    shell:
        """
        printf '\nMaking FOFN\n'
        echo {input.beds} > {output.fofn}
        
        printf '\nMake header\n'
        ((cat $(cat {output.fofn}) | grep "^#" | head -n 1) || true) \
            | bgzip -@ {threads} \
            > {output.bed}

        printf '\nConcatenating\n'
        cat $(cat {output.fofn}) \
            | grep -v "^#" \
            | bgzip -@ {threads} \
        >> {output.bed}

        printf '\nIndexing\n'
        tabix -f -p bed {output.bed}
        """


rule fdr_track_filtered:
    input:
        bed=rules.fdr_track.output.bed,
        minimum=rules.coverage.output.minimum,
        maximum=rules.coverage.output.maximum,
    output:
        bed="results/{sm}/FDR-peaks/FDR.track.coverage.filtered.bed.gz",
        tbi="results/{sm}/FDR-peaks/FDR.track.coverage.filtered.bed.gz.tbi",
    threads: 8
    conda:
        DEFAULT_ENV
    shell:
        """
        MIN=$(cat {input.minimum})
        MAX=$(cat {input.maximum})

        ( \
            zcat {input.bed} | head -n 1 || true; \
            zcat {input.bed} | bioawk -tc hdr -v MAX=$MAX -v MIN=$MIN  '$coverage > MIN && $coverage < MAX' \
        ) \
            | bgzip -@ {threads} \
            > {output.bed}
        tabix -f -p bed {output.bed}
        """


rule helper_fdr_peaks_by_fire_elements:
    input:
        bed=rules.fdr_track.output.bed,
        tbi=rules.fdr_track.output.tbi,
        fire=rules.fire_sites.output.bed,
        fire_tbi=rules.fire_sites_index.output.tbi,
    output:
        bed=temp("temp/{sm}/FDR-peaks/{chrom}-FDR-FIRE-peaks.bed.gz"),
    threads: 2
    conda:
        DEFAULT_ENV
    params:
        max_peak_fdr=MAX_PEAK_FDR,
        min_per_acc_peak=MIN_PER_ACC_PEAK,
    shell:
        """
        HEADER=$(zcat {input.bed} | head -n 1 || true)
        NC=$(echo $HEADER | awk '{{print NF}}' || true)
        FIRE_CT=$((NC+1))
        FIRE_ST=$((NC+2))
        FIRE_EN=$((NC+3))
        FIRE_SIZE=$((NC+4))
        FIRE_ID=$((NC+5))

        OUT_HEADER=$(printf "$HEADER\\tpeak_chrom\\tpeak_start\\tpeak_end\\tFIRE_IDs\\tFIRE_size_mean\\tFIRE_size_ssd\\tFIRE_start_ssd\\tFIRE_end_ssd")
        echo $OUT_HEADER
        
        ( \
            printf "$OUT_HEADER\\n"; \
            tabix -h {input.bed} {wildcards.chrom} \
                | rg -w "#chrom|True" \
                | csvtk filter -tT -C '$' -f "FDR<={params.max_peak_fdr}" \
                | csvtk filter -tT -C '$' -f "fire_coverage>1" \
                | bioawk -tc hdr 'NR==1 || ($fire_coverage/$coverage>={params.min_per_acc_peak})' \
                | bedtools intersect -wa -wb -sorted -a - \
                    -b <(tabix {input.fire} {wildcards.chrom} \
                            | cut -f 1-3 \
                            | awk -v OFMT="%f" '{{print $0"\t"$3-$2"\t"NR}}' \
                        ) \
                | bedtools groupby -g 1-$NC \
                    -o first,median,median,collapse,mean,sstdev,sstdev,sstdev \
                    -c $FIRE_CT,$FIRE_ST,$FIRE_EN,$FIRE_ID,$FIRE_SIZE,$FIRE_SIZE,$FIRE_ST,$FIRE_EN \
        ) \
            | hck -f 1,$FIRE_ST,$FIRE_EN,2-$NC,$FIRE_SIZE- \
            | csvtk round -tT -C '$' -n 0 -f 2,3 \
            | bedtools sort -header -i - \
            | bgzip -@ {threads} \
            > {output.bed}
        """


rule fdr_peaks_by_fire_elements_chromosome:
    input:
        bed=rules.helper_fdr_peaks_by_fire_elements.output.bed,
        minimum=rules.coverage.output.minimum,
        maximum=rules.coverage.output.maximum,
    output:
        bed=temp("temp/{sm}/FDR-peaks/grouped-{chrom}-FDR-FIRE-peaks.bed.gz"),
    threads: 8
    conda:
        "../envs/python.yaml"
    params:
        script=workflow.source_path("../scripts/merge_fire_peaks.py"),
        min_frac_accessible=MIN_FRAC_ACCESSIBLE,
    shell:
        """
        zcat {input.bed} \
            | python {params.script} -v 1 \
                --max-cov $(cat {input.maximum}) \
                --min-cov $(cat {input.minimum}) \
                --min-frac-accessible {params.min_frac_accessible} \
            | bgzip -@ {threads} \
        > {output.bed}
        """


rule fdr_peaks_by_fire_elements:
    input:
        beds=expand(
            rules.fdr_peaks_by_fire_elements_chromosome.output.bed,
            chrom=get_chroms(),
            allow_missing=True,
        ),
    output:
        fofn=temp("temp/{sm}/FDR-peaks/FDR-FIRE-peaks.fofn"),
        bed="results/{sm}/FDR-peaks/FDR-FIRE-peaks.bed.gz",
        tbi="results/{sm}/FDR-peaks/FDR-FIRE-peaks.bed.gz.tbi",
    threads: 8
    conda:
        DEFAULT_ENV
    shell:
        """
        printf "\nMaking FOFN\n"
        echo {input.beds} > {output.fofn}

        printf "\nMaking header\n"        
        ((cat $(cat {output.fofn}) | bgzip -cd | grep "^#" | head -n 1) || true) \
            | bgzip -@ {threads} > {output.bed}

        printf "\nConcatenating\n"
        cat $(cat {output.fofn}) | bgzip -cd -@ {threads} | grep -v "^#" \
            | bgzip -@ {threads} >> {output.bed}
        
        printf "\nIndexing\n"
        tabix -f -p bed {output.bed}
        """


rule wide_fdr_peaks:
    input:
        bed=rules.fdr_peaks_by_fire_elements.output.bed,
        track=rules.fdr_track.output.bed,
        fai=ancient(FAI),
    output:
        bed="results/{sm}/FDR-peaks/FDR-wide-peaks.bed.gz",
        tbi="results/{sm}/FDR-peaks/FDR-wide-peaks.bed.gz.tbi",
        bb="results/{sm}/trackHub/bb/FDR-wide-peaks.bb",
    conda:
        DEFAULT_ENV
    threads: 4
    params:
        nuc_size=config.get("nucleosome_size", 147),
        max_peak_fdr=MAX_PEAK_FDR,
        min_frac_acc=max(MIN_FRAC_ACCESSIBLE, MIN_PER_ACC_PEAK),
    shell:
        """
        FILE={output.bed}
        TMP="${{FILE%.*}}"
        echo $TMP

        ( \
            zcat {input.bed}; \
            bioawk -tc hdr 'NR==1 || $FDR<={params.max_peak_fdr}' {input.track} \
                | bioawk -tc hdr 'NR==1 || $coverage>0' \
                | bioawk -tc hdr 'NR==1 || ($fire_coverage/$coverage>={params.min_frac_acc})' \
        ) \
            | cut -f 1-3 \
            | bedtools sort \
            | bedtools merge -d {params.nuc_size} \
        > $TMP
        bedToBigBed $TMP {input.fai} {output.bb}        
        bgzip -f -@ {threads} $TMP
        tabix -p bed {output.bed}
        """


rule one_percent_fdr_peaks:
    input:
        bed=rules.fdr_peaks_by_fire_elements.output.bed,
        track=rules.fdr_track.output.bed,
    output:
        bed="results/{sm}/FDR-peaks/one-percent-FDR/FDR-01-FIRE-peaks.bed.gz",
        tbi="results/{sm}/FDR-peaks/one-percent-FDR/FDR-01-FIRE-peaks.bed.gz.tbi",
        wide="results/{sm}/FDR-peaks/one-percent-FDR/FDR-01-FIRE-wide-peaks.bed.gz",
        wide_tbi="results/{sm}/FDR-peaks/one-percent-FDR/FDR-01-FIRE-wide-peaks.bed.gz.tbi",
    threads: 4
    conda:
        DEFAULT_ENV
    params:
        nuc_size=config.get("nucleosome_size", 147),
    shell:
        """
        zcat {input.bed} \
            | csvtk filter -tT -C '$' -f "FDR<=0.01" \
            | bgzip -@ {threads} \
            > {output.bed}
        tabix -f -p bed {output.bed}

        ( \
            zcat {output.bed}; \
            bioawk -tc hdr '$FDR<=0.01' {input.track} \
        ) \
            | cut -f 1-3 \
            | bedtools sort \
            | bedtools merge -d {params.nuc_size} \
            | bgzip -@ {threads} \
        > {output.wide}
        tabix -f -p bed {output.wide}
        """


rule peaks_vs_percent:
    input:
        bed=rules.fdr_peaks_by_fire_elements.output.bed,
    output:
        fig1=report(
            "results/{sm}/FDR-peaks/{sm}.peaks-vs-percent.pdf",
            category="Peak calls",
        ),
    threads: 8
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/peaks-vs-percent.R"
