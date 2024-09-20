rule filtered_and_shuffled_fiber_locations_chromosome:
    input:
        filtered=rules.fiber_locations.output.filtered,
        exclude=rules.exclude_from_shuffle.output.bed,
        fai=ancient(FAI),
    output:
        shuffled=temp("temp/{sm}/shuffle/{chrom}.fiber-locations-shuffled.bed.gz"),
    threads: 4
    conda:
        DEFAULT_ENV
    shell:
        """
        tabix {input.filtered} {wildcards.chrom} \
            | bioawk -t '{{print $1,$2,$3,$4,$2}}' \
            | bedtools shuffle -chrom -seed 42 \
                -excl {input.exclude} \
                -i - \
                -g {input.fai} \
            |  sort -k1,1 -k2,2n -k3,3n -k4,4 \
            | bgzip -@ {threads} \
        > {output.shuffled}
        """


rule shuffled_pileup_chromosome:
    input:
        cram=rules.merged_fire_bam.output.cram,
        shuffled=rules.filtered_and_shuffled_fiber_locations_chromosome.output.shuffled,
    output:
        bed=temp("temp/{sm}/shuffle/{chrom}.pileup.bed.gz"),
    threads: 4
    conda:
        DEFAULT_ENV
    shell:
        """
        {FT_EXE} pileup {input.cram} {wildcards.chrom} -t {threads} \
            --fiber-coverage --shuffle {input.shuffled} \
            --no-msp --no-nuc \
            | bgzip -@ {threads} \
        > {output.bed}    
        """


rule shuffled_pileup:
    input:
        beds=expand(
            rules.shuffled_pileup_chromosome.output.bed,
            chrom=get_chroms(),
            allow_missing=True,
        ),
    output:
        bed=temp("temp/{sm}/shuffle/shuffled-pileup.bed.gz"),
        tbi=temp("temp/{sm}/shuffle/shuffled-pileup.bed.gz.tbi"),
    threads: 4
    conda:
        DEFAULT_ENV
    shell:
        """
        cat {input.beds} > {output.bed}
        tabix -p bed {output.bed}
        """


#
# FIRE sites and FDR tracks
#
rule fdr_table:
    input:
        shuffled=rules.shuffled_pileup.output.bed,
    output:
        tbl="results/{sm}/fire-peaks/{sm}-fire-score-to-fdr.tbl",
    conda:
        "../envs/python.yaml"
    params:
        script=workflow.source_path("../scripts/fdr-table.py"),
    resources:
        mem_mb=get_mem_mb,
    shell:
        """
        python {params.script} -v 1 {input.shuffled} {output.tbl}
        """


# Colnames made by this
# #chrom  start   end
# coverage  fire_coverage   score   nuc_coverage    msp_coverage
# coverage_H1  fire_coverage_H1  score_H1   nuc_coverage_H1 msp_coverage_H1
# coverage_H2  fire_coverage_H2  score_H2   nuc_coverage_H2 msp_coverage_H2
rule pileup_chromosome:
    input:
        bam=rules.merged_fire_bam.output.cram,
    output:
        bed=temp("temp/{sm}/{chrom}.pileup.bed.gz"),
    threads: 12
    conda:
        DEFAULT_ENV
    shell:
        """
        {FT_EXE} pileup -t {threads} \
            --haps --fiber-coverage \
            {input.bam} {wildcards.chrom} \
            | bgzip -@ {threads} \
            > {output.bed}
        """


rule fdr_track_chromosome:
    input:
        pileup=rules.pileup_chromosome.output.bed,
        fdr_tbl=rules.fdr_table.output.tbl,
    output:
        bed=temp("temp/{sm}/fire-peaks/{chrom}-FDR.track.bed"),
    threads: 8
    conda:
        "../envs/python.yaml"
    params:
        script=workflow.source_path("../scripts/fdr-table.py"),
    resources:
        mem_mb=get_mem_mb_xl,
    shell:
        """
        python {params.script} -v 1 \
            --fdr-table {input.fdr_tbl} \
            {input.pileup} {output.bed}
        """


rule pileup:
    input:
        beds=expand(
            rules.fdr_track_chromosome.output.bed,
            chrom=get_chroms(),
            allow_missing=True,
        ),
    output:
        fofn=temp("temp/{sm}/fire/fire-pileup.fofn"),
        bed="results/{sm}/{sm}-fire-pileup.bed.gz",
        tbi="results/{sm}/{sm}-fire-pileup.bed.gz.tbi",
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


rule helper_fdr_peaks_by_fire_elements:
    input:
        bed=rules.pileup.output.bed,
        tbi=rules.pileup.output.tbi,
        fire=rules.fire_sites.output.bed,
        fire_tbi=rules.fire_sites_index.output.tbi,
    output:
        bed=temp("temp/{sm}/fire-peaks/{chrom}-fire-peaks.bed.gz"),
    threads: 2
    conda:
        DEFAULT_ENV
    params:
        max_peak_fdr=MAX_PEAK_FDR,
        min_per_acc_peak=MIN_PER_ACC_PEAK,
    shell:
        """
        HEADER=$(bgzip -cd {input.bed} | head -n 1 || true)
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
                | bioawk -tc hdr '(NR==1)||($is_local_max=="true")' \
                | csvtk filter -tT -C '$' -f "FDR<={params.max_peak_fdr}" \
                | csvtk filter -tT -C '$' -f "fire_coverage>1" \
                | bioawk -tc hdr '(NR==1)||($fire_coverage/$coverage>={params.min_per_acc_peak})' \
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
        bed=temp("temp/{sm}/fire-peaks/grouped-{chrom}-fire-peaks.bed.gz"),
    threads: 8
    conda:
        "../envs/python.yaml"
    params:
        script=workflow.source_path("../scripts/merge_fire_peaks.py"),
        min_frac_accessible=MIN_FRAC_ACCESSIBLE,
    shell:
        """
        bgzip -cd {input.bed} \
            | python {params.script} -v 1 \
                --max-cov $(cat {input.maximum}) \
                --min-cov $(cat {input.minimum}) \
                --min-frac-accessible {params.min_frac_accessible} \
            | bgzip -@ {threads} \
        > {output.bed}
        """


rule fire_peaks:
    input:
        beds=expand(
            rules.fdr_peaks_by_fire_elements_chromosome.output.bed,
            chrom=get_chroms(),
            allow_missing=True,
        ),
    output:
        fofn=temp("temp/{sm}/fire-peaks/{sm}-fire-peaks.fofn"),
        bed="results/{sm}/{sm}-fire-peaks.bed.gz",
        tbi="results/{sm}/{sm}-fire-peaks.bed.gz.tbi",
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


rule wide_fire_peaks:
    input:
        bed=rules.fire_peaks.output.bed,
        track=rules.pileup.output.bed,
        fai=ancient(FAI),
    output:
        bed="results/{sm}/fire-peaks/{sm}-fire-wide-peaks.bed.gz",
        tbi="results/{sm}/fire-peaks/{sm}-fire-wide-peaks.bed.gz.tbi",
        bb="results/{sm}/trackHub/bb/fire-wide-peaks.bb",
    conda:
        DEFAULT_ENV
    threads: 4
    params:
        nuc_size=config.get("nucleosome_size", 147),
        max_peak_fdr=MAX_PEAK_FDR,
        min_frac_acc=max(MIN_FRAC_ACCESSIBLE, MIN_PER_ACC_PEAK),
        bed3_as=workflow.source_path("../templates/bed3.as"),
    shell:
        """
        ( \
            bgzip -cd {input.bed}; \
            bioawk -tc hdr 'NR==1 || $FDR<={params.max_peak_fdr}' {input.track} \
                | bioawk -tc hdr 'NR==1 || $coverage>0' \
                | bioawk -tc hdr 'NR==1 || ($fire_coverage/$coverage>={params.min_frac_acc})' \
        ) \
            | cut -f 1-3 \
            | bedtools sort \
            | bedtools merge -d {params.nuc_size} \
            | bgzip -@ {threads} \
        > {output.bed}
        
        bgzip -cd -@ 16 {output.bed} \
            | bigtools bedtobigbed \
                -s start -a {params.bed3_as} \
                - {input.fai} {output.bb}        
        
        tabix -p bed {output.bed}
        """


rule one_percent_fire_peaks:
    input:
        bed=rules.fire_peaks.output.bed,
        track=rules.pileup.output.bed,
    output:
        bed="results/{sm}/fire-peaks/one-percent-FDR/{sm}-01-fire-peaks.bed.gz",
        tbi="results/{sm}/fire-peaks/one-percent-FDR/{sm}-01-fire-peaks.bed.gz.tbi",
        wide="results/{sm}/fire-peaks/one-percent-FDR/{sm}-01-fire-wide-peaks.bed.gz",
        wtbi="results/{sm}/fire-peaks/one-percent-FDR/{sm}-01-fire-wide-peaks.bed.gz.tbi",
    threads: 4
    conda:
        DEFAULT_ENV
    params:
        nuc_size=config.get("nucleosome_size", 147),
    shell:
        """
        bgzip -cd {input.bed} \
            | csvtk filter -tT -C '$' -f "FDR<=0.01" \
            | bgzip -@ {threads} \
            > {output.bed}
        tabix -f -p bed {output.bed}

        ( \
            bgzip -cd {output.bed}; \
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
        bed=rules.fire_peaks.output.bed,
    output:
        fig1=report(
            "results/{sm}/fire-peaks/{sm}.peaks-vs-percent.pdf",
            category="Peak calls",
        ),
    threads: 8
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/peaks-vs-percent.R"
