#
# fiber locations and coverages
#
rule fiber_locations_chromosome:
    input:
        bam=lambda wc: data.loc[wc.sm, "bam"],
    output:
        bed=temp("temp/{sm}/coverage/{chrom}.fiber-locations.bed.gz"),
    threads: 8
    conda:
        conda
    shell:
        """
        # get fiber locations
        (samtools view -@ {threads} -F 2308 -u {input.bam} {wildcards.chrom} \
            | ft extract -t {threads} -s --all - \
            | hck -F '#ct' -F st -F en -F fiber -F strand -F HP ) \
            | (grep -v "^#" || true) \
            | bgzip -@ {threads} \
        > {output.bed}
        """


rule fiber_locations:
    input:
        fibers=expand(
            rules.fiber_locations_chromosome.output.bed,
            chrom=get_chroms(),
            allow_missing=True,
        ),
    output:
        bed="results/{sm}/coverage/fiber-locations.bed.gz",
        bed_tbi="results/{sm}/coverage/fiber-locations.bed.gz.tbi",
    threads: 1
    conda:
        conda
    shell:
        """
        cat {input.fibers} > {output.bed}
        tabix -f -p bed {output.bed}
        """


rule filtered_and_shuffled_fiber_locations_chromosome:
    input:
        bed=rules.fiber_locations_chromosome.output.bed,
        fai=ancient(f"{ref}.fai"),
        # required for the coverage function to work
        bg=rules.genome_bedgraph.output.bg,
    output:
        bed=temp("temp/{sm}/coverage/{chrom}.fiber-locations-filtered.bed.gz"),
        bg=temp("temp/{sm}/coverage/{chrom}.fiber-locations-filtered.coverage.bed.gz"),
        shuffled=temp("temp/{sm}/coverage/{chrom}.fiber-locations-shuffled.bed.gz"),
    threads: 4
    params:
        min_cov=get_min_coverage,
        max_cov=get_max_coverage,
    conda:
        conda
    shell:
        """
        # check if file is empty
        if [ -n "$(gunzip <{input.bed} | head -c 1 | tr '\\0\\n' __)" ]; then
            echo "input is not empty"
        else
            echo "input is empty"
            bgzip -c <(printf "") > {output.bed}
            bgzip -c <(printf "") > {output.bg}
            bgzip -c <(printf "") > {output.shuffled}
            exit 0
        fi

        # get fiber locations
        bedtools intersect -v -f 0.2 \
            -a {input.bed} \
            -b <(zcat {input.bg} | awk '$4 <= {params.min_cov} || $4 >= {params.max_cov}') \
        | bgzip -@ {threads} \
        > {output.bed}

        # get bedgraph
        bedtools genomecov -bga -i {output.bed} -g {input.fai} | bgzip -@ {threads} > {output.bg}

        # make shuffled fiber locations
        bedtools shuffle -chrom \
            -excl <(zcat {output.bg} | awk '$4 == 0') \
            -i {output.bed} \
            -g {input.fai} \
            |  sort -k1,1 -k2,2n -k3,3n -k4,4 \
            | bgzip -@ {threads} \
        > {output.shuffled}
        """


rule filtered_and_shuffled_fiber_locations:
    input:
        bed=expand(
            rules.filtered_and_shuffled_fiber_locations_chromosome.output.bed,
            chrom=get_chroms(),
            allow_missing=True,
        ),
        shuffled=expand(
            rules.filtered_and_shuffled_fiber_locations_chromosome.output.shuffled,
            chrom=get_chroms(),
            allow_missing=True,
        ),
    output:
        bed="results/{sm}/FDR-peaks/filtered-for-fdr/fiber-locations.bed.gz",
        shuffled="results/{sm}/FDR-peaks/filtered-for-fdr/fiber-locations-shuffled.bed.gz",
    threads: 1
    conda:
        conda
    shell:
        """
        cat {input.bed} > {output.bed}
        cat {input.shuffled} > {output.shuffled}
        """


#
# FIRE sites and FDR tracks
#
rule fdr_table:
    input:
        fire=rules.fire_sites.output.bed,
        fiber_locations=rules.filtered_and_shuffled_fiber_locations.output.bed,
        shuffled=rules.filtered_and_shuffled_fiber_locations.output.shuffled,
        fai=ancient(f"{ref}.fai"),
    output:
        tbl="results/{sm}/FDR-peaks/FIRE.score.to.FDR.tbl",
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


rule fdr_track:
    input:
        fire=rules.fire_sites.output.bed,
        fiber_locations=rules.fiber_locations.output.bed,
        fai=ancient(f"{ref}.fai"),
        fdr_tbl=rules.fdr_table.output.tbl,
    output:
        bed="results/{sm}/FDR-peaks/FDR.track.bed.gz",
    threads: 8
    conda:
        "../envs/python.yaml"
    params:
        script=workflow.source_path("../scripts/fire-null-distribution.py"),
    resources:
        mem_mb=get_mem_mb,
    shell:
        """
        OUT={output.bed}
        TMP_OUT="${{OUT%.*}}"
        echo $TMP_OUT
        python {params.script} \
            -v 1 {input.fire} {input.fiber_locations} {input.fai} -f {input.fdr_tbl} \
            -o $TMP_OUT
        
        bgzip -@ {threads} $TMP_OUT
        """


rule fdr_track_filtered:
    input:
        bed=rules.fdr_track.output.bed,
    output:
        bed="results/{sm}/FDR-peaks/FDR.track.coverage.filtered.bed.gz",
        tbi="results/{sm}/FDR-peaks/FDR.track.coverage.filtered.bed.gz.tbi",
    threads: 8
    conda:
        conda
    params:
        min_cov=get_min_coverage,
        max_cov=get_max_coverage,
    shell:
        """
        zcat {input.bed} \
            | csvtk -tT filter -f '$5 >= {params.min_cov} && $5 <= {params.max_cov}' \
            | bgzip -@ {threads} \
            > {output.bed}
        tabix -f -p bed {output.bed}
        
        #| awk '$5 >= {params.min_cov} && $5 <= {params.max_cov}' \
        """
