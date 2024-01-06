#
# Applying the model
#
rule fire:
    input:
        bam=ancient(lambda wc: data.loc[wc.sm, "bam"]),
    output:
        bam=temp("temp/{sm}/fire/{chrom}.fire.bam"),
    threads: 8
    resources:
        mem_mb=8 * 1024,
    conda:
        "../envs/env.yaml"
    shell:
        """
        samtools view -u -@ {threads} {input.bam} {wildcards.chrom} \
            | ft fire -t {threads} --skip-no-m6a - {output.bam}
        """


rule merged_fire_bam:
    input:
        bams=expand(rules.fire.output.bam, chrom=get_chroms(), allow_missing=True),
    output:
        bam="results/{sm}/fire/{sm}.fire.bam",
        bai="results/{sm}/fire/{sm}.fire.bam.bai",
    threads: 16
    resources:
        mem_mb=8 * 1024,
    conda:
        "../envs/env.yaml"
    shell:
        """
        samtools merge -@ {threads} -o {output.bam} {input.bams}
        samtools index -@ {threads} {output.bam}
        """


rule extract_from_fire:
    input:
        bam=rules.fire.output.bam,
    output:
        bed=temp("temp/{sm}/chrom/{chrom}.sorted.bed"),
    threads: 4
    conda:
        conda
    resources:
        mem_mb=16 * 1024,
    priority: 10
    shell:
        """
        ft fire -t {threads} --extract {input.bam} \
            | LC_ALL=C sort \
                --parallel={threads} \
                -k1,1 -k2,2n -k3,3n -k4,4 \
            > {output.bed}
        """


rule merge_model_results:
    input:
        beds=expand(
            rules.extract_from_fire.output.bed, chrom=get_chroms(), allow_missing=True
        ),
    output:
        bed="results/{sm}/fiber-calls/model.results.bed.gz",
    benchmark:
        "benchmarks/{sm}/merge_model_results.tsv"
    threads: 8
    conda:
        conda
    params:
        n_chunks=len(get_chroms()) + 1,
    priority: 20
    shell:
        """
        LC_ALL=C sort \
            --parallel={threads} \
            --batch-size={params.n_chunks} \
            -k1,1 -k2,2n -k3,3n -k4,4 -m -u \
            {input.beds} \
          | bgzip -@ {threads} \
          > {output.bed}
        """


rule index_model_results:
    input:
        bed=rules.merge_model_results.output.bed,
    output:
        tbi=rules.merge_model_results.output.bed + ".tbi",
    conda:
        conda
    shell:
        """
        tabix -p bed {input.bed}
        """


rule fire_sites:
    input:
        bed=rules.merge_model_results.output.bed,
    output:
        bed="results/{sm}/fiber-calls/FIRE.bed.gz",
    threads: 8
    conda:
        conda
    params:
        min_fdr=min_fire_fdr,
    shell:
        """
        bgzip -cd -@{threads} {input.bed} \
            | bioawk -tc hdr '$10<={params.min_fdr}' \
            | bgzip -@{threads} \
            > {output.bed}
        """


rule fire_sites_index:
    input:
        bed=rules.fire_sites.output.bed,
    output:
        tbi="results/{sm}/fiber-calls/FIRE.bed.gz.tbi",
    threads: 1
    conda:
        conda
    shell:
        """
        tabix -p bed {input.bed}
        """


rule split_by_hap_per_chrom:
    input:
        bed=rules.merge_model_results.output.bed,
        tbi=rules.index_model_results.output.tbi,
        fai=f"{ref}.fai",
    output:
        #bed=temp("temp/{sm}/coverage/{hp}/{el_type}_coverage_{hp}_{chrom}.bed.gz"),
        both=pipe("temp/{sm}/coverage/all/{chrom}.bed"),
        H1=pipe("temp/{sm}/coverage/hap1/{chrom}.bed"),
        H2=pipe("temp/{sm}/coverage/hap2/{chrom}.bed"),
    conda:
        conda
    resources:
        disk_mb=100,
        time=240,
        mem_mb=4 * 1024,
    shell:
        """
        tabix {input.bed} {wildcards.chrom} | tee \
            >( (rg -w H1 || true) > {output.H1} ) \
            >( (rg -w H2 || true) > {output.H2} ) \
            > {output.both}
        """


# el_types = ["fire", "linker", "nucleosome"]
rule split_hap_by_element_type_per_chrom:
    input:
        bed="temp/{sm}/coverage/{hp}/{chrom}.bed",
        fai=f"{ref}.fai",
    output:
        fire=temp("temp/{sm}/coverage/{hp}/fire_{chrom}.bed.gz"),
        link=temp("temp/{sm}/coverage/{hp}/linker_{chrom}.bed.gz"),
        nuc=temp("temp/{sm}/coverage/{hp}/nucleosome_{chrom}.bed.gz"),
    params:
        min_fire_fdr=min_fire_fdr,
    threads: 2
    resources:
        disk_mb=100,
        mem_mb=8 * 1024,
    shell:
        """
        cat {input.bed} \
            | tee \
                >( \
                    awk '$10<={params.min_fire_fdr}' \
                        | bedtools genomecov -bg -i - -g {input.fai} \
                        | bgzip > {output.fire} \
                ) \
                >( \
                    (awk '$10<=1.0 && $10>{params.min_fire_fdr}') \
                        | bedtools genomecov -bg -i - -g {input.fai} \
                        | bgzip > {output.link} \
                ) \
            | awk '$10>1.0' \
            | bedtools genomecov -bg -i - -g {input.fai} \
            | bgzip > {output.nuc}
        """

rule element_coverages_per_chrom:
    input:
        beds=expand(
            "temp/{sm}/coverage/{hp}/{el_type}_{chrom}.bed.gz",
            el_type=el_types,
            allow_missing=True,
        ),
    output:
        bed=temp("temp/{sm}/coverage/{hp}_{chrom}_element_coverages.bed.gz"),
    conda:
        conda
    params:
        names="\t".join(el_types),
    resources:
        time=300,
    threads: 4
    shell:
        """
        HAS_LINES=$(zcat {input.beds} | grep -cv '^#') || true
        if [ $HAS_LINES -eq 0 ]; then
            echo "No element coverages found for {wildcards.sm} {wildcards.hp} {wildcards.chrom}"
            printf "#chrom\\tstart\\tend\\t{params.names}\\n{wildcards.chrom}\\t0\\t1\\t0\\t0\\t0\\n" \
                | bgzip -@{threads} \
                > {output.bed}
        else
            bedtools unionbedg -header -i {input.beds} -names {params.names} \
                | sed 's/^chrom/#chrom/' \
                | bgzip -@ {threads} \
            > {output.bed}
        fi
        """

rule element_coverages:
    input:
        beds=expand(
            rules.element_coverages_per_chrom.output.bed,
            chrom=get_chroms(),
            allow_missing=True,
        ),
    output:
        bed="results/{sm}/coverage/{hp}_element_coverages.bed.gz",
        tbi="results/{sm}/coverage/{hp}_element_coverages.bed.gz.tbi",
    conda:
        conda
    threads: 1
    shell:
        """
        cat {input.beds} > {output.bed}
        tabix -p bed {output.bed}
        """
