#
# Applying the model
#
rule bed_chunks:
    input:
        ref=ref,
        fai=f"{ref}.fai",
    output:
        beds=temp(
            expand(
                "temp/{sm}/chunks/{chunk}.chunk.bed",
                chunk=chunks,
                allow_missing=True,
            )
        ),
    threads: 1
    conda:
        "../envs/fibertools.yaml"
    params:
        keep_chrs="|".join(get_chroms()),
    shell:
        """
        fibertools split \
          -g <(grep -Pw '{params.keep_chrs}' {input.fai} | sort -k1,1 -k2,2n -k3,3n -k4,4) \
          -o {output.beds}
        """


rule extract_and_split:
    input:
        bam=ancient(lambda wc: data.loc[wc.sm, "bam"]),
        bed="temp/{sm}/chunks/{chunk}.chunk.bed",
    output:
        bed=temp("temp/{sm}/{chunk}.extract.all.bed.gz"),
    threads: 4
    resources:
        mem_mb=get_mem_mb,
    conda:
        conda
    shell:
        """
        samtools view \
          -F 2308 -b -M \
          -L {input.bed} \
          -@ {threads} {input.bam} \
          | ft -t {threads} extract --all - \
          | bgzip -@ {threads} \
          > {output.bed}
        """


def get_model(wc):
    if "model" in config:
        return ancient(config["model"])
    return ancient(rules.make_model.output.model)


rule apply_model:
    input:
        bed=rules.extract_and_split.output.bed,
        model=get_model,
    output:
        haps=temp("temp/{sm}/all/chunks/{chunk}.bed"),
    benchmark:
        "benchmarks/{sm}/chunks/apply_model_{chunk}.tsv"
    threads: 1
    resources:
        mem_mb=get_mem_mb,
    conda:
        "../envs/fibertools.yaml"
    priority: 0
    shell:
        """
        fibertools -v model -m {input.model} {input.bed} \
            -o {output.haps} 
        """


rule sort_model:
    input:
        bed=rules.apply_model.output.haps,
    output:
        bed=temp("temp/{sm}/chunks/{chunk}.sorted.bed"),
    threads: 4
    conda:
        conda
    resources:
        mem_mb=get_mem_mb_small,
    priority: 10
    shell:
        """
        LC_ALL=C sort \
            --parallel={threads} \
            -k1,1 -k2,2n -k3,3n -k4,4 \
            {input.bed} -o {output.bed}
        """


rule merge_model_results:
    input:
        beds=expand(rules.sort_model.output.bed, chunk=chunks, allow_missing=True),
    output:
        bed="results/{sm}/fiber-calls/model.results.bed.gz",
    benchmark:
        "benchmarks/{sm}/merge_model_results.tsv"
    threads: 8
    conda:
        conda
    params:
        n_chunks=n_chunks + 10,
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


rule element_coverages_by_type:
    input:
        bed=rules.merge_model_results.output.bed,
        fai=f"{ref}.fai",
    output:
        bed=temp("temp/{sm}/coverage/{hp}/{el_type}_coverage_{hp}.bed.gz"),
    benchmark:
        "benchmarks/{sm}/element_coverages/{el_type}_{hp}.tsv"
    conda:
        conda
    params:
        filter_cmd=grep_command_for_el_type,
        filter_hap=hap_grep_term,
    threads: 4
    shell:
        """
        ( \
            printf "#chrom\\tstart\\tend\\tcoverage\\n"; \
            bgzip -cd -@{threads} {input.bed} \
                | (rg -w {params.filter_hap} || true) \
                | {params.filter_cmd} \
                | bedtools genomecov -bg -i - -g {input.fai} \
        ) \
            | bgzip -@{threads} \
            > {output.bed}
        """


rule element_coverages:
    input:
        beds=expand(
            rules.element_coverages_by_type.output.bed,
            el_type=el_types,
            allow_missing=True,
        ),
    output:
        bed="results/{sm}/coverage/{hp}_element_coverages.bed.gz",
        tbi="results/{sm}/coverage/{hp}_element_coverages.bed.gz.tbi",
    conda:
        conda
    params:
        names=" ".join(el_types),
    threads: 4
    shell:
        """
        HAS_LINES=$(zcat {input.beds} | grep -cv '^#')
        if [ $HAS_LINES -eq 0 ]; then
            echo "No element coverages found for {wildcards.sm} {wildcards.hp}"
            touch {output.bed}
        else
            bedtools unionbedg -header -i {input.beds} -names {params.names} \
                | sed 's/^chrom/#chrom/' \
                | bgzip -@ {threads} \
            > {output.bed}
        fi
        tabix -p bed {output.bed}
        """
