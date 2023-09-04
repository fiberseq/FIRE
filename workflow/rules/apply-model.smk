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
        conda
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
        #hap1=temp("temp/{sm}/hap1/chunks/{chunk}.bed"),
        #hap2=temp("temp/{sm}/hap2/chunks/{chunk}.bed"),
        #unk=temp("temp/{sm}/unk/chunks/{chunk}.bed"),
    benchmark:
        "benchmarks/{sm}/chunks/apply_model_{chunk}.tsv"
    threads: 1
    resources:
        mem_mb=get_mem_mb,
    conda:
        conda
    priority: 0
    shell:
        """
        fibertools -v model -m {input.model} {input.bed} \
            -o {output.haps} 
        """
        #--haps {output.hap1} {output.hap2} {output.unk}


rule sort_model:
    input:
        bed=rules.apply_model.output.haps,
        #bed="temp/{sm}/{hp}/chunks/{chunk}.bed",
    output:
        bed=temp("temp/{sm}/chunks/{chunk}.sorted.bed"),
        #bed=temp("temp/{sm}/{hp}/chunks/{chunk}.sorted.bed"),
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
            | awk '$10<={params.min_fdr}' \
            | bgzip -@{threads} \
            > {output.bed}
        """


rule nucleosome_and_linker_coverages:
    input:
        bed=rules.merge_model_results.output.bed,
        fai=f"{ref}.fai",
    output:
        bed="results/{sm}/fiber-calls/{el_type}_coverage_{hp}.bed.gz",
    threads: 8
    conda:
        conda
    params:
        get_color=lambda wc: "230,230,230" if wc.el_type == "nucleosome" else "147,112,219",
        hap_grep=lambda wc: "" if wc.hp == "all" else wc.hp,
    shell:
        """
        bgzip -cd -@{threads} {input.bed} \
            | grep -w "{params.hap_grep}" \
            | grep "^#\|{params.get_color}"  \
            | bedtools genomecov -bga -i - -g {input.fai} \
            | bgzip -@{threads} \
            > {output.bed}
        """
