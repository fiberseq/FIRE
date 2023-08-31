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


rule fiber_locations:
    input:
        bam=lambda wc: data.loc[wc.sm, "bam"],
        coverage=rules.genome_bedgraph.output.bg,
        fai=ancient(f"{ref}.fai"),
    output:
        bed="results/coverage/fiber-locations.bed.gz",
        shuffled="results/coverage/fiber-locations-shuffled.bed.gz",
    threads: 8
    conda:
        conda
    shell:
        """
        samtools view -@ {threads} -F 2308 -u {input.bam} \
            | bedtools bamtobed -i - \
            | bgzip -@ {threads} \
        > {output.bed}

        bedtools shuffle -chrom \
            -excl <(zcat {input.coverage} | '$4 == 0') \
            -i {output.bed} \
            -g {input.fai} \
            | bgzip -@ {threads} \
        > {output.shuffled}
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
        hap1=temp("temp/{sm}/hap1/chunks/{chunk}.bed"),
        hap2=temp("temp/{sm}/hap2/chunks/{chunk}.bed"),
        unk=temp("temp/{sm}/unk/chunks/{chunk}.bed"),
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
            -o {output.haps} \
            --haps {output.hap1} {output.hap2} {output.unk}
        """


rule sort_model:
    input:
        bed="temp/{sm}/{hp}/chunks/{chunk}.bed",
    output:
        bed=temp("temp/{sm}/{hp}/chunks/{chunk}.sorted.bed"),
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
        bed="results/{sm}/{hp}/acc.model.results.bed.gz",
    benchmark:
        "benchmarks/{sm}/{hp}/merge_model_results.tsv"
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
