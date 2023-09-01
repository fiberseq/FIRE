#
# Coverage calculations
#
rule genome_bedgraph:
    input:
        bam=ancient(lambda wc: data.loc[wc.sm, "bam"]),
        fai=ancient(f"{ref}.fai"),
    output:
        d4="results/{sm}/coverage/{sm}.d4",
        bg="results/{sm}/coverage/{sm}.bed.gz",
        median="results/{sm}/coverage/{sm}.median.chromosome.coverage.bed",
    threads: 16
    conda:
        conda
    shell:
        """ 
        d4tools create -F 2308 -t {threads} -Azr {input.fai} {input.bam} {output.d4}
        d4tools view {output.d4} | bgzip -@ {threads} > {output.bg}
        zcat {output.bg} \
            | awk '$4>0' \
            | datamash -g 1 min 2 max 3 median 4 \
        > {output.median}
        """


rule average_coverage:
    input:
        median=rules.genome_bedgraph.output.median,
    output:
        cov="results/{sm}/coverage/{sm}.median.coverage.txt",
    run:
        find_median_coverage(input["median"], outfile=output["cov"])


#
# fiber locations and coverages
#
rule fiber_locations_chromosome:
    input:
        bam=lambda wc: data.loc[wc.sm, "bam"],
    output:
        bed=temp("temp/{sm}/coverage/{chrom}.fiber-locations.bed.gz"),
    threads: 4
    conda:
        conda
    shell:
        """
        # get fiber locations
        samtools view -@ {threads} -F 2308 -u {input.bam} {wildcards.chrom} \
            | bedtools bamtobed -i - \
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
        min_coverage=get_min_coverage,
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
            -b <(zcat {input.bg} | awk '$4 <= {params.min_coverage}') \
        | bgzip -@ {threads} \
        > {output.bed}

        # get bedgraph
        bedtools genomecov -bg -i {output.bed} -g {input.fai} | bgzip -@ {threads} > {output.bg}

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
        bed="results/{sm}/coverage/filtered-for-fdr/fiber-locations.bed.gz",
        shuffled="results/{sm}/coverage/filtered-for-fdr/fiber-locations-shuffled.bed.gz",
    threads: 1
    conda:
        conda
    shell:
        """
        cat {input.bed} > {output.bed}
        cat {input.shuffled} > {output.shuffled}
        """


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

#
# FIRE sites and tracks
#
rule fire_sites:
    input:
        bed=expand(rules.merge_model_results.output.bed, hp="all", allow_missing=True),
    output:
        bed="results/{sm}/FIRE.bed.gz",
    threads: 8
    conda:
        conda
    params:
        min_fdr=min_fire_fdr,
    shell:
        """
        bgzip -cd -@{threads} {input.bed} \
            | awk '$5<={params.min_fdr}' \
            | bgzip -@{threads} \
            > {output.bed}
        """


rule fire_tracks:
    input:
        fire=rules.fire_sites.output.bed,
        fiber_locations=rules.filtered_and_shuffled_fiber_locations.output.bed,
        shuffled=rules.filtered_and_shuffled_fiber_locations.output.shuffled,
        fai=ancient(f"{ref}.fai"),
    output:
        tbl="results/{sm}/FIRE.to.FDR.tbl",
    threads: 8
    conda:
        conda
    params:
        script=workflow.source_path("../scripts/fire-null-distribution.py"),
    shell:
        """
        python {params.script} -v 1 \
            {input.fire} \
            {input.fiber_locations} \
            {input.fai} \
            -s {input.shuffled} \
            {output.tbl}
        """
