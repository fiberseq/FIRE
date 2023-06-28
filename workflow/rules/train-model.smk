rule dhs_null:
    input:
        fai=ancient(f"{ref}.fai"),
        dhs=dhs,
        exclude=excludes,
    output:
        bed="results/{sm}/dhs_with_null.bed.gz",
        exclude=temp("temp/{sm}/exclues.bed"),
    threads: 2
    conda:
        conda
    shell:
        """
        less {input.exclude} | cut -f 1-3 | bedtools sort -i - | bedtools merge -i - > {output.exclude}
        (
            bedtools intersect -v -a {input.dhs} -b {output.exclude} \
                | cut -f 1-3 | sed 's/$/\tDHS/g'; 
            bedtools shuffle -noOverlapping \
                -excl <( less {input.dhs} {output.exclude} | cut -f 1-3 | bedtools sort -i - | bedtools merge -i -) \
                -i <(bedtools intersect -v -a {input.dhs} -b {output.exclude}) \
                -g {input.fai} |
                cut -f 1-3 | sed 's/$/\tNULL/g' 
        ) |
            sort -k 1,1 -k2,2n --parallel={threads} -S 5G |
            grep -vw "chrM" |
            grep -vw "chrY" |
            grep -vw "chrEBV" |
            grep -vw "chrX" |
            grep -v "_" |
            bgzip -@ {threads} > {output.bed}
        """


rule model_bam:
    input:
        bed=rules.dhs_null.output.bed,
        bam=lambda wc: data.loc[wc.sm, "bam"],
    output:
        bam=temp("temp/{sm}/dhs_and_null.bam"),
    threads: 16
    conda:
        conda
    shell:
        """
        samtools view -F 2308 -@ {threads} -M -L {input.bed} --write-index {input.bam} -o {output.bam}
        """


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
        d4tools create -t {threads} -Azr {input.fai} {input.bam} {output.d4}
        d4tools view -t {threads} {output.d4} | bgzip -@ {threads} > {output.bg}
        d4tools stat -t {threads} -s median {output.d4} > {output.median}
        """


rule filter_model_input_by_coverage:
    input:
        fai=ancient(f"{ref}.fai"),
        bed=rules.dhs_null.output.bed,
        bam=lambda wc: data.loc[wc.sm, "bam"],
        bg=rules.genome_bedgraph.output.bg,
        d4=rules.genome_bedgraph.output.d4,
    output:
        bed="results/{sm}/dhs_with_null_cov_filtered.bed",
    threads: 8
    conda:
        conda
    params:
        chrom=get_chroms()[0],
    shell:
        """
        median=$(d4tools stat -t {threads} -s median {input.d4} | grep -w {params.chrom} | cut -f 4)
        min=$(echo "$median" | awk '{{print $1-3*sqrt($1)}}')
        max=$(echo "$median" | awk '{{print $1+3*sqrt($1)}}')
        echo $median $min $max 
        bedmap --ec --delim '\t' --echo --mean \
            <(zcat {input.bed} | sort -k 1,1 -k2,2n) \
            <(zcat {input.bg} | awk '{{print $1"\t"$2"\t"$3"\t"$4"\t"$4}}' | sort -k 1,1 -k2,2n) \
            | awk -v min="$min" -v max="$max" '$5 > min && $5 < max' \
            | cut -f 1-4 \
            > {output.bed}
        head {output.bed}
        """


rule model_input:
    input:
        bam=lambda wc: data.loc[wc.sm, "bam"],
        dhs=rules.filter_model_input_by_coverage.output.bed,
    output:
        bed=temp("temp/{sm}/small.extract.all.bed.gz"),
    threads: 16
    params:
        n=500_000,
        sample_rate=0.05,
    conda:
        conda
    shell:
        """ 
        (samtools view -F 2308 -u -M -@ {threads} \
                -L {input.dhs} \
                {input.bam} -s {params.sample_rate} \
          | ft -t {threads} extract --all - \
          | bgzip -@ {threads} \
          > {output.bed} ) || echo "random head error"
        """


rule make_model:
    input:
        bed=rules.model_input.output.bed,
        dhs=rules.filter_model_input_by_coverage.output.bed,
    output:
        model=protected("results/{sm}/model.dat"),
    benchmark:
        "benchmarks/{sm}/make_model.tsv"
    threads: 60
    conda:
        conda
    params:
        n=100_000,
    shell:
        """
        fibertools -t {threads} -v model \
            -o {output.model} \
            --dhs {input.dhs} \
            {input.bed} 
        """
