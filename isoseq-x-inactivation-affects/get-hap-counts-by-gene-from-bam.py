import pysam
import tqdm


tags = ["HP", "gn", "in", "tn"]
print("chrom\tstart\tend\tsample\tHP\tGene\tisoform\ttranscript")

for sample, f in [
    ("GM12878", "./gm12878.haplotagged.bam"),
    ("UDN318336", "./mat.pat.old.data.tagged.bam"),
]:
    bam = pysam.AlignmentFile(f, threads=16)

    for r in tqdm.tqdm(bam.fetch(until_eof=True)):
        if r.is_secondary or r.is_supplementary or r.is_unmapped:
            continue
        rtn = f"{r.reference_name}\t{r.reference_start}\t{r.reference_end}\t{sample}"
        for tag in tags:
            if r.has_tag(tag):
                res = str(r.get_tag(tag))
            else:
                res = "NA"
            rtn += f"\t{res}"
        print(rtn)
