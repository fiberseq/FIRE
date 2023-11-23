source("Rscripts/utils.R")
library(regioneR)
library(BSgenome.Hsapiens.UCSC.hg38)
if(F){
    in_file="results/GM12878_FDR/FDR-peaks/FDR-FIRE-peaks.bed.gz"
    out_file="Figures/peaks-vs-percentage.pdf"
    encode=my_read_bed("data/ENCODE3_consensus_DHS_ENCFF503GCK.tsv")
    sds=my_read_bed("workflow/annotations/SDs.merged.hg38.bed.gz")
    tss=my_read_bed("data/gencode.v42.annotation_TSS.gff3")
    imprinted=my_read_bed("data/lcl_dmr_coordinates_Akbari.bed")

    dnase_peaks=my_read_bed("data/ENCFF762CRQ_DNase_peaks.bed.gz")
    dnase=my_read_bed("data/bedgraph_annotations/ENCFF960FMM_dnase_signal.bed")
    colnames(dnase)[4] = "dnase_sig"
    dnase
    
    atac_peaks = my_read_bed("scATAC/10X_GM12878_peaks_max_cov.bed")
    atac = my_read_bed("data/ATAC/10X_GM12878_aggr_scATAC.bg.gz")
    colnames(atac)[4] = "atac_sig"
    atac
    

    df=fread(in_file)
    df$peak_cov = df$coverage
    df$acc_percent = df$fire_coverage/df$peak_cov
    colnames(df)[1:5]=c("chrom", "start","end","ostart","oend")
    df$ID = paste0(df$chrom, "_", df$start, "_", df$end)

    df = df %>%
        arrange(-acc_percent) %>%
        mutate(
            n=seq(n()),
            min_percent_acc=min(acc_percent)
        ) %>%
        arrange(-n) %>%
        mutate(
            Mbp=cumsum(end-start)/1e6,
        ) %>%
        bed_map(encode,
            encode_count=length(core_end),
            encode_anno = paste0(unique(component), collapse=";")
        ) %>%
        bed_map(sds, sd_count=length(end)) %>%
        bed_map(dnase, dnase_max=max(dnase_sig)) %>%
        bed_map(tss, TSS=length(V4)) %>%
        bed_map(dnase_peaks, 
            is_dnase_peak = n()>0,
        ) %>%
        bed_map(imprinted, imprinted=(length(V4) > 0)) %>%
        bed_map(atac, atac_max = max(atac_sig)) %>%
        bed_map(atac_peaks, 
            is_atac_peak = n()>0,
        ) %>%
        replace_na(
            list(
                encode_count = 0,
                sd_count=0,
                TSS=0,
                is_atac_peak=F,
                imprinted=F
            )
        ) %>%
        arrange(-acc_percent, -n) %>%
        mutate(
            group = case_when(
                floor(acc_percent * 20) / 20 == 1 ~ 0.95,
                TRUE ~ floor(acc_percent * 20) / 20,
            )
        ) %>%
        data.table()

    #Build a GRanges from your matrix
    ranges <- toGRanges(df[,c("chrom", "start", "end")])

    #Get the sequences and compute the GC content
    freqs = alphabetFrequency(getSeq(BSgenome.Hsapiens.UCSC.hg38, ranges))
    df$GC_frac = (freqs[,'C'] + freqs[,'G'])/rowSums(freqs)

    # GC correction
    df$psize=df$end-df$start
    fit = lm(acc_percent ~ `GC_frac` * log10(dnase_max), data=df[is_dnase_peak==T,])
    summary(fit)
    df$dnase_max_gc_corrected = predict(fit, newdata=df)

    fit = lm(acc_percent ~ log10(atac_max), data=df[is_atac_peak==T,])
    summary(fit)
    fit = lm(acc_percent ~ `GC_frac` * log10(atac_max), data=df[is_atac_peak==T,])
    summary(fit)
    df$atac_max_gc_corrected = predict(fit, newdata=df) 

    fire_df = df



    system("mkdir -p Rdata")
    con <- pipe("pigz -p8 > Rdata/df.fire-peaks.gz", "wb")
    save(
        fire_df,
        encode,
        dnase_peaks,
        atac_peaks,
        sds,
        tss,
        file = con
    ); close(con)
} else{
    load("Rdata/df.fire-peaks.gz")
}

