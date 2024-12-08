library(tidyverse)
library(data.table)
library(scales)
library(ggforce)
library(cowplot)
library(dplyr)
#library(splitstackshape)
#library(ggridges)
#library(IRanges)
library(ggrepel)
#library(ggnewscale)
#library(ggside)
library(glue)
#library("tidylog", warn.conflicts = FALSE)
#library(patchwork)
#library(ggh4x)
library(tools)
#library(purrr)
#library(reticulate)
library(ggpubr)
#library(weights)
#library(karyoploteR)
#library(zoo)
library(scales)
library(ggplot2)
library(ggforce)
library(grid)
options(scipen = 999)


Red="#c1272d"
Indigo="#0000a7"
Yellow="#eecc16"
Teal="#008176"
Gray="#b3b3b3"
FONT_SIZE=8
MODEL_COLORS = c(PacBio=Indigo, 
    CNN=Red,  
    XGB=Yellow,
    GMM=Teal,
    IPD=Gray,
    SEMI="purple",
    Revio="#f41c90" # black
)

my_ggsave <- function(file, ...){
    file = glue(file)
    print(file)
    ext = file_ext(file)
    file_without_ext = tools::file_path_sans_ext(file)
    ggsave(glue("tmp.{ext}"), bg='transparent', ...)
    cmd = glue("cp tmp.{ext} {file}")
    fwrite(last_plot()$data, file=file_without_ext + ".tbl.gz", sep="\t")
    print(cmd)
    system(cmd)
}

my_grid = function(...){
    theme_minimal_grid(font_size=FONT_SIZE, ...)
} 

my_hgrid = function(...){
    theme_minimal_hgrid(font_size=FONT_SIZE, ...)
} 

my_vgrid = function(...){
    theme_minimal_vgrid(font_size=FONT_SIZE, ...)
} 

reverselog_trans <- function(base = exp(1)) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv, 
              log_breaks(base = base), 
              domain = c(1e-100, Inf))
}

scientific_10 <- function(x) {
    is_one = as.numeric(x) == 1
    text = gsub("e", " %*% 10^", scales::scientific_format()(x))
    print(text)
    text = str_remove(text, "^1 %\\*% ") # remove leading one 
    print(text)
    text[is_one] = "10^0"
    rtn = parse(text=text)
    rtn
}

#
#
# SCRIPT
#
#
p_threshold=0.05
in_file=snakemake@input[[1]]
out_file_1=snakemake@output[[1]]
out_file_2=snakemake@output[[2]]
out_file_3=snakemake@output[[3]]
out_file_4=snakemake@output[[4]]


df=fread(in_file) %>%
    mutate_at(
        c("fire_coverage_H1","fire_coverage_H2","coverage_H1","coverage_H2"),
        as.numeric
    ) %>%
    filter(pass_coverage) %>%
    filter(
        !is.na(fire_coverage_H1),
        !is.na(fire_coverage_H2),
        !is.na(coverage_H1),
        !is.na(coverage_H2),
    ) %>%
    mutate(
        st=peak_start,
        en=peak_end,
        `#ct` = `#chrom`,
    ) %>%
    data.table()

print("Before coverage filter")
print(nrow(df))

# continue 
df$hap1_frac_acc = df$fire_coverage_H1/df$coverage_H1
df$hap2_frac_acc = df$fire_coverage_H2/df$coverage_H2
df$autosome = "Autosome"
df[`#chrom` == "chrY"]$autosome = "chrY"
df[`#chrom` == "chrX"]$autosome = "chrX"

# filter by coverage
sd = 3
cov = median(df$coverage)
my_min_cov = max(cov*0.5 - sd * sqrt(cov*0.5), 10)
my_max_cov = cov*0.5 + sd * sqrt(cov*0.5)
#my_min_cov = 10
#my_max_cov = median(df$coverage)

#print(glue("min_cov={my_min_cov} max_cov={my_max_cov} cov={cov}"))

pdf = df %>%
    filter(coverage_H1 > 0 & coverage_H2 > 0) %>%
    mutate(
        min_cov = my_min_cov,
        max_cov = my_max_cov,
    ) %>%
    filter(autosome != "chrY" ) %>%
    filter(coverage_H1 > min_cov & coverage_H2 > min_cov) %>%
    filter(coverage_H1 < max_cov & coverage_H2 < max_cov) %>%
    mutate(
        hap1_acc = fire_coverage_H1,
        hap2_acc = fire_coverage_H2,
        hap1_nacc = coverage_H1 - fire_coverage_H1,
        hap2_nacc = coverage_H2 - fire_coverage_H2,
    )

print("After coverage filter")
print(nrow(pdf))

if(nrow(pdf)== 0){
    system(glue("touch {out_file_1}"))
    system(glue("touch {out_file_2}"))
    system(glue("touch {out_file_3}"))
    system(glue("touch {out_file_4}"))
    quit()
}

pdf = pdf %>%
    filter(
      hap1_acc >= 0,
      hap1_nacc >= 0,
      hap2_acc >= 0,
      hap2_nacc >= 0,
    ) %>%
    rowwise() %>%
    mutate(
        p_value=fisher.test(matrix(c(hap1_acc, hap1_nacc, hap2_acc, hap2_nacc),nrow=2))$p.value
    ) %>%
    ungroup() %>%
    group_by(autosome) %>%
    mutate(
        p_adjust = p.adjust(p_value, method="BH"),
    ) %>%
    select(!starts_with("V")) %>%
    mutate(
        diff = hap1_frac_acc - hap2_frac_acc,
    ) %>%
    data.table()

print("After fisher test")
print(nrow(pdf))

print("After p-value adjustment")
print(nrow(pdf %>% filter(p_adjust <= p_threshold)))

# make the plots
tdf = pdf 

tdf %>%
    ggplot(aes(x=hap1_frac_acc, y=hap2_frac_acc)) +
    stat_cor(size=2) +
    geom_hex(bins=75) +
    geom_abline(aes(intercept=0, slope=1), linetype="dashed")+
    scale_fill_distiller("", palette = "Spectral", trans="log10") +
    scale_x_continuous("Paternal accessibility", labels=percent) +
    scale_y_continuous("Maternal accessibility", labels=percent) +
    #annotation_logticks(sides="lb") +
    facet_wrap(~autosome, ncol=2)+
    my_grid()
my_ggsave(out_file_1, height=3, width=6)

cor_p_threshold = max(tdf[p_adjust <= p_threshold & !is.na(p_value) & !is.na(p_adjust) ]$p_value)
print(cor_p_threshold)
y_lim = ceiling(max(max(-log10(tdf$p_value)), max(-log10(tdf$p_adjust)))) 
y_by = 1 
if(y_lim > 10){
    y_by = 2
}
# add p-value col, volcano plot
n=comma(nrow(tdf))
p = tdf %>%
    ggplot(aes(x=diff, y=p_value)) +
    geom_hex(bins=100) + scale_fill_distiller("", palette = "Spectral", trans="log10") +
    geom_hline(aes(yintercept=(p_threshold)), linetype="dashed", color="darkblue")+
    geom_hline(aes(yintercept=(cor_p_threshold)), linetype="dashed", color="darkred")+
    facet_wrap(~autosome, ncol=2)+
    scale_x_continuous("Difference between paternal and maternal accessibility", labels=percent) +
    scale_y_continuous(
        glue("p-value   (n = {n})"), 
        #labels=comma,
        breaks=10**(-seq(0, y_lim, y_by)),
        minor_breaks=10**(-seq(0, y_lim, 0.1)),
        trans=reverselog_trans(10),
        labels=scientific_10,
    ) + 
    my_grid()
my_ggsave(out_file_2, height=3, width=5)

# save the table 
fwrite(tdf[tdf$p_value <= 1], out_file_3, sep="\t")
# save the bed9
# 
tdf %>%
    mutate(
        score = pmin(round(1/p_value), 1000),
        name = paste0(
            round(p_value,5), "_",
            round(100*hap1_frac_acc,2), "_",
            round(100*hap2_frac_acc,2)
        ),
        tst=st+1,
        ten=case_when(
            p_value<=p_threshold ~ en - 1,
            TRUE ~ st+2
        ),
        strand=".",
        color = case_when(
            hap1_frac_acc > hap2_frac_acc & p_value > p_threshold ~ "0,0,100",
            hap1_frac_acc > hap2_frac_acc ~ "0,0,255",
            hap2_frac_acc > hap1_frac_acc & p_value > p_threshold ~ "100,0,0",
            hap2_frac_acc > hap1_frac_acc ~ "255,0,0",
            TRUE ~ "200,200,200"
        )
    ) %>%
    select(c("#chrom", "st", "en", "name", "score", "strand", "tst", "ten", "color")) %>%
    fwrite(out_file_4, sep="\t", quote=F, na=".", row.names=F)
