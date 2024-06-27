library(tidyverse)
library(data.table)
library(scales)
library(ggforce)
library(cowplot)
library(dplyr)
library(ggrepel)
library(glue)
library(patchwork)
library(ggpubr)

FONT_SIZE=6
my_grid = function(...){
    theme_minimal_grid(font_size=FONT_SIZE, ...)
} 

my_hgrid = function(...){
    theme_minimal_hgrid(font_size=FONT_SIZE, ...)
} 

my_vgrid = function(...){
    theme_minimal_vgrid(font_size=FONT_SIZE, ...)
} 

theme_no_x = function(...){
    theme(
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()
    )
}
in_file=snakemake@input[[1]]
out_file=snakemake@output[[1]]

# chrom  peak_start      peak_end        start   end     
# fire_coverage   coverage        score   FDR     log_FDR 
# fire_coverage_H1        coverage_H1     score_H1        FDR_H1  log_FDR_H1
# fire_coverage_H2        coverage_H2  score_H2 FDR_H2  log_FDR_H2      
# FIRE_size_mean  FIRE_size_ssd   FIRE_start_ssd  FIRE_end_ssd    local_max_count peak_length


df=fread(in_file)
df$peak_cov = df$coverage
df$acc_percent = df$fire_coverage/df$peak_cov
df = df %>%
    filter(pass_coverage) %>%
    arrange(-acc_percent) %>%
    mutate(
        count=row_number(),
        min_percent_acc=min(acc_percent)
    )

pecdf=df %>%
    arrange(-acc_percent, -count) %>%
    #filter((n+1)%%100==0) %>%
    ggplot(aes(x=acc_percent, y=count)) +
    geom_line()+
    geom_text_repel(
        data = df %>% tail(1),
        aes(
            x=min_percent_acc,
            label=paste(
                "limit of detection:", percent(min_percent_acc, accuracy=0.01),
                "\n# peaks", comma(count)
            ),
        ),
        min.segment.length=0,
        size=2,
        nudge_y=-1.5,
        nudge_x=0.1,
    ) +
    scale_y_continuous(
        "# of regulatory elements in the genome",
        trans="log10",
        label=comma
    ) + 
    annotation_logticks(side="l")+
    scale_x_continuous(
        "Minimum % of fibers that are accessible",
        breaks=seq(5,100,5)/100,
        label=percent,
        #guide = guide_axis(n.dodge = 2),
    ) +
    my_grid() +
    coord_cartesian(xlim=c(0,1)) 

by_5_per = df %>%
    mutate(
        group = floor(20 * acc_percent ) / 20 * 100 
    ) %>%
    #filter(group %% 5 == 0) %>%
    group_by(group) %>%
    filter(count == max(count)) %>%
    slice_max(order_by = count, n = 1) %>%
    select(group, acc_percent, count) 

print(by_5_per, nrow=25)


p5hist=by_5_per %>%
    ggplot(aes(x=(group+2.5)/100, y=count)) +
    geom_bar(stat="identity")+
    geom_text_repel(
        aes(
            #x=acc_percent-2.5/100,
            label=paste(
                comma(count)
            ),
        ),
        min.segment.length=0,
        direction="y",
        nudge_y=10,
        size=1
    ) +
    scale_x_continuous(
        "",
        #"Minimum % of fibers that are accessible", 
        breaks=seq(5,100,5)/100,
        label=rep("",20),
    )+
    scale_y_continuous("", label=comma) +
    my_grid() +
    coord_cartesian(xlim=c(0,1))+
    theme_no_x()


z=plot_grid(p5hist,pecdf, ncol=1, align="v", rel_heights=c(1,4))
ggsave(out_file, height=3, width=5)