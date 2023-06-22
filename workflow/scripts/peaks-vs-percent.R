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

df=fread(in_file)
df$peak_cov = df$peak_acc + df$peak_nuc + df$peak_link
df$acc_percent = df$peak_acc/df$peak_cov
df = df %>%
    group_by(sample) %>%
    arrange(-acc_percent) %>%
    mutate(
        min_cov = cov - 5*sqrt(cov),
        max_cov = cov + 5*sqrt(cov),
    ) %>%
    filter(
        peak_cov > min_cov,
        peak_cov < max_cov,
    ) %>%
    mutate(
        n=seq(n()),
        min_percent_acc=min(acc_percent)
    )

pecdf=df %>%
    #filter((n+1)%%100==0) %>%
    ggplot(aes(x=acc_percent, y=n)) +
    geom_line()+
    geom_text_repel(
        data = df %>% tail(1),
        aes(
            x=min_percent_acc,
            label=paste(
                "limit of detection:", percent(min_percent_acc, accuracy=0.01),
                "\n# peaks", comma(n)
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
        group = floor(acc_percent * 100) 
    ) %>%
    filter(group %% 5 == 0) %>%
    group_by(group) %>%
    slice_max(order_by = n, n = 1)

p5hist=by_5_per %>%
    ggplot(aes(x=acc_percent, y=n)) +
    geom_bar(stat="identity")+
    geom_text_repel(
        aes(
            x=acc_percent,
            label=paste(
                comma(n)
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