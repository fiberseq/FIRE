library(tidyverse)
library(data.table)
library(scales)
library(ggforce)
library(cowplot)
library(dplyr)
library(splitstackshape)
library(ggridges)
library(IRanges)
library(ggrepel)
library(ggnewscale)
library(ggside)
library(glue)
library("tidylog", warn.conflicts = FALSE)
library(patchwork)
library(ggh4x)
library(tools)
library("purrr")
library("reticulate")
library(ggpubr)
library(ggbeeswarm)
library(weights)
library(karyoploteR)
library(zoo)
library("scales")
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
library(valr)
library(pals)


Red="#c1272d"
Indigo="#0000a7"
Yellow="#eecc16"
Teal="#008176"
Gray="#b3b3b3"

MODEL_COLORS = c(PacBio=Indigo, 
    CNN=Red,  
    XGB=Yellow,
    GMM=Teal,
    IPD=Gray,
    SEMI="purple",
    Revio="#f41c90" # black
)

read_m6a = function(file, my_tag = "", min_ml = 200, nrows=Inf, ref=TRUE){
    tmp = fread(glue(file), nrows=nrows)  %>%
        filter(en - st > 0.5 * fiber_length | (en == 0 & st == 0)) %>%
        filter(ec > 3.9)
    if ("sam_flag" %in% colnames(tmp)){
        print("filtering by sam flag")
        tmp = tmp %>% filter(sam_flag <= 16)
    }
        #[tmp$fiber %in% sample(unique(tmp$fiber), 500)] 
    tmp = tmp %>% 
        mutate(
            fake_end = dplyr::case_when(
                en == 0 ~ 11,
                TRUE ~ as.numeric(en)
            ),
            st=as.numeric(st),
            index = row_number(),
            bin = IRanges::disjointBins(
                IRanges(st+1, fake_end) + 150 #+ (fake_end - st+1) / 10
            ),
        )
    if(ref){
        m6a = tmp %>%
            select(ref_m6a, fiber, fiber_length, RG, m6a_qual, bin, st, en, strand) %>%
            filter(ref_m6a!=".") %>%
            cSplit(c("m6a_qual", "ref_m6a"), direction = "long") %>%
            filter(m6a_qual > min_ml) %>%
            mutate(
                type = "m6A",
                start = ref_m6a,
                end = ref_m6a + 1,
                alpha = 1.0,
                size = 2,
            )

        nuc = tmp %>%
            select(fiber, fiber_length, RG, ref_nuc_starts,  ref_nuc_lengths, bin, st, en, strand) %>%
            filter(ref_nuc_starts!=".") %>%
            cSplit(c("ref_nuc_starts", "ref_nuc_lengths"), direction = "long") %>%
            mutate(
                type = "Nucleosome",
                start = ref_nuc_starts,
                end = ref_nuc_starts + ref_nuc_lengths,
                alpha = 0.8,
                size = 1,
            )

        msp = tmp %>%
            select(fiber, fiber_length, RG, ref_msp_starts, ref_msp_lengths, bin, st, en, strand) %>%
            filter(ref_msp_starts!=".") %>%
            cSplit(c("ref_msp_starts", "ref_msp_lengths"), direction = "long") %>%
            mutate(
                type = "MSP",
                start = ref_msp_starts,
                end = ref_msp_starts + ref_msp_lengths,
                alpha = 0.8,
                size = 1,
            )
    } else {
        m6a = tmp %>%
            select(m6a, fiber, fiber_length, RG, m6a_qual, bin, st, en, strand) %>%
            filter(m6a!=".") %>%
            cSplit(c("m6a_qual", "m6a"), direction = "long") %>%
            filter(m6a_qual > min_ml) %>%
            mutate(
                type = "m6A",
                start = m6a,
                end = m6a + 1,
                alpha = 1.0,
                size = 2,
            )

        nuc = tmp %>%
            select(fiber, fiber_length, RG, nuc_starts, nuc_lengths, bin, st, en, strand) %>%
            filter(nuc_starts!=".") %>%
            cSplit(c("nuc_starts", "nuc_lengths"), direction = "long") %>%
            mutate(
                type = "Nucleosome",
                start = nuc_starts,
                end = nuc_starts + nuc_lengths,
                alpha = 0.8,
                size = 1,
            ) 
        
        msp = tmp %>%
            select(fiber, fiber_length, RG, msp_starts, msp_lengths, bin, st, en, strand) %>%
            filter(msp_starts!=".") %>%
            cSplit(c("msp_starts", "msp_lengths"), direction = "long") %>%
            mutate(
                type = "MSP",
                start = msp_starts,
                end = msp_starts + msp_lengths,
                alpha = 0.8,
                size = 1,
            )
    }

    print(my_tag)
    print(dim(tmp))
    print(dim(m6a))
    bind_rows(list(nuc, m6a, msp)) %>%
        mutate(tag = my_tag) %>%
        filter(start != -1) %>%
        group_by(type, fiber, RG) %>%
        arrange(start) %>%
        mutate(
            dist = start - lag(start)
        ) %>%
        data.table()
}


get_maxima_from_plot = function(p, n_keep=3 ){
    dens <- layer_data(p, 1)
    dens <- dens %>% arrange(colour, fill, x)

    # Run length encode the sign of difference
    rle <- rle(diff(as.vector(dens$y)) > 0)
    # Calculate startpoints of runs
    starts <- cumsum(rle$lengths) - rle$lengths + 1
    # Take the points where the rle is FALSE (so difference goes from positive to negative) 
    maxima_id <- starts[!rle$values]
    maxima <- dens[maxima_id,]
    maxima %>% group_by(colour, fill) %>% slice_max(order_by = y, n = n_keep)
}


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

library(ggplot2)
library(ggforce)
library(grid)

# define facet_zoom2 function to use FacetZoom2 instead of FacetZoom
# (everything else is the same as facet_zoom)
facet_zoom2 <- function(x, y, xy, zoom.data, xlim = NULL, ylim = NULL, 
                        split = FALSE, horizontal = TRUE, zoom.size = 2, 
                        show.area = TRUE, shrink = TRUE) {
  x <- if (missing(x)) if (missing(xy)) NULL else lazyeval::lazy(xy) else lazyeval::lazy(x)
  y <- if (missing(y)) if (missing(xy)) NULL else lazyeval::lazy(xy) else lazyeval::lazy(y)
  zoom.data <- if (missing(zoom.data)) NULL else lazyeval::lazy(zoom.data)
  if (is.null(x) && is.null(y) && is.null(xlim) && is.null(ylim)) {
    stop("Either x- or y-zoom must be given", call. = FALSE)
  }
  if (!is.null(xlim)) x <- NULL
  if (!is.null(ylim)) y <- NULL
  ggproto(NULL, FacetZoom2,
          shrink = shrink,
          params = list(
            x = x, y = y, xlim = xlim, ylim = ylim, split = split, zoom.data = zoom.data,
            zoom.size = zoom.size, show.area = show.area,
            horizontal = horizontal
          )
  )
}

# define FacetZoom as a ggproto object that inherits from FacetZoom,
# with a modified draw_panels function. the compute_layout function references
# the version currently on GH, which is slightly different from the CRAN
# package version.
FacetZoom2 <- ggproto(
  "FacetZoom2",
  ggforce::FacetZoom,

  compute_layout = function(data, params) {
    layout <- rbind( # has both x & y dimension
      data.frame(name = 'orig', SCALE_X = 1L, SCALE_Y = 1L),
      data.frame(name = 'x', SCALE_X = 2L, SCALE_Y = 1L),
      data.frame(name = 'y', SCALE_X = 1L, SCALE_Y = 2L),
      data.frame(name = 'full', SCALE_X = 2L, SCALE_Y = 2L),
      data.frame(name = 'orig_true', SCALE_X = 1L, SCALE_Y = 1L),
      data.frame(name = 'zoom_true', SCALE_X = 1L, SCALE_Y = 1L)
    )
    if (is.null(params$y) && is.null(params$ylim)) { # no y dimension
      layout <- layout[c(1,2, 5:6),]
    } else if (is.null(params$x) && is.null(params$xlim)) { # no x dimension
      layout <- layout[c(1,3, 5:6),]
    }
    layout$PANEL <- seq_len(nrow(layout))
    layout
  },

  draw_panels = function(panels, layout, x_scales, y_scales, ranges, coord,
                         data, theme, params) {

    if (is.null(params$x) && is.null(params$xlim)) {
      params$horizontal <- TRUE
    } else if (is.null(params$y) && is.null(params$ylim)) {
      params$horizontal <- FALSE
    }
    if (is.null(theme[['zoom']])) {
      theme$zoom <- theme$strip.background
    }
    if (is.null(theme$zoom.x)) {
      theme$zoom.x <- theme$zoom
    }
    if (is.null(theme$zoom.y)) {
      theme$zoom.y <- theme$zoom
    }
    axes <- render_axes(ranges, ranges, coord, theme, FALSE)
    panelGrobs <- ggforce:::create_panels(panels, axes$x, axes$y)
    panelGrobs <- panelGrobs[seq_len(length(panelGrobs) - 2)]
    if ('full' %in% layout$name && !params$split) {
      panelGrobs <- panelGrobs[c(1, 4)]
    }

    # changed coordinates in indicator / lines to zoom from 
    # the opposite horizontal direction
    if ('y' %in% layout$name) {
      if (!inherits(theme$zoom.y, 'element_blank')) {
        zoom_prop <- scales::rescale(
          y_scales[[2]]$dimension(ggforce:::expansion(y_scales[[2]])),
          from = y_scales[[1]]$dimension(ggforce:::expansion(y_scales[[1]])))
        indicator <- polygonGrob(
          x = c(0, 0, 1, 1), # was x = c(1, 1, 0, 0), 
          y = c(zoom_prop, 1, 0), 
          gp = gpar(col = NA, fill = alpha(theme$zoom.y$fill, 0.5)))
        lines <- segmentsGrob(
          x0 = c(1, 1), x1 = c(0, 0), # was x0 = c(0, 0), x1 = c(1, 1)
          y0 = c(0, 1), y1 = zoom_prop,
          gp = gpar(col = theme$zoom.y$colour,
                    lty = theme$zoom.y$linetype,
                    lwd = theme$zoom.y$size,
                    lineend = 'round'))
        indicator_h <- grobTree(indicator, lines)
      } else {
        indicator_h <- zeroGrob()
      }
    }

    if ('x' %in% layout$name) {
      if (!inherits(theme$zoom.x, 'element_blank')) {
        zoom_prop <- scales::rescale(x_scales[[2]]$dimension(ggforce:::expansion(x_scales[[2]])),
                                     from = x_scales[[1]]$dimension(ggforce:::expansion(x_scales[[1]])))
        indicator <- polygonGrob(c(zoom_prop, 1, 0), c(1, 1, 0, 0), 
                                 gp = gpar(col = NA, fill = alpha(theme$zoom.x$fill, 0.5)))
        lines <- segmentsGrob(x0 = c(0, 1), y0 = c(0, 0), x1 = zoom_prop, y1 = c(1, 1), 
                              gp = gpar(col = theme$zoom.x$colour,
                                        lty = theme$zoom.x$linetype,
                                        lwd = theme$zoom.x$size,
                                        lineend = 'round'))
        indicator_v <- grobTree(indicator, lines)
      } else {
        indicator_v <- zeroGrob()
      }
    }

    if ('full' %in% layout$name && params$split) {
      space.x <- theme$panel.spacing.x
      if (is.null(space.x)) space.x <- theme$panel.spacing
      space.x <- unit(5 * as.numeric(convertUnit(space.x, 'cm')), 'cm')
      space.y <- theme$panel.spacing.y
      if (is.null(space.y)) space.y <- theme$panel.spacing
      space.y <- unit(5 * as.numeric(convertUnit(space.y, 'cm')), 'cm')

      # change horizontal order of panels from [zoom, original] to [original, zoom]
      # final <- gtable::gtable_add_cols(panelGrobs[[3]], space.x)
      # final <- cbind(final, panelGrobs[[1]], size = 'first')
      # final_tmp <- gtable::gtable_add_cols(panelGrobs[[4]], space.x)
      # final_tmp <- cbind(final_tmp, panelGrobs[[2]], size = 'first')
      final <- gtable::gtable_add_cols(panelGrobs[[1]], space.x)
      final <- cbind(final, panelGrobs[[3]], size = 'first')
      final_tmp <- gtable::gtable_add_cols(panelGrobs[[2]], space.x)
      final_tmp <- cbind(final_tmp, panelGrobs[[4]], size = 'first')

      final <- gtable::gtable_add_rows(final, space.y)
      final <- rbind(final, final_tmp, size = 'first')
      final <- gtable::gtable_add_grob(final, list(indicator_h, indicator_h),
                                       c(2, 6), 3, c(2, 6), 5,
                                       z = -Inf, name = "zoom-indicator")
      final <- gtable::gtable_add_grob(final, list(indicator_v, indicator_v), 
                                       3, c(2, 6), 5, 
                                       z = -Inf, name = "zoom-indicator")
      heights <- unit.c(
        unit(max_height(list(axes$x[[1]]$top, axes$x[[3]]$top)), 'cm'),
        unit(1, 'null'),
        unit(max_height(list(axes$x[[1]]$bottom, axes$x[[3]]$bottom)), 'cm'),
        space.y,
        unit(max_height(list(axes$x[[2]]$top, axes$x[[4]]$top)), 'cm'),
        unit(params$zoom.size, 'null'),
        unit(max_height(list(axes$x[[2]]$bottom, axes$x[[4]]$bottom)), 'cm')
      )

      # swop panel width specifications according to the new horizontal order
      widths <- unit.c(
        # unit(max_width(list(axes$y[[3]]$left, axes$y[[4]]$left)), 'cm'),
        # unit(params$zoom.size, 'null'),
        # unit(max_height(list(axes$y[[3]]$right, axes$y[[4]]$right)), 'cm'),
        # space.x,
        # unit(max_width(list(axes$y[[1]]$left, axes$y[[2]]$left)), 'cm'),
        # unit(1, 'null'),
        # unit(max_height(list(axes$y[[1]]$right, axes$y[[2]]$right)), 'cm')        
        unit(max_width(list(axes$y[[1]]$left, axes$y[[2]]$left)), 'cm'),
        unit(1, 'null'),
        unit(max_height(list(axes$y[[1]]$right, axes$y[[2]]$right)), 'cm'),
        space.x,
        unit(max_width(list(axes$y[[3]]$left, axes$y[[4]]$left)), 'cm'),
        unit(params$zoom.size, 'null'),
        unit(max_height(list(axes$y[[3]]$right, axes$y[[4]]$right)), 'cm')

      )
      final$heights <- heights
      final$widths <- widths
    } else {
      if (params$horizontal) {
        space <- theme$panel.spacing.x
        if (is.null(space)) space <- theme$panel.spacing
        space <- unit(5 * as.numeric(convertUnit(space, 'cm')), 'cm')
        heights <- unit.c(
          unit(max_height(list(axes$x[[1]]$top, axes$x[[2]]$top)), 'cm'),
          unit(1, 'null'),
          unit(max_height(list(axes$x[[1]]$bottom, axes$x[[2]]$bottom)), 'cm')
        )

        # change horizontal order of panels from [zoom, original] to [original, zoom]
        # first <- gtable::gtable_add_cols(panelGrobs[[2]], space)
        # first <- cbind(final, panelGrobs[[1]], size = 'first')
        final <- gtable::gtable_add_cols(panelGrobs[[1]], space) 
        final <- cbind(final, panelGrobs[[2]], size = "first") 

        final$heights <- heights

        # swop panel width specifications according to the new horizontal order
        # unit(c(params$zoom.size, 1), 'null')
        final$widths[panel_cols(final)$l] <- unit(c(1, params$zoom.size), 'null') 

        final <- gtable::gtable_add_grob(final, indicator_h, 2, 3, 2, 5, 
                                         z = -Inf, name = "zoom-indicator")
      } else {
        space <- theme$panel.spacing.y
        if (is.null(space)) space <- theme$panel.spacing
        space <- unit(5 * as.numeric(convertUnit(space, 'cm')), 'cm')
        widths <- unit.c(
          unit(max_width(list(axes$y[[1]]$left, axes$y[[2]]$left)), 'cm'),
          unit(1, 'null'),
          unit(max_height(list(axes$y[[1]]$right, axes$y[[2]]$right)), 'cm')
        )
        final <- gtable::gtable_add_rows(panelGrobs[[1]], space)
        final <- rbind(final, panelGrobs[[2]], size = 'first')
        final$widths <- widths
        final$heights[panel_rows(final)$t] <- unit(c(1, params$zoom.size), 'null')
        final <- gtable::gtable_add_grob(final, indicator_v, 3, 2, 5, 
                                         z = -Inf, name = "zoom-indicator")
      }
    }
    final
  }
)


load_ipd_from_npz = function(file){
    np <- import("numpy")
    #npz <- np$load("tmp.npz")
    #npz <- np$load("ml_data/PS00075_2.val.npz")
    print(file)
    npz = np$load(file)
    fibers = npz["fibers"]
    labels = npz["labels"]
    n_rows = min(5e6, length(labels))
    print(n_rows)
    f = npz["features"]
    print("done reading features")
    ipd = f[,5,]
    raw_data = data.table(ipd)
    raw_data$label = labels
    raw_data$fiber = fibers
    raw_data %>% 
        pivot_longer(
            colnames(raw_data)[1:15]
        ) %>%
        mutate(
            name = factor(name, levels =paste0("V", seq(15)))
        )
}


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

my_read_bed = function(...){
    df=fread(...)
    names = colnames(df)
    names[names == "start"] = "start.other"
    names[names == "end"] = "end.other"
    names[1:3] = c("chrom", "start", "end")
    colnames(df) = names
    df
}

