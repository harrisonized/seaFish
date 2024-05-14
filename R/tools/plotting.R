wd = dirname(dirname(this.path::here()))
import::here(rlang, 'sym')
import::here(ggplot2,
    'ggplot', 'aes', 'aes_string', 'theme', 'theme_bw', 'labs', 
    'geom_violin', 'geom_boxplot', 'geom_jitter', 'geom_tile', 'geom_text',
    'geom_point', 'geom_hline', 'geom_vline', 'geom_segment', 'geom_bar',
    'guides', 'guide_axis', 'guide_legend', 'guide_colorbar',
    'xlim', 'ylim', 'coord_fixed', 'scale_fill_gradient', 'scale_radius',
    'scale_x_log10', 'scale_y_log10', 'scale_x_discrete', 'scale_y_discrete',
    'scale_color_gradient', 'scale_colour_gradientn',
    'annotate', 'element_text', 'element_blank', 'margin', 'unit')
import::here(scales, 'trans_breaks', 'trans_format', 'math_format')
import::here(cowplot, 'theme_cowplot')
import::here(patchwork, 'wrap_plots', 'plot_layout')
import::here(RColorBrewer, 'brewer.pal')
import::here(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'filter_list_for_match', .character_only=TRUE)
import::here(file.path(wd, 'R', 'tools', 'df_tools.R'),
    'rev_df', 'smelt', .character_only=TRUE)


## Functions
## plot_bar
## plot_dotplot
## plot_scatter
## plot_violin
## plot_volcano
## plot_waterfall
## plot_heatmap


#' Plot Bar
#' 
#' @description Plot a simple bar plot.
#' 
plot_bar <- function(
    df,
    x='cell_type',
    y='value',
    group.by=NULL,  # gene of interest
    fill=NULL,  # steelblue, overwrites group.by
    xlabel=NULL,
    ylabel="Number of Genes",
    title="Number of Cells",
    xaxis_angle=45,
    legend_title=NULL,
    legend_position='bottom',
    sort=TRUE
) {

    # color
    if (is.null(group.by)) { color <- NULL } else { color <- sym(group.by) }
    if (is.null(fill)) {
        bar <- geom_bar(stat="identity")
    } else {
        bar <- geom_bar(stat="identity", fill=fill)
    }

    # plot base
    if (sort) {
        base_plot <- ggplot(data=df,
            aes(x=reorder(.data[[x]], .data[[y]], decreasing=TRUE),
                y=.data[[y]], fill=!!color ))
    } else {
        base_plot <- ggplot(data=df,
            aes(x=.data[[x]], y=.data[[y]], fill=!!color ))
    }

    # legend
    if (!is.null(legend_title)) {
        guide <- guides( fill=guide_legend(title=legend_title) )
    } else {
        guide <- list()
    }

    fig <- base_plot +
        bar +
        labs(x=xlabel, y=ylabel, title=title) +
        scale_x_discrete( guide=guide_axis(angle = xaxis_angle) ) +
        theme(legend.position=legend_position) +
        guide

    return(fig)
}


#' Plot a Dot Plot
#' 
#' @description Plot a Dot Plot
#' 
#' @references
#' \href{https://satijalab.org/seurat/reference/dotplot}{Seurat}
#' 
plot_dotplot <- function(df,
    x='id',
    y='cell_type',
    size='pct_cells_pos',
    color='avg_expression',
    xlabel=NULL,
    ylabel=NULL,
    title="Expression",
    xaxis_angle=90,
    dot_scale=3,
    xtick_size=6,
    ytick_size=6,
    title_size=12,
    legend_size=7,
    legend_title='Percent Expressed',
    colorbar_title='Average Expression'
) {

    fig <- ggplot(data=df, 
            aes(x=.data[[x]],
                y =.data[[y]]) ) +
        geom_point(mapping=aes(size=.data[[size]], color=.data[[color]])) +
        scale_radius(range=c(0, dot_scale)) +
        labs(x=xlabel, y=ylabel, title=title) +
        scale_x_discrete(guide=guide_axis(angle=xaxis_angle)) +
        scale_y_discrete(limits=rev) +
        scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu"))) +
        # scale_color_gradient(low="lightgrey", high="blue") +
        theme_bw() +
        theme(axis.text.x = element_text(size=xtick_size, color="black"),
              axis.text.y = element_text(size=ytick_size, color="black"),
              plot.title = element_text(size=title_size),
              legend.title = element_text(size = legend_size),
              legend.text = element_text(size = legend_size),
              legend.key.size = unit(5, 'mm'),
              legend.margin=margin(c(20, 0, -5, 0)),
              panel.border = element_blank(), 
              panel.grid.major = element_blank(), ) +
        guides(size = guide_legend(title = legend_title),
               color = guide_colorbar(title = colorbar_title))

    return(fig)
}


#' Plot Scatter
#' 
#' @description Plot a simple 2D scatterplot
#' 
plot_scatter <- function(
    df,
    x='nCount_RNA',
    y='nFeature_RNA',
    color='sample_name',  # NULL for no groups
    xlabel="Read Counts",
    ylabel="Number of Genes",
    title="Read Counts vs Sequencing Depth",
    alpha=0.7,
    point_size=0.5,
    log_x=FALSE,
    log_y=FALSE,
    legend_large_circle=TRUE,
    jitter_height=0.05,
    jitter=FALSE
) {

    # jitter
    if (jitter) {
        geom_func <- geom_jitter(alpha=alpha, size=point_size, height=jitter_height, na.rm=TRUE)
    } else {
        geom_func <- geom_point(alpha=alpha, size=point_size, na.rm=TRUE)
    }

    # group.by
    if (is.null(color)) { colour <- NULL } else { colour <- sym(color) }

    # scale axes
    if (log_x) {
        scale_x <- scale_x_log10(
            breaks = trans_breaks("log10", function(x) 10^x),
            labels = trans_format("log10", math_format(10^.x))
        )
    } else { scale_x <- list() }
    if (log_y) {
        scale_y <- scale_y_log10(
            breaks = trans_breaks("log10", function(x) 10^x),
            labels = trans_format("log10", math_format(10^.x))
        )
    } else { scale_y <- list() }

    # legend
    if (legend_large_circle) {
        guide <- guides( colour=guide_legend(override.aes=list(size=4L,alpha=1)) )
    } else {
        guide <- list()
    }

    # plot
    fig <- ggplot(df,
              aes(x=.data[[x]], y=.data[[y]],
                  colour=!!colour) ) +
       geom_func +
       labs(x=xlabel, y=ylabel, title=title) + 
       guide +
       scale_x +
       scale_y

    return(fig)
}


#' Plot Violin
#' 
#' @description Similar to [Seurat::VlnPlot()], but prettier.
#' 
plot_violin <- function(
    seurat_obj,
    cols=c('nCount_RNA', 'nFeature_RNA', 'percent.mt'),
    group.by='orig.ident',
    threshold_data=data.frame(
        'sample_name'=c(rep('bm1', each=5), rep('bm2', each=5)),
        'col'=rep(c('nCount_RNA', 'nCount_RNA', 'nFeature_RNA', 'nFeature_RNA', 'percent.mt'), n=2),
        'threshold'=c(c(1000, 20000, 300, 2500, 5),
                      c(500, 10000, 500, 4000, 10))
    ),
    box=TRUE,
    pt.size=0.05,
    alpha=0.1,
    angle=45,
    xlabel=NULL,
    ylabel=NULL,
    title=NULL,
    title.size=12,
    showlegend=FALSE,
    legend.position='right'
) {
    subfigs <- new.env()
    for (i in 1:length(cols)) {
        col <- cols[i]

        # boxplot
        if (box) {
            boxplot <- geom_boxplot(width=0.1, alpha=1, outlier.size = pt.size)
        } else {
            boxplot <- list()
        }

        # thresholds
        if (is.null(threshold_data)) {
            segment_plot <- list()
        } else {

            threshold_data['group_num'] <- as.numeric(factor(threshold_data[[group.by]]))
            threshold_data[['x']] = threshold_data[['group_num']]-0.4
            threshold_data[['xend']] = threshold_data[['group_num']]+0.4
            threshold_data[['y']] = threshold_data[['threshold']]
            threshold_data[['yend']] = threshold_data[['threshold']]

            groups <- unname(unlist(unique( seurat_obj[[group.by]] )))

            segment_plot <- geom_segment(
                data = threshold_data[
                    (threshold_data[, group.by] %in% groups) &
                    (threshold_data['col']==col), 
                ],
                aes(x = x, y = y, xend = xend, yend = yend),
                inherit.aes=TRUE
            )
        }

        subfig <- ggplot(
            seurat_obj@meta.data,
            aes(x=.data[[group.by]],
                y=.data[[col]],
                fill=.data[[group.by]])) +
        geom_violin() +
        boxplot +
        geom_jitter(size = pt.size, alpha = alpha, na.rm=TRUE) +
        segment_plot +
        labs(x = ifelse(!is.null(xlabel), xlabel, col),
             y = ifelse(!is.null(ylabel), ylabel, col),
             title = ifelse(!is.null(title), title, col),
             fill=group.by) +
        cowplot::theme_cowplot() +
        theme(plot.title = element_text(size=title.size, hjust=0.5, vjust=0.5),
              axis.text.x = element_text(angle = angle, hjust=1),
              legend.title = element_text(size=12, face = "bold")) +
        ylim(0, NA)

        subfigs[[paste(i)]] <- subfig
    }

    fig <- patchwork::wrap_plots(as.list(subfigs)) +
        patchwork::plot_layout(guides = "collect") &
        theme(legend.position = ifelse(showlegend, legend.position, "none"))

    return(fig)
}


#' Plot Volcano
#' 
#' @description
#' 
plot_volcano <- function(
    df,
    x='avg_log2FC',
    y='p_val_adj',
    gene='gene',
    title="Volcano plot"
) {

    fig <- ggplot(df,
        aes(x=.data[['avg_log2FC']], y=-log10(.data[['p_val_adj']]),
            text=paste("Symbol:", .data[['gene']]))
    ) +
        geom_point(size=0.5, na.rm=TRUE) +
        labs(title=title) +
        theme_bw() +
        geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", linewidth=1) +
        geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=1) +
        geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=1) +
        annotate("rect", xmin = 1, xmax = 2, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#BE684D") +
        annotate("rect", xmin = -1, xmax = -2, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#2C467A")
    
    return(fig)
}


#' Plot Waterfall
#' 
#' @description Thin wrapper around plot_scatter
#' 
#' @references
#' \href{SingleCellPlus - HKU Workshop}{https://sydneybiox.github.io/SingleCellPlus/qc.html#3_qc1:_waterfall_plot}
#' 
plot_waterfall <- function(
    barcode_ranks,
    x='rank',
    y='total',
    xlabel="Barcode Rank",
    ylabel="UMI Counts",
    title="Barcode Rank Plot",
    show_thresholds=TRUE
) {

    if (show_thresholds) {
        thresholds <- data.frame(
            type = c("knee",
                     "inflection"),
            value = c(barcode_ranks@metadata$knee,
                      barcode_ranks@metadata$inflection)
        )
        hline <- geom_hline(
            data=thresholds,
            aes(yintercept = value, colour = type),
            linetype = 2
        )
    } else { hline <- list() }

    fig <- plot_scatter(
        as.data.frame(barcode_ranks),
        x=x, y=y, color=NULL,
        xlabel=xlabel, ylabel=ylabel, title=title,
        log_x=TRUE, log_y=TRUE,
        legend_large_circle = FALSE
    ) +
    hline +
    guides( color=guide_legend( override.aes=list(shape=c(NA, NA)) ))
    
    return(fig)
}


#' Plot a heatmap
#'
#' @usage
#' plot_heatmap(
#'   df,
#'   title="Raw Data",
#'   show_xlabel=TRUE,
#'   show_ylabel=TRUE,
#'   annotations=TRUE,
#'   scientific_notation=FALSE,
#'   digits=0
#' )
#' 
#' @section Vignette:
#' See `vignettes/plot_heatmap.Rmd`
#' 
#' @export
plot_heatmap <- function(
    df,
    x='col',
    y='row',
    fill='val',
    xlabel=NULL,
    ylabel=NULL,
    title=NULL,
    xlabel_order=NULL,
    annotations=FALSE,
    scientific_notation=FALSE,
    digits=1,
    xaxis_angle=0
) {
    
    tab <- smelt(rev_df(df))  # reshape
    if (!is.null(xlabel_order)) {
        tab[['col']] <- factor(tab[['col']], levels = xlabel_order)
    }
    
    # axis labels
    if (!is.null(xlabel)) {
        xtitle = element_text()
    } else {
        xtitle = element_blank()
    }
    if (!is.null(ylabel)) {
        ytitle = element_text()
    } else {
        ytitle = element_blank()
    }

    # annotations
    if (annotations) {
        label = 'label'
        if (scientific_notation) {
            tab['label'] = lapply(
                tab['val'], 
                function(x) formatC(x, format='e', digits=2)
            )
            tab['label']
        } else {
            tab['label'] = lapply(
                tab['val'], 
                function(x) round(x, digits)
            )
        }
        # tab[['label']] <- as.character(tab[['label']])
    } else {
        label = NULL
    }

    fig <- ggplot(tab, aes(x=.data[[x]], y=.data[[y]], fill=.data[[fill]])) +
        geom_tile(color="white", lwd=0.3, linetype=1, na.rm=FALSE) +
        scale_y_discrete(limits=rev) +
        coord_fixed(expand=TRUE) +
        labs(title = title, x=xlabel, y=ylabel) +
        scale_x_discrete(guide=guide_axis(angle=xaxis_angle)) +
        theme(plot.title = element_text(size = 10),
              axis.title.x = xtitle,
              axis.title.y = ytitle) +
        scale_fill_gradient(low="#FFF8F8", high="#A50026") +
        if (annotations) {
            geom_text(
                aes_string(label=label),
                color = 'black',
                size = 2,
                na.rm=TRUE
            )
        }

    return(fig)
}
