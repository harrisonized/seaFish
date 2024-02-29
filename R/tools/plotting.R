wd = dirname(dirname(this.path::here()))
import::here(rlang, 'sym')
import::here(ggplot2,
    'ggplot', 'aes', 'theme', 'theme_bw', 'labs', 
    'geom_violin', 'geom_boxplot', 'geom_jitter',
    'geom_point', 'geom_hline', 'geom_vline', 'geom_segment', 'geom_bar',
    'guides', 'guide_axis', 'guide_legend',
    'xlim', 'ylim', 'scale_x_log10', 'scale_y_log10', 'scale_x_discrete',
    'annotate', 'element_text', 'element_blank')
import::here(scales, 'trans_breaks', 'trans_format', 'math_format')
import::here(cowplot, theme_cowplot)
import::here(patchwork, 'wrap_plots', 'plot_layout')
import::here(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'filter_list_for_match', .character_only=TRUE)

## Functions
## plot_bar
## plot_scatter
## plot_violin
## plot_volcano
## plot_waterfall


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
