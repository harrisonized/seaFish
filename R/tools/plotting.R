# library('Matrix')
wd = dirname(dirname(this.path::here()))
import::here(ggplot2, 'ggsave')
import::here(grid, 'grid.newpage', 'grid.draw')
import::here(patchwork, 'wrap_plots', 'plot_layout')
import::here(cowplot, theme_cowplot)
import::here(scales, 'trans_breaks', 'trans_format', 'math_format')
import::here(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'filter_list_for_match', .character_only=TRUE)

## Functions
## save_fig
## plot_violin
## plot_waterfall


#' Save Figure
#' 
#' @description Switch case to reduce the number of lines in the main script
#' 
savefig <- function(
    filepath,
    fig=NULL,
    height=800, width=1200, dpi=300, units="px", scaling=0.5,
    makedir=FALSE,
    troubleshooting=FALSE,
    lib='ggplot'  # choose: ggplot, grid
) {
    if (!troubleshooting) {

        # make directory
        dirpath <- dirname(filepath)
        if (makedir && !dir.exists(dirpath)) {
            dir.create(dirpath, recursive=TRUE)
        }

        if (lib=='ggplot') {
            ggsave(
                filepath,
                height=height, width=width, dpi=dpi, units=units, scaling=scaling
            )
        } else if (lib=='grid') {
            png(filepath,
                height=height, width=width, res=dpi, units=units
            )
            grid.newpage()
            grid.draw(fig$gtable)
            dev.off()
        } else {
            warning(paste0("lib='", lib, "' not found"))
        }
    }
}


#' Plot Violin
#' 
#' @description
#' Similar to [Seurat::VlnPlot()], but prettier.
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
    pt.size = 0.05,
    alpha = 0.1,
    angle = 45,
    title.size = 12,
    showlegend=FALSE,
    legend.position = 'right'
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
        geom_jitter(size = pt.size, alpha = alpha) +
        segment_plot +
        labs(title = col,
             x = NULL,
             y = NULL,
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


#' Plot Waterfall
#' 
#' @references
#' \href{https://sydneybiox.github.io/SingleCellPlus/qc.html#3_qc1:_waterfall_plot}
#'  
plot_waterfall <- function(barcode_ranks) {

    barcode_points = data.frame(
        type = c("knee", "inflection"),
        value = c(barcode_ranks@metadata$knee, barcode_ranks@metadata$inflection)
    )

    fig <- ggplot(data = as.data.frame(barcode_ranks), aes(x = rank, y = total)) +
        geom_point() +
        geom_hline(data = barcode_points,
                   aes(yintercept = value,
                       colour = type), linetype = 2) +
        scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x))) +
        scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x))) +
        labs(title = "Counts Per Cell",
             x = "Cell Rank",
             y = "Read Counts")

    return(fig)
}
