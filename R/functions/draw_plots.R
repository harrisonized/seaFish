## Functions that draws groups of figures

import::here(magrittr, '%>%')
import::here(dplyr, 'select')
import::here(tidyselect, 'all_of')
import::here(tidyr, 'pivot_longer')
import::here(Seurat, 'DimPlot', 'FeaturePlot', 'RidgePlot')
import::here(DropletUtils, 'barcodeRanks')
import::here(SingleR, 'plotScoreHeatmap')
import::here(ggplot2, 'ggplot', 'aes',  'theme', 'labs',
    'xlim', 'ggtitle', 'element_text', 'element_blank')

import::here(file.path(wd, 'R', 'tools', 'file_io.R'),
    'savefig', .character_only=TRUE)
import::here(file.path(wd, 'R', 'tools', 'plotting.R'),
    'plot_violin', 'plot_scatter', 'plot_waterfall', 'plot_bar', 
    .character_only=TRUE)
import::here(file.path(wd, 'R', 'functions', 'computations.R'),
    'compute_cell_counts', .character_only=TRUE)

## Functions
## draw_qc_plots
## draw_predictions
## draw_clusters
## draw_gene_of_interest


#' Draw QC plots
#' 
#' @description
#' 
draw_qc_plots <- function(
    seurat_obj,
    dirpath,
    prefix=NULL,
    suffix=NULL,
    sample_name='SeuratProject',
    group.by='sample_name',
    threshold_data = NULL,
    troubleshooting=FALSE,
    showfig=FALSE
) {

    # ----------------------------------------------------------------------
    # Figure 1. Standard QC Metrics
    
    fig <- plot_violin(seurat_obj, group.by=group.by, threshold_data=threshold_data)
    if (showfig) { print(fig) }
    savefig(file.path(dirpath, paste0('violin-counts_features_mt-', prefix, sample_name, suffix, '.png')),
            makedir=TRUE, troubleshooting=troubleshooting)

    # ----------------------------------------------------------------------
    # Figure 2. Waterfall Plot of Read Counts
    # See: https://sydneybiox.github.io/SingleCellPlus/qc.html#3_qc1:_waterfall_plot
    
    barcode_ranks <- barcodeRanks(seurat_obj)
    fig <- plot_waterfall(barcode_ranks)
    if (showfig) { print(fig) }
    savefig(file.path(dirpath, paste0('waterfall-count_vs_rank-', prefix, sample_name, suffix, '.png')),
            troubleshooting=troubleshooting)

    # ----------------------------------------------------------------------
    # Figure 3. Genes vs. Read Counts per cell

    fig <- plot_scatter(seurat_obj@meta.data, color=group.by)
    if (showfig) { print(fig) }
    savefig(file.path(dirpath, paste0('scatter-features_vs_counts-', prefix, sample_name, suffix, '.png')),
            troubleshooting=troubleshooting)
   
    # ----------------------------------------------------------------------
    # Figure 4. Waterfall Plot of Gene Representation
    # See: https://ucdavis-bioinformatics-training.github.io/2017_2018-single-cell-RNA-sequencing-Workshop-UCD_UCB_UCSF/day2/scRNA_Workshop-PART2.html
    
    cells_per_gene <- data.frame(
        num_cells=sort(rowSums(seurat_obj[['RNA']]@counts>=2), decreasing=TRUE)
    )
    cells_per_gene['rank'] = 1:nrow(seurat_obj[['RNA']]@counts)

    fig <- plot_scatter(
        cells_per_gene,
        x='rank',  y='num_cells', color=NULL,
        xlabel="Number of Cells", ylabel="Gene Rank",
        title="Cells Per Gene ( >=2 )",
        point_size=0.2)
    if (showfig) { print(fig) }
    savefig(file.path(dirpath, paste0('histogram-cells_per_gene-', prefix, sample_name, suffix, '.png')),
            troubleshooting=troubleshooting)
}


#' Draw Predictions QC
#' 
#' @description
#' 
draw_predictions <- function(
    predictions,
    dirpath,
    prefix=NULL,
    suffix=NULL,
    group_name='SeuratProject',
    troubleshooting=FALSE,
    showfig=FALSE
) {

    # ----------------------------------------------------------------------
    # Figure 1. Plot CellDex QC Heatmap
    
    fig <- plotScoreHeatmap(predictions)
    if (showfig) { print(fig) }
    savefig(file.path(dirpath, paste0('heatmap-prediction_score-', prefix, group_name, suffix, '.png')),
            fig=fig, lib='grid', height=1600, width=2400, dpi=300,
            makedir=TRUE, troubleshooting=troubleshooting)

    # ----------------------------------------------------------------------
    # Figure 2. Number of cells per label

    num_cells_per_label <- table(predictions['labels'])
    num_cells_per_label <- num_cells_per_label[sort.list(num_cells_per_label, decreasing=TRUE)]

    plot_bar(
        data.frame(num_cells_per_label),
        x='labels',
        y='Freq',
        fill='steelblue',
        group.by=NULL,  # gene of interest
        xlabel=NULL,
        ylabel="Number of Cells",
        title="Number of Cells Per Label"
    )
    if (showfig) { print(fig) }
    savefig(file.path(dirpath, paste0('histogram-cell_type-', prefix, group_name, suffix, '.png')),
            height=800, width=1200, dpi=400, troubleshooting=troubleshooting)
}


#' Draw Clusters
#' 
#' @description
#' Plot UMAPs after clustering
#' 
draw_clusters <- function(
    seurat_obj,
    dirpath,
    prefix=NULL,
    suffix=NULL,
    celldex_dataset='ImmGen',
    group_name='SeuratProject',
    split.by=NULL,
    troubleshooting=FALSE,
    showfig=FALSE
) {

    num_samples <- max(length(unname(unlist(unique( seurat_obj[[split.by]] )))), 1)

    # ----------------------------------------------------------------------
    # Figure 1. Seurat Clusters

    fig <- DimPlot(seurat_obj,
            reduction = "umap",
            group.by = "seurat_clusters",
            split.by = split.by,
            label = TRUE) +
        theme(strip.text.x = element_blank(),
              plot.title = element_text(hjust = 0.5)) +
        ggtitle("Seurat Annotations")
    if (showfig) { print(fig) }
    savefig(file.path(dirpath, paste0('umap-seurat_labeled-', prefix, group_name, suffix, '.png')),
            width=1000*num_samples,
            makedir=TRUE, troubleshooting=troubleshooting)

    # ----------------------------------------------------------------------
    # Figure 2. CellDex Clusters

    fig <- DimPlot(seurat_obj,
            reduction = "umap", 
            group.by = "cell_type",
            split.by = split.by,
            label = TRUE
        ) + ggtitle(group_name)
    if (showfig) { print(fig) }
    savefig(file.path(dirpath, paste0('umap-', tolower(celldex_dataset), '_labeled-', prefix, group_name, suffix, '.png')),
            width=1000*num_samples, troubleshooting=troubleshooting)
}


#' Draw Gene of Interest
#' 
#' @description
#' 
draw_gene_of_interest <- function(
    seurat_obj,
    gene,
    dirpath,
    prefix=NULL,
    suffix=NULL,
    group_name='SeuratProject',
    troubleshooting=FALSE,
    showfig=FALSE
) {

    seurat_obj[[gene]] <- seurat_obj[["RNA"]]@data[gene, ]

    # ----------------------------------------------------------------------
    # Figure 1. UMAP

    fig <- FeaturePlot(seurat_obj,
                reduction = "umap", features = gene,
                pt.size = 0.4, min.cutoff = 'q10', order = TRUE, label = FALSE) +
        ggtitle( paste(opt[['gene-of-interest']], 'in', group_name) )
    if (showfig) { print(fig) }
    savefig(file.path(dirpath, paste0('umap-', prefix, group_name, suffix, '-', tolower(gene), '.png')),
            height=800, width=800,
            troubleshooting=troubleshooting)

    # ----------------------------------------------------------------------
    # Figure 2. Violin

    fig <- plot_violin(seurat_obj,
                cols=c(gene), group.by='cell_type',
                threshold_data=NULL, alpha=0.5)
    if (showfig) { print(fig) }
    savefig(file.path(dirpath, paste0('violin-', prefix, group_name, suffix, '-', tolower(gene), '.png')),
            height=800, width=800,
            troubleshooting=troubleshooting)

    # ----------------------------------------------------------------------
    # Figure 3. Ridge without 0

    fig <- RidgePlot(seurat_obj, gene) + xlim(5e-5, NA)
    if (showfig) { print(fig) }
    savefig(file.path(dirpath, paste0('ridge-', prefix, group_name, suffix, '-', tolower(gene), '.png')),
            height=800, width=800,
            troubleshooting=troubleshooting)

    # ----------------------------------------------------------------------
    # Figure 4. Bar plot

    cell_counts <- compute_cell_counts(seurat_obj, gene=gene, ident='cell_type')

    value_cols = c('num_cells_neg', 'num_cells_pos')
    cell_counts_long <- cell_counts %>%
        dplyr::select(tidyselect::all_of(c('cell_type', value_cols))) %>%
        tidyr::pivot_longer(cols=value_cols)
    cell_counts_long[gene] <- sapply(cell_counts_long['name'], function(x) gsub('num_cells_', '', x))
    
    fig <- plot_bar(
        cell_counts_long[order(-cell_counts_long$value, decreasing = FALSE), ],
        x='cell_type',
        y='value',
        group.by=gene,  # gene of interest
        xlabel=NULL,
        ylabel="Number of Genes",
        title=paste0('Number of ', gene, '+ Cells')
    )
    if (showfig) { print(fig) }   
    savefig(file.path(dirpath, paste0('histogram-cell_type-', prefix, group_name, suffix, '-', tolower(gene), '.png')),
            height=800, width=1200, dpi=400, troubleshooting=troubleshooting)
}
