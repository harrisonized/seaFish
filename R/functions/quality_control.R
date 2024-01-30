## Subroutines for the main scripts

import::here(ggplot2, 'ggplot', 'aes', 'geom_point', 'labs')
import::here(file.path(wd, 'R', 'tools', 'plotting.R'),
    'savefig', 'plot_scatter', 'plot_violin', 'plot_waterfall', .character_only=TRUE)


## Functions
## compute_thresholds
## draw_qc_plots
## draw_predictions
## draw_clusters


#' Compute Thresholds
#'
compute_thresholds <- function(seurat_obj) {

    thresholds <- new.env()  # upper thresholds
    for (col in c('nCount_RNA', 'nFeature_RNA', 'percent.mt')) {
        thresholds[[col]]  <- median( unlist(tmp_seurat_obj[[col]]) ) +
            3*sd( unlist(tmp_seurat_obj[[col]]) )
    }

    df <- data.frame(
        sample_name=sample_name,
        col=c('nCount_RNA', 'nCount_RNA',
              'nFeature_RNA', 'nFeature_RNA',
              'percent.mt'),
        type=c('min_nCount_RNA', 'max_nCount_RNA',
               'min_nFeature_RNA', 'max_nFeature_RNA',
               'max_percent.mt'),
        threshold=c(
            1000, thresholds[['nCount_RNA']],
            300, thresholds[['nFeature_RNA']],
            min(thresholds[['percent.mt']], 5)
        )
    )
    return(
        list(df,
            thresholds[['nCount_RNA']],
            thresholds[['nFeature_RNA']],
            thresholds[['percent.mt']])
    )
}


#' Draw QC plots
#' 
#' @description
#' 
draw_qc_plots <- function(
    seurat_obj,
    dirpath,
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
    savefig(file.path(dirpath, paste0('violin-qc-', sample_name, '.png')),
            makedir=TRUE, troubleshooting=troubleshooting)

    
    # ----------------------------------------------------------------------
    # Figure 2. Waterfall Plot of Read Counts
    # See: https://sydneybiox.github.io/SingleCellPlus/qc.html#3_qc1:_waterfall_plot
    
    barcode_ranks <- barcodeRanks(seurat_obj)
    fig <- plot_waterfall(barcode_ranks)
    if (showfig) { print(fig) }
    savefig(file.path(dirpath, paste0('waterfall-counts_per_cell-', sample_name, '.png')),
            troubleshooting=troubleshooting)

   
    # ----------------------------------------------------------------------
    # Figure 3. Waterfall Plot of Gene Representation
    # See: https://ucdavis-bioinformatics-training.github.io/2017_2018-single-cell-RNA-sequencing-Workshop-UCD_UCB_UCSF/day2/scRNA_Workshop-PART2.html
    
    cells_per_gene <- data.frame(
        num_cells=sort(
            rowSums(seurat_obj[['RNA']]@counts>=2),
            decreasing=TRUE)
    )
    cells_per_gene['rank'] = 1:nrow(seurat_obj[['RNA']]@counts)

    fig <- ggplot(cells_per_gene, aes(x = rank, y=num_cells)) +
        geom_point(size = 0.2) +
        labs(title = "Cells Per Gene ( >=2 )",
             x = "Gene Rank",
             y = "Number of Cells")
    if (showfig) { print(fig) }
    savefig(file.path(dirpath, paste0('histogram-cells_per_gene-', sample_name, '.png')),
            troubleshooting=troubleshooting)

    
    # ----------------------------------------------------------------------
    # Figure 4. Genes vs. Read Counts per cell

    fig <- plot_scatter(seurat_obj, group.by=group.by)
    if (showfig) { print(fig) }
    savefig(file.path(dirpath, paste0('scatter-reads_vs_depth-', sample_name, '.png')),
            troubleshooting=troubleshooting)

}



#' Draw Predictions QC
#' 
#' @description
#' Plot UMAPs after clustering
#' 
draw_predictions <- function(
    predictions,
    dirpath,
    group_name='SeuratProject',
    troubleshooting=FALSE,
    showfig=FALSE
) {

    # ----------------------------------------------------------------------
    # Figure 1. Plot CellDex QC Heatmap
    
    fig <- plotScoreHeatmap(predictions)
    if (showfig) { print(fig) }
    savefig(file.path(dirpath, paste0('heatmap-predictions-', group_name, '.png')),
            fig=fig, lib='grid', height=1600, width=2400, dpi=300, troubleshooting=troubleshooting)

    # ----------------------------------------------------------------------
    # Figure 2. Number of cells per label

    num_cells_per_label <- table(predictions['labels'])
    num_cells_per_label <- num_cells_per_label[sort.list(num_cells_per_label, decreasing=TRUE)]
    
    fig <- ggplot(data=data.frame(num_cells_per_label), aes(x=labels, y=Freq)) +
        geom_bar(stat="identity", fill="steelblue") +
        labs(title = "Number of Cells Per Label",
             x = NULL,
             y = "Number of Cells") +
        scale_x_discrete(guide = guide_axis(angle = 45))
    if (showfig) { print(fig) }
    savefig(file.path(dirpath, paste0('histogram-cells_type_dist-', group_name, '.png')),
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
    group_name='SeuratProject',
    split.by=NULL,
    troubleshooting=FALSE,
    showfig=FALSE
) {

    num_samples <- max(length(unname(unlist(unique( seurat_obj[[split.by]] )))), 1)

    # ----------------------------------------------------------------------
    # Seurat FindClusters Result

    fig <- DimPlot(seurat_obj,
            reduction = "umap",
            group.by = "seurat_clusters",
            split.by = split.by,
            label = TRUE) +
        theme(strip.text.x = element_blank(),
              plot.title = element_text(hjust = 0.5)) +
        ggtitle("Seurat FindClusters")
    if (showfig) { print(fig) }
    savefig(file.path(dirpath, paste0('umap-integrated-unlabeled-', group_name, '.png')),
            width=1000*num_samples, makedir=TRUE, troubleshooting=troubleshooting)


    # ----------------------------------------------------------------------
    # Plot CellDex labeled clusters

    fig <- DimPlot(seurat_obj,
            reduction = "umap", 
            group.by = "cell_type",
            split.by = split.by,
            label = TRUE
        ) + ggtitle(group_name)
    if (showfig) { print(fig) }
    savefig(file.path(dirpath, paste0('umap-integrated-labeled-', group_name, '.png')),
            width=1000*num_samples, troubleshooting=troubleshooting)
    
}
