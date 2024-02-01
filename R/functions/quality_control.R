## Subroutines for the main scripts

import::here(dplyr, 'count')
import::here(ggplot2, 'ggplot', 'aes', 'geom_point', 'labs')
import::here(file.path(wd, 'R', 'tools', 'plotting.R'),
    'savefig', 'plot_scatter', 'plot_violin', 'plot_waterfall', .character_only=TRUE)


## Functions
## compute_thresholds
## compute_cell_counts
## draw_qc_plots
## draw_predictions
## draw_clusters
## draw_gene_of_interest


#' Compute Thresholds
#'
compute_thresholds <- function(
    seurat_obj,
    sample_name = 'SeruatProject'
) {

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


#' Compute Cell Counts
#'
compute_cell_counts <- function(seurat_obj, gene, ident='cell_type') {

    seurat_obj[[gene]] <- seurat_obj[["RNA"]]@data[gene, ]

    num_total_cells <- count(seurat_obj@meta.data, .data[[ident]], name='num_cells')
    num_pos_cells <- count(
        filter(seurat_obj@meta.data, (.data[[gene]] > 5e-5)),
        cell_type, name='num_cells')
    cell_counts <- merge(
        num_total_cells, num_pos_cells,
        by = 'cell_type',
        all.x=TRUE,  # left join
        suffixes = c("_total", '_pos')
    )
    cell_counts <- fillna(cell_counts, cols=c('num_cells_pos'), 0)

    cell_counts['num_cells_neg'] <- cell_counts['num_cells_total'] - cell_counts['num_cells_pos']
    cell_counts['pct_cells_pos'] <- cell_counts['num_cells_pos'] / cell_counts['num_cells_total']

    return(cell_counts)
}


#' Draw QC plots
#' 
#' @description
#' 
draw_qc_plots <- function(
    seurat_obj,
    dirpath,
    prefix=NULL,
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
    savefig(file.path(dirpath, paste0('violin-counts_features_mt-', prefix, sample_name, '.png')),
            makedir=TRUE, troubleshooting=troubleshooting)

    
    # ----------------------------------------------------------------------
    # Figure 2. Waterfall Plot of Read Counts
    # See: https://sydneybiox.github.io/SingleCellPlus/qc.html#3_qc1:_waterfall_plot
    
    barcode_ranks <- barcodeRanks(seurat_obj)
    fig <- plot_waterfall(barcode_ranks)
    if (showfig) { print(fig) }
    savefig(file.path(dirpath, paste0('waterfall-count_vs_rank-', prefix, sample_name, '.png')),
            troubleshooting=troubleshooting)

    # ----------------------------------------------------------------------
    # Figure 3. Genes vs. Read Counts per cell

    fig <- plot_scatter(seurat_obj, group.by=group.by)
    if (showfig) { print(fig) }
    savefig(file.path(dirpath, paste0('scatter-features_vs_counts-', prefix, sample_name, '.png')),
            troubleshooting=troubleshooting)
   
    # ----------------------------------------------------------------------
    # Figure 4. Waterfall Plot of Gene Representation
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
    savefig(file.path(dirpath, paste0('histogram-cells_per_gene-', prefix, sample_name, '.png')),
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
    group_name='SeuratProject',
    troubleshooting=FALSE,
    showfig=FALSE
) {

    # ----------------------------------------------------------------------
    # Figure 1. Plot CellDex QC Heatmap
    
    fig <- plotScoreHeatmap(predictions)
    if (showfig) { print(fig) }
    savefig(file.path(dirpath, paste0('heatmap-prediction_score-', prefix, group_name, '.png')),
            fig=fig, lib='grid', height=1600, width=2400, dpi=300,
            makedir=TRUE, troubleshooting=troubleshooting)

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
    savefig(file.path(dirpath, paste0('histogram-cell_type-', prefix, group_name, '.png')),
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
    group_name='SeuratProject',
    split.by=NULL,
    troubleshooting=FALSE,
    showfig=FALSE
) {

    num_samples <- max(length(unname(unlist(unique( seurat_obj[[split.by]] )))), 1)

    # ----------------------------------------------------------------------
    # Figure 1. Seurat FindClusters Result

    fig <- DimPlot(seurat_obj,
            reduction = "umap",
            group.by = "seurat_clusters",
            split.by = split.by,
            label = TRUE) +
        theme(strip.text.x = element_blank(),
              plot.title = element_text(hjust = 0.5)) +
        ggtitle("Seurat FindClusters")
    if (showfig) { print(fig) }
    savefig(file.path(dirpath, paste0('umap-integrated-unlabeled-', prefix, group_name, '.png')),
            width=1000*num_samples,
            makedir=TRUE, troubleshooting=troubleshooting)


    # ----------------------------------------------------------------------
    # Figure 2. Plot CellDex labeled clusters

    fig <- DimPlot(seurat_obj,
            reduction = "umap", 
            group.by = "cell_type",
            split.by = split.by,
            label = TRUE
        ) + ggtitle(group_name)
    if (showfig) { print(fig) }
    savefig(file.path(dirpath, paste0('umap-integrated-labeled-', prefix, group_name, '.png')),
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
    group_name='SeuratProject',
    troubleshooting=FALSE,
    showfig=FALSE
) {

    seurat_obj[[gene]] <- seurat_obj[["RNA"]]@data[gene, ]

    # ----------------------------------------------------------------------
    # Figure 1. UMAP

    FeaturePlot(seurat_obj,
                reduction = "umap", features = gene,
                pt.size = 0.4, min.cutoff = 'q10', order = TRUE, label = FALSE) +
        ggtitle( paste(opt[['gene-of-interest']], 'in', group_name) )
    if (!troubleshooting) {
        savefig(file.path(wd, figures_dir, 'integrated', group_name, 'expression', tolower(gene),
                          paste0('umap-integrated-', group_name,  '-', prefix, tolower(gene), '.png')),
                height=800, width=800,
                troubleshooting=troubleshooting)
    }


    # ----------------------------------------------------------------------
    # Figure 2. Violin

    plot_violin(seurat_obj,
                cols=c(gene), group.by='cell_type',
                threshold_data=NULL, alpha=0.5)
    if (!troubleshooting) {
        savefig(file.path(wd, figures_dir, 'integrated', group_name, 'expression', tolower(gene),
                          paste0('violin-integrated-', group_name, '-', prefix, tolower(gene), '.png')),
                height=800, width=800,
                troubleshooting=troubleshooting)
    }


    # ----------------------------------------------------------------------
    # Figure 3. Ridge without 0

    RidgePlot(seurat_obj, gene) + xlim(5e-5, NA)
    if (!troubleshooting) {
        savefig(file.path(wd, figures_dir, 'integrated', group_name, 'expression', tolower(gene),
                          paste0('ridge-integrated-', group_name, '-', prefix, tolower(gene), '.png')),
                height=800, width=800,
                troubleshooting=troubleshooting)
    }

}
