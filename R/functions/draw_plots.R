## Functions that draw groups of figures

import::here(magrittr, '%>%')
import::here(dplyr, 'select')
import::here(tidyselect, 'all_of')
import::here(tidyr, 'pivot_longer', 'pivot_wider')
import::here(EnhancedVolcano, 'EnhancedVolcano')
import::here(Seurat, 'DimPlot', 'FeaturePlot', 'RidgePlot')
import::here(DropletUtils, 'barcodeRanks')
import::here(SingleR, 'plotScoreHeatmap')
import::here(ggplot2, 'ggplot', 'aes',  'theme', 'labs',
    'xlim', 'ggtitle', 'element_text', 'element_blank')

import::here(file.path(wd, 'R', 'tools', 'df_tools.R'),
    'set_index', .character_only=TRUE)
import::here(file.path(wd, 'R', 'tools', 'file_io.R'),
    'savefig', .character_only=TRUE)
import::here(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'items_in_a_not_b', .character_only=TRUE)
import::here(file.path(wd, 'R', 'tools', 'text_tools.R'),
    'snake_to_title_case', 'title_to_snake_case', .character_only=TRUE)
import::here(file.path(wd, 'R', 'functions', 'computations.R'),
    'compute_gene_expression', 'compute_gene_labels', 'compute_p_cutoff',
    .character_only=TRUE)
import::here(file.path(wd, 'R', 'tools', 'plotting.R'),
    'plot_violin', 'plot_scatter', 'plot_waterfall', 'plot_bar', 'plot_heatmap',
    .character_only=TRUE)
import::here(file.path(wd, 'R', 'tools', 'single_cell_tools.R'),
    'FindDESeq2Markers', .character_only=TRUE)

## Functions
## draw_qc_plots
## draw_clusters
## draw_gene_of_interest
## draw_predictions
## draw_differential_genes


#' Draw QC plots
#' 
#' @description
#' 
draw_qc_plots <- function(
    seurat_obj,
    group.by='sample_name',
    threshold_data = NULL,  # dataframe
    dirpath,
    file_basename='SeuratProject',
    troubleshooting=FALSE,
    showfig=FALSE
) {

    # ----------------------------------------------------------------------
    # Figure 1. Standard QC Metrics

    fig <- plot_violin(seurat_obj, group.by=group.by, threshold_data=threshold_data)
    if (showfig) { print(fig) }
    savefig(file.path(dirpath, paste0('violin-counts_features_mt-', file_basename, '.png')),
            makedir=TRUE, troubleshooting=troubleshooting)

    # ----------------------------------------------------------------------
    # Figure 2. Waterfall Plot of Read Counts
    # See: https://sydneybiox.github.io/SingleCellPlus/qc.html#3_qc1:_waterfall_plot
    
    barcode_ranks <- barcodeRanks(seurat_obj)
    fig <- plot_waterfall(barcode_ranks)
    if (showfig) { print(fig) }
    savefig(file.path(dirpath, paste0('waterfall-count_vs_rank-', file_basename, '.png')),
            troubleshooting=troubleshooting)

    # ----------------------------------------------------------------------
    # Figure 3. Genes vs. Read Counts per cell

    fig <- plot_scatter(seurat_obj@meta.data, color=group.by)
    if (showfig) { print(fig) }
    savefig(file.path(dirpath, paste0('scatter-features_vs_counts-', file_basename, '.png')),
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
    savefig(file.path(dirpath, paste0('histogram-cells_per_gene-', file_basename, '.png')),
            troubleshooting=troubleshooting)
}


#' Draw Clusters
#' 
#' @description
#' Plot UMAPs after clustering
#' 
draw_clusters <- function(
    seurat_obj,
    celldex_dataset='ImmGen',
    split.by=NULL,
    dirpath,
    file_basename='SeuratProject',
    title='SeuratProject',
    include_heatmap=FALSE,
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
        ggtitle(title)
    if (showfig) { print(fig) }
    savefig(file.path(dirpath, paste0('umap-seurat_labeled-', file_basename, '.png')),
            width=4000*num_samples+800, height=4000, dpi=1600,
            makedir=TRUE, troubleshooting=troubleshooting)

    # ----------------------------------------------------------------------
    # Figure 2. CellDex Clusters

    fig <- DimPlot(seurat_obj,
            reduction = "umap", 
            group.by = "cell_type",
            split.by = split.by,
            label = TRUE
        ) + ggtitle(title)
    if (showfig) { print(fig) }
    savefig(file.path(dirpath, paste0('umap-', tolower(celldex_dataset), '_labeled-', file_basename, '.png')),
            width=4000*num_samples+800, height=4000, dpi=1600,
            troubleshooting=troubleshooting)

    # ----------------------------------------------------------------------
    # Figure 3. Seurat vs. Celldex

    if (include_heatmap) {
        seurat_vs_celldex <- as.data.frame(table(
                seurat_obj@meta.data[, c('seurat_clusters', 'cell_type')])) %>%
            pivot_wider(
                id_cols = 'cell_type',
                names_from='seurat_clusters',
                values_from='Freq') %>%
            as.data.frame()
        
        fig <- plot_heatmap(
                set_index(seurat_vs_celldex, 'cell_type'),
                xlabel='Seurat Cluster',
                ylabel=paste(celldex_dataset, 'Label'),
                annotations=TRUE
            )
        if (showfig) { print(fig) }
        savefig(file.path(dirpath, paste0('heatmap-', 'seurat_vs_', tolower(celldex_dataset), '-', file_basename, '.png')),
                height=800, width=1200, troubleshooting=troubleshooting)
    }
}


#' Draw Gene of Interest
#' 
#' @description
#' 
draw_gene_of_interest <- function(
    seurat_obj,
    gene,
    dirpath,
    file_basename='SeuratProject',
    troubleshooting=FALSE,
    showfig=FALSE
) {

    count_gene <- paste0("nCount_", gene)
    seurat_obj[[count_gene]] <- seurat_obj[["RNA"]]@counts[gene, ]

    lognorm_gene <- paste0('log1p_cp10k_', gene)
    seurat_obj[[lognorm_gene]] <- seurat_obj[["RNA"]]@data[gene, ]

    # ----------------------------------------------------------------------
    # Figure 1. QC: Count vs. Total Counts per Cell

    for (color in c('sample_name', 'cell_type')) {

        fig <- plot_scatter(
            seurat_obj@meta.data,
            x='nCount_RNA',
            y=paste0("nCount_", gene),
            color=color,
            xlabel='Total Reads',
            ylabel=paste(gene, 'Counts'),
            title=paste(gene, 'Counts vs. Sequencing Depth per Sample'),
            point_size=0.05,
            jitter_height=0.1,
            jitter=TRUE
        )

        if (showfig) { print(fig) }
        savefig(file.path(dirpath, 'qc', paste0('scatter-qc-', color, '-', file_basename, '-', tolower(gene), '.png')),
                height=800, width=1200,
                troubleshooting=troubleshooting)
    }

    # ----------------------------------------------------------------------
    # Figure 2. UMAP

    fig <- FeaturePlot(seurat_obj,
            reduction = "umap", features = lognorm_gene,
            pt.size = 0.4, min.cutoff = 'q10', order = TRUE, label = FALSE) +
        ggtitle( opt[['gene-of-interest']] )
    if (showfig) { print(fig) }
    savefig(file.path(dirpath, paste0('umap-', file_basename, '-', tolower(gene), '.png')),
            height=4000, width=4000, dpi=1600,
            troubleshooting=troubleshooting)

    # ----------------------------------------------------------------------
    # Figure 3. Violin

    withCallingHandlers({
        fig <- plot_violin(seurat_obj,
            cols=c(lognorm_gene), group.by='cell_type',
            threshold_data=NULL, alpha=0.5,
            xlabel='',
            ylabel='LogNorm Expression',
            title=paste(gene, "Expression")
        )
    }, warning = function(w) {
        if ( any(grepl("fewer than two data points", w)) ) {
            invokeRestart("muffleWarning")
        }
    })
    if (showfig) { print(fig) }
    savefig(file.path(dirpath, paste0('violin-', file_basename, '-', tolower(gene), '.png')),
            height=800, width=800,
            troubleshooting=troubleshooting)

    # ----------------------------------------------------------------------
    # Figure 4. Ridge without 0

    withCallingHandlers({
        fig <- suppressMessages(
            # Scale for x is already present.
            # Adding another scale for x, which will replace the existing scale.
            RidgePlot(seurat_obj, lognorm_gene) +
                xlim(5e-5, NA) +
                labs(x='LogNorm Expression\n( log(1+CP10K) )',
                     y=NULL,
                     title=paste(gene, "Expression") )
        )
    }, warning = function(w) {
        if ( any(grepl("rows containing non-finite values", w)) ) {
            invokeRestart("muffleWarning")
        }
    })
    if (showfig) { print(fig) }
    suppressMessages(
        # Picking joint bandwidth of ...
        savefig(file.path(dirpath, paste0('ridge-', file_basename, '-', tolower(gene), '.png')),
                height=800, width=800,
                troubleshooting=troubleshooting)
    )

    # ----------------------------------------------------------------------
    # Figure 5. Bar plot

    cell_counts <- compute_gene_expression(seurat_obj, gene=gene, ident='cell_type')

    value_cols = c('num_cells_neg', 'num_cells_pos')
    cell_counts_long <- cell_counts %>%
        dplyr::select(all_of(c('cell_type', value_cols))) %>%
        tidyr::pivot_longer(cols=value_cols)
    cell_counts_long[lognorm_gene] <- sapply(cell_counts_long['name'], function(x) gsub('num_cells_', '', x))
    
    fig <- plot_bar(
        cell_counts_long[order(-cell_counts_long$value, decreasing = FALSE), ],
        x='cell_type',
        y='value',
        group.by=lognorm_gene,  # gene of interest
        xlabel=NULL,
        ylabel="Number of Cells",
        title=paste0('Number of ', gene, '+ Cells'),
        legend_title=gene
    )
    if (showfig) { print(fig) }   
    savefig(file.path(dirpath, paste0('histogram-cell_type-', file_basename, '-', tolower(gene), '.png')),
            height=800, width=1200, dpi=400, troubleshooting=troubleshooting)
}


#' Draw Predictions QC
#' 
#' @description
#' 
draw_predictions <- function(
    predictions,
    dirpath,
    file_basename='SeuratProject',
    troubleshooting=FALSE,
    showfig=FALSE
) {

    # ----------------------------------------------------------------------
    # Figure 1. Plot CellDex QC Heatmap
    
    fig <- plotScoreHeatmap(predictions)
    if (showfig) { print(fig) }
    savefig(file.path(dirpath, paste0('heatmap-prediction_score-', file_basename, '.png')),
            fig=fig, lib='grid', height=1600, width=2400, dpi=300,
            makedir=TRUE, troubleshooting=troubleshooting)

    # ----------------------------------------------------------------------
    # Figure 2. Number of cells per label

    num_cells_per_label <- table(predictions['labels'])
    num_cells_per_label <- num_cells_per_label[sort.list(num_cells_per_label, decreasing=TRUE)]

    fig <- plot_bar(
        data.frame(num_cells_per_label),
        x='labels',
        y='Freq',
        fill='steelblue',
        group.by=NULL,  # gene of interest
        xlabel=NULL,
        ylabel="Number of Cells",
        title="Histogram of Cell Types"
    )
    if (showfig) { print(fig) }
    savefig(file.path(dirpath, paste0('histogram-cell_type-', file_basename, '.png')),
            height=800, width=1200, dpi=400, troubleshooting=troubleshooting)

}


#' Draw Differential Genes
#' 
draw_differential_genes <- function(
    seurat_obj,
    gene,
    dirpath,
    file_basename='SeuratProject',
    include_pseudo_bulk=FALSE,
    troubleshooting=FALSE,
    showfig=FALSE
) {
    gene_pos <- paste0(tolower(gene), '_pos')
    gene_neg <- paste0(tolower(gene), '_neg')

    df <- compute_gene_labels(seurat_obj, gene=gene, min_reads=1, sep=';')
    orig_ident <- Idents(seurat_obj)
    seurat_obj[['gene_cell_type']] <- df[['gene_cell_type']]
    Idents(seurat_obj) <- df[['gene_cell_type']]
    
    expr_tbl <- compute_gene_expression(seurat_obj, gene=gene, min_reads=1, ident='cell_type')

    # ----------------------------------------------------------------------
    # Figure 1. DEG analysis on single cell

    # only export groups with sufficient cells
    cell_types <- expr_tbl[
        (expr_tbl['num_cells_pos'] >= 50), 'cell_type'
    ]
    for (cell_type in cell_types) {

        markers <- FindMarkers(
            seurat_obj[items_in_a_not_b(rownames(seurat_obj[['RNA']]), gene), ],
            ident.1 = paste(cell_type, snake_to_title_case(gene_pos), sep=';'),
            ident.2 = paste(cell_type, snake_to_title_case(gene_neg), sep=';'),
            logfc.threshold = 0.1,  # note: smaller values take significantly more time to compute
            min.pct = 0.01,
            test.use = "wilcox",  # note: LR gives similar results
            only.pos = FALSE
        )
        p_cutoff <- compute_p_cutoff(
            markers, n_genes=15,
            fc_col='avg_log2FC', pval_col='p_val'
        )

        withCallingHandlers({
            fig <- EnhancedVolcano(
                markers[(rownames(markers)!=gene), ],
                x='avg_log2FC',
                y="p_val",  # p_val_adj scales everything using the bonferroni correction
                lab=rownames(markers),
                title = paste(gene, 'Positive vs. Negative', cell_type),
                # ylab =  bquote(~-Log[10] ~ italic(adjusted P)),
                subtitle = NULL,
                pCutoff = 10^(-p_cutoff),  # default
                FCcutoff = 1
            )
        }, warning = function(w) {
            if ( any(grepl("One or more p-values is 0.", w)) ) {
                invokeRestart("muffleWarning")
            }
        })

        if (showfig) { print(fig) }
        savefig(
            file.path(dirpath,
            paste0('volcano-', title_to_snake_case(cell_type), '-', tolower(gene), '-', file_basename, '.png')),
            height=2000, width=3000, dpi=300, troubleshooting=troubleshooting
        )
    }

    # ----------------------------------------------------------------------
    # Figure 2. DEG analysis on pseudo-bulk data

    if (include_pseudo_bulk) {

        pseudo_bulk <- AggregateExpression(
            seurat_obj[items_in_a_not_b(rownames(seurat_obj[['RNA']]), gene), ],
            assays = "RNA",
            return.seurat = TRUE,
            group.by = c("sample_name", "gene_cell_type"),
            normalization.method = "LogNormalize"  # nothing to do with this
        )

        # this may be problematic if cell labels include ':'
        Idents(pseudo_bulk) <- unname(sapply(
            rownames(pseudo_bulk@meta.data),
            function(x) strsplit(
                stringi::stri_replace_last_regex(str=x, pattern = "_", replacement = ":"),
                split=':'
            )[[1]][2]
        ))

        for (cell_type in cell_types) {

            markers <- FindDESeq2Markers(
                object = pseudo_bulk,
                ident.1 = paste(cell_type, snake_to_title_case(gene_pos), sep=';'),
                ident.2 = paste(cell_type, snake_to_title_case(gene_neg), sep=';'),
                verbose = FALSE
            )
            p_cutoff <- compute_p_cutoff(
                markers, n_genes=15,
                fc_col='log2FoldChange', pval_col='pvalue'  # default
            )

            withCallingHandlers({
                fig <- EnhancedVolcano(
                    markers[(rownames(markers)!=gene), ],
                    x='log2FoldChange',
                    y="pvalue",  # padj scales everything using the bonferroni correction
                    lab=rownames(markers[(rownames(markers)!=gene), ]),
                    title = paste(gene, 'Positive vs. Negative', cell_type),
                    subtitle = NULL,
                    # ylab =  bquote(~-Log[10] ~ italic(adjusted P)),
                    pCutoff = 10^(-p_cutoff),  # default
                    FCcutoff = 1
                )
            }, warning = function(w) {
                if ( any(grepl("One or more p-values is 0.", w)) ) {
                    invokeRestart("muffleWarning")
                }
            })

            if (showfig) { print(fig) }
            savefig(
                file.path(dirpath,
                paste0('pseudo_bulk-', title_to_snake_case(cell_type), '-', tolower(gene), '-', file_basename, '.png')),
                height=2000, width=3000, dpi=300, troubleshooting=troubleshooting
            )
        }
    }
    
    Idents(seurat_obj) <- orig_ident
}
