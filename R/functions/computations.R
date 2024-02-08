## Compute data from seurat_obj for plotting

import::here(dplyr, 'count')
import::here(file.path(wd, 'R', 'tools', 'df_tools.R'),
    'fillna', .character_only=TRUE)

## Functions
## compute_thresholds
## compute_cell_counts


#' Compute Thresholds
#'
compute_thresholds <- function(
    seurat_obj,
    sample_name = 'SeruatProject'
) {

    thresholds <- new.env()  # upper thresholds
    for (col in c('nCount_RNA', 'nFeature_RNA', 'percent.mt')) {
        thresholds[[col]]  <- median( unlist(seurat_obj[[col]]) ) +
            3*sd( unlist(seurat_obj[[col]]) )
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

    num_total_cells <- dplyr::count(seurat_obj@meta.data, .data[[ident]], name='num_cells')
    num_pos_cells <- dplyr::count(
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
