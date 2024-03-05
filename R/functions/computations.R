## Compute data from seurat_obj for plotting

# import::here(dplyr, 'count', 'filter')
import::here(magrittr, '%>%')
import::here(tidyr, 'pivot_wider')
import::here(file.path(wd, 'R', 'tools', 'df_tools.R'),
    'fillna', .character_only=TRUE)
import::here(file.path(wd, 'R', 'tools', 'text_tools.R'),
    'snake_to_title_case', .character_only=TRUE)

## Functions
## compute_thresholds
## compute_cell_counts
## compute_gene_labels


#' Compute Thresholds
#'
compute_thresholds <- function(
    seurat_obj,
    sample_name = 'SeuratProject'
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

    gene_pos <- paste0(tolower(gene), '_pos')
    seurat_obj[[gene_pos]] <- ifelse(
        (seurat_obj[["RNA"]]@counts[gene, ] > 0),
        'num_cells_pos', 'num_cells_neg'
    )
    cell_counts <- as.data.frame(table(
        seurat_obj@meta.data[, c('cell_type', gene_pos)])) %>%
        pivot_wider(
            id_cols = 'cell_type',
            names_from=gene_pos,
            values_from='Freq') %>%
        as.data.frame()

    cell_counts['num_cells_total'] <- cell_counts['num_cells_neg'] + cell_counts['num_cells_pos']
    cell_counts['pct_cells_pos'] <- cell_counts['num_cells_pos'] / cell_counts['num_cells_total']
    
    # sort
    cell_counts <- cell_counts[order(cell_counts[['num_cells_total']], decreasing = TRUE),]
    return(cell_counts)
}


#' Compute Gene Labels
#' 
compute_gene_labels <- function(
    seurat_obj,
    gene='Dnase1l1',
    sep=';'
) {

    # colnames
    gene_pos <- paste0(tolower(gene), '_pos')
    gene_neg <- paste0(tolower(gene), '_neg')
    is_gene_pos <- paste('is', tolower(gene), 'pos', sep='_')
    
    # assignment
    df <- data.frame(setNames(
        list(
            seurat_obj@meta.data[, 'cell_type'],
            as.integer(seurat_obj[["RNA"]]@counts[gene, ] > 0),
            ifelse((seurat_obj[["RNA"]]@counts[gene, ] > 0),
                   snake_to_title_case(gene_pos), snake_to_title_case(gene_neg))
        ),  # data
        c('cell_type', is_gene_pos, gene_pos)  # colnames
    ))
    
    df[['gene_cell_type']] <- do.call(
        paste, c(df[c('cell_type', gene_pos)], sep=sep)
    )

    return(df) 
}
