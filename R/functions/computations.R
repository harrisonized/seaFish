## Compute data from seurat_obj for plotting

# import::here(dplyr, 'count', 'filter')
import::here(rlang, 'sym')
import::here(magrittr, '%>%')
import::here(tidyr, 'pivot_wider')
import::here(file.path(wd, 'R', 'tools', 'df_tools.R'),
    'fillna', 'reset_index', 'rename_columns', .character_only=TRUE)
import::here(file.path(wd, 'R', 'tools', 'text_tools.R'),
    'snake_to_title_case', .character_only=TRUE)

## Functions
## compute_thresholds
## compute_gene_expression
## compute_gene_labels
## compute_p_cutoff


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


#' Compute Gene Expresion
#' 
#' @description
#' Compute number of cells that express the gene and the average expression level
#' of positive cells. [Seurat::AverageExpression()] as a function is extremely convoluted,
#' so I coded the formula to be more explicit.
#' 
#' @examples
#' avg_expr <- setNames(
#'     reset_index(data.frame(t(
#'         log1p(AverageExpression(
#'             seurat_obj_subset[gene, ],
#'             group.by='cell_type')$RNA )
#'    ))),
#'    c('cell_type', 'avg_expression'))
#' 
#' @seealso \href{https://github.com/satijalab/seurat/issues/2423}{ AverageExpression vs FindMarkers: difference in the avg_logFC #2423 }
#' @seealso \href{https://github.com/satijalab/seurat/issues/1719}{ Standard Error or Deviation in AverageExpression #1719 }
#' 
compute_gene_expression <- function(
    seurat_obj,
    gene='Dnase1l1',
    min_reads = 1,
    ident='cell_type'
) {

    is_gene_pos <- paste('is', tolower(gene), 'pos', sep='_')


    # ----------------------------------------------------------------------
    # Cell Counts
    
    tmp <- data.frame(setNames(
        list(  # data
            seurat_obj@meta.data[, 'cell_type'],
            sapply(seurat_obj[["RNA"]]@counts[gene, ] >= min_reads, as.integer)
        ),
        c('cell_type', is_gene_pos)  # colnames
    ))
    df <- as.data.frame(table(tmp)) %>%
        pivot_wider(
            id_cols = 'cell_type',
            names_from=is_gene_pos,
            values_from='Freq')
    df <- rename_columns(
        df,
        columns=c('0'='num_cells_neg',
                  '1'='num_cells_pos'))
    if (!('num_cells_pos' %in% colnames(df))) {
        df[['num_cells_pos']] <- 0
    }

    df['num_cells_total'] <- df['num_cells_neg'] + df['num_cells_pos']
    df['pct_cells_pos'] <- df['num_cells_pos'] / df['num_cells_total']


    # ----------------------------------------------------------------------
    # Average Expression
    
    seurat_obj[[is_gene_pos]] <- as.integer(
        seurat_obj[["RNA"]]@counts[gene, ] >= min_reads
    )
    seurat_obj_subset <- subset(seurat_obj, subset = (!!sym(is_gene_pos) == 1))

    tmp <- data.frame(setNames(
        list(  # data
            seurat_obj_subset@meta.data[, 'cell_type'],
            seurat_obj_subset[["RNA"]]@data[gene, ]
        ),
        c('cell_type', 'expression')  # colnames
    ))
    avg_expr <- setNames(aggregate(
            x=tmp$expression,
            by=list(tmp$cell_type),
            FUN=function(x) log(1+mean(exp(x)-1))
        ), c('cell_type', 'avg_expression')
    )
    # Note about stdev:
    # expression is already log-transformed
    # this is a good approximation as long as expression is not highly skewed
    stdev_expr <- setNames(aggregate(
            x=tmp$expression,
            by=list(tmp$cell_type),
            FUN=function(x) sd(x)
        ), c('cell_type', 'stdev_expression')
    )

    df <- merge(df, avg_expr, by='cell_type', all.x=TRUE, all.y=FALSE)
    df <- merge(df, stdev_expr, by='cell_type', all.x=TRUE, all.y=FALSE)
    df <- fillna(df, cols=c('avg_expression', 'stdev_expression'), val=0)

    # sort
    df <- df[order(df[['num_cells_total']], decreasing = TRUE),]

    return(df)
}


#' Compute Gene Labels
#' 
#' @description
#' Concatenate the cell type and whether it's gene positive or negative into one column.
#' I may revisit this function in the future to make it more general.
#' 
compute_gene_labels <- function(
    seurat_obj,
    gene='Dnase1l1',
    min_reads=1,
    sep=';'
) {

    # colnames
    gene_pos <- paste0(tolower(gene), '_pos')
    gene_neg <- paste0(tolower(gene), '_neg')
    is_gene_pos <- paste('is', tolower(gene), 'pos', sep='_')
    
    # assignment
    df <- data.frame(setNames(
        list(  # data
            seurat_obj@meta.data[, 'cell_type'],
            as.integer(seurat_obj[["RNA"]]@counts[gene, ] >= min_reads),
            ifelse((seurat_obj[["RNA"]]@counts[gene, ] >= min_reads),
                   snake_to_title_case(gene_pos), snake_to_title_case(gene_neg))
        ),
        c('cell_type', is_gene_pos, gene_pos)  # colnames
    ))
    
    df[['gene_cell_type']] <- do.call(
        paste, c(df[c('cell_type', gene_pos)], sep=sep)
    )

    return(df) 
}


#' Compute pCutoff
#'
#' @description
#' Adjust pCutoff so that we get around 15 genes on each side of the volcano.
#' Note: pvalues are base 10, fold changes are base 2
#'
compute_p_cutoff <- function(
    markers,
    n_genes=15,
    fc_col='log2FoldChange',
    pval_col='pvalue',
    p_max=5,
    p_min=-log(0.05, base=10)  # 1.301
) {

    markers[['log2fc_gte_1']] = as.integer(markers[fc_col] >= 1)
    markers[is.na(markers['log2fc_gte_1']), ] <- 0
    num_degs <- sum(markers['log2fc_gte_1'])
    if (num_degs == 0) {
        return(p_min)
    }

    p_cutoff = -log(
        quantile(
            markers[(markers[['log2fc_gte_1']]==1), pval_col],
            probs = min(n_genes / num_degs, 1)
        ), base=10)-0.1

    if (p_cutoff > p_max) {p_cutoff <- p_max}
    if (p_cutoff < p_min) {p_cutoff <- p_min}

    return(p_cutoff)
}
