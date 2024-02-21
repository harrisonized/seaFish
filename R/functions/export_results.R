import::here(file.path(wd, 'R', 'functions', 'computations.R'),
    'compute_cell_counts', .character_only=TRUE)
import::here(file.path(wd, 'R', 'functions', 'draw_plots.R'),
    'draw_qc_plots', 'draw_clusters', 'draw_gene_of_interest',
    .character_only=TRUE)

## Functions
## export_analysis_results


#' Export Analysis Results
#' 
#' @description
#' Given a seurat_obj, exports three groups of figures:
#' qc plots, umaps of labeled clusters, and gene_of_interest figures
#' Can also export the raw counts
#' 
export_analysis_results <- function(
    seurat_obj,
    gene='Dnase1l1',
    threshold_data=NULL,
    group.by='sample_name',
    data_dir='data/output',
    figures_dir='figures/output',
    multiplicity='integrated',
    sample_name='SeuratProject',
    figures_subdir='',
    file_basename='SeuratProject',
    plot_qc=TRUE,
    export_counts=TRUE,
    troubleshooting=FALSE,
    showfig=FALSE
) {

    # gene-of-interest by cell-type csv
    if (export_counts) {
        cell_counts <- compute_cell_counts(seurat_obj, gene=gene, ident='cell_type')
        filepath=file.path(wd, data_dir, 'expression', tolower(gene), multiplicity,
            paste0('cell_type-', file_basename, '-', tolower(gene), '.csv'))
        if (!troubleshooting) {
            if ( !dir.exists(dirname(filepath)) ) { dir.create(dirname(filepath), recursive=TRUE) }
            write.table(cell_counts, file = filepath, row.names = FALSE, sep = ',')
        }
    }

    if (plot_qc) {
        draw_qc_plots(
            seurat_obj,
            group.by=group.by,
            threshold_data=threshold_data,
            dirpath=file.path(wd, figures_dir, multiplicity, sample_name, 'qc', figures_subdir),
            file_basename=file_basename,
            troubleshooting=troubleshooting,
            showfig=showfig
        )
    }

    draw_clusters(
        seurat_obj,
        dirpath=file.path(wd, figures_dir, multiplicity, sample_name, 'expression'),
        file_basename=file_basename,
        title=sample_name,
        troubleshooting=troubleshooting,
        showfig=showfig
    )

    draw_gene_of_interest(
        seurat_obj,
        gene=gene,
        dirpath=file.path(wd, figures_dir, multiplicity, sample_name, 'expression', tolower(gene), figures_subdir),
        file_basename=file_basename,
        troubleshooting=troubleshooting,
        showfig=showfig
    )
}
