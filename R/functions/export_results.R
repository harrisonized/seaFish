import::here(file.path(wd, 'R', 'functions', 'computations.R'),
    'compute_cell_counts', .character_only=TRUE)
import::here(file.path(wd, 'R', 'functions', 'draw_plots.R'),
    'draw_qc_plots', 'draw_clusters', 'draw_gene_of_interest',
    .character_only=TRUE)

## Functions
## export_results


#' Export Results
#' 
#' @description
#' Given a seurat_obj, exports all required data and figures
#' 
export_results <- function(
    seurat_obj,
    threshold_data=NULL,
    group_name='orig.ident',
    gene='Dnase1l1',
    output_dir='data/output',
    figures_dir='figures/output',
    subdir='',
    prefix=NULL,
    suffix=NULL,
    integrated_label='integrated',
    group.by='sample_name',
    export_counts=TRUE,
    plot_qc=TRUE,
    troubleshooting=FALSE
) {

    if (is.null(prefix)) {
        prefix <- ifelse(integrated_label=='integrated', 'integrated-', '')
    }

    if (export_counts) {
         # expression csv
        cell_counts <- compute_cell_counts(seurat_obj, gene=gene, ident='cell_type')
        filepath=file.path(wd, output_dir, 'expression', tolower(gene), integrated_label,
            paste0('cell_type-', prefix, group_name, '-', tolower(gene), '.csv'))
        if (!troubleshooting) {
            if ( !dir.exists(dirname(filepath)) ) { dir.create(dirname(filepath), recursive=TRUE) }
            write.table(cell_counts, file = filepath, row.names = FALSE, sep = ',')
        }
    }

    if (plot_qc) {
        draw_qc_plots(
            seurat_obj,
            dirpath=file.path(wd, figures_dir, integrated_label, group_name, 'qc', subdir),
            prefix=prefix,
            suffix=suffix,
            sample_name=group_name,
            group.by=group.by,
            threshold_data=threshold_data,
            troubleshooting=troubleshooting,
            showfig=TRUE
        )
    }

    draw_clusters(
        seurat_obj,
        dirpath=file.path(wd, figures_dir, integrated_label, group_name, 'expression'),
        prefix=prefix,
        suffix=suffix,
        group_name=group_name,
        troubleshooting=troubleshooting,
        showfig=TRUE
    )

    draw_gene_of_interest(
        seurat_obj,
        gene=gene,
        dirpath=file.path(wd, figures_dir, integrated_label, group_name, 'expression', tolower(gene), subdir),
        prefix=prefix,
        suffix=suffix,
        group_name=group_name,
        troubleshooting=troubleshooting,
        showfig=TRUE
    )
}
