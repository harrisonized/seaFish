import::here(file.path(wd, 'R', 'functions', 'computations.R'),
    'compute_cell_counts', .character_only=TRUE)
import::here(file.path(wd, 'R', 'functions', 'draw_plots.R'),
    'draw_qc_plots', 'draw_clusters', 'draw_gene_of_interest',
    'draw_differential_genes', .character_only=TRUE)

## Functions
## export_clustering_results
## export_gene_of_interest


#' Export Clustering Results
#' 
export_clustering_results <- function(
    seurat_obj,
    threshold_data=NULL,
    group.by='sample_name',
    sample_name='SeuratProject',
    multiplicity='integrated',
    data_dir='data/output',
    figures_dir='figures/output',
    figures_subdir='',
    file_basename='SeuratProject',
    include_heatmap=TRUE,
    troubleshooting=FALSE,
    showfig=FALSE
) {

    draw_qc_plots(
        seurat_obj,
        group.by=group.by,
        threshold_data=threshold_data,
        dirpath=file.path(wd, figures_dir, multiplicity, sample_name, 'qc', figures_subdir),
        file_basename=file_basename,
        troubleshooting=troubleshooting,
        showfig=showfig
    )

    draw_clusters(
        seurat_obj,
        dirpath=file.path(wd, figures_dir, multiplicity, sample_name, 'expression'),
        file_basename=file_basename,
        title=sample_name,
        include_heatmap=include_heatmap,
        troubleshooting=troubleshooting,
        showfig=showfig
    )
}


#' Export Gene of Interest
#' 
export_gene_of_interest <- function(
    seurat_obj,
    gene='Dnase1l1',
    sample_name='SeuratProject',
    multiplicity='integrated',
    data_dir='data/output',
    figures_dir='figures/output',
    figures_subdir='',
    file_basename='SeuratProject',
    include_volcano=TRUE,
    troubleshooting=FALSE,
    showfig=FALSE
) {

    # gene-of-interest by cell-type csv
    cell_counts <- compute_cell_counts(seurat_obj, gene=gene, ident='cell_type')
    filepath=file.path(
        wd, data_dir, 'expression', tolower(gene), multiplicity,
        paste0('cell_type-', file_basename, '-', tolower(gene), '.csv')
    )
    if (!troubleshooting) {
        if ( !dir.exists(dirname(filepath)) ) { dir.create(dirname(filepath), recursive=TRUE) }
        write.table(cell_counts, file = filepath, row.names = FALSE, sep = ',')
    }

    draw_gene_of_interest(
        seurat_obj,
        gene=gene,
        dirpath=file.path(wd, figures_dir, multiplicity, sample_name, 'expression', tolower(gene), figures_subdir),
        file_basename=file_basename,
        troubleshooting=troubleshooting,
        showfig=showfig
    )

    if (include_volcano) {
        draw_differential_genes(
            seurat_obj,
            gene=gene,
            dirpath=file.path(wd, figures_dir, multiplicity, sample_name, 'expression', tolower(gene), figures_subdir, 'degs'),
            file_basename=file_basename,
            troubleshooting=troubleshooting,
            showfig=troubleshooting
        )
    }
}
