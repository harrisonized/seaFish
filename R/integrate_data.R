## Code used to analyze scRNAseq data
## Integrates the datasets based on the information in config.json

## This script reads in the datasets individually
## During the preprocessesing step, filter first, ten normalize
## For the filters, thresholds are adjusted automatically
## Lower bound: low quality cells
## Upper bound: potential doublets

## References:
## https://satijalab.org/seurat/archive/v4.3/integration_introduction
## https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
## https://www.10xgenomics.com/analysis-guides/common-considerations-for-quality-control-filters-for-single-cell-rna-seq-data


wd = dirname(this.path::here())  # wd = '~/github/R/seaFish'
suppressPackageStartupMessages(library('Seurat'))
library("ggplot2")
library('rjson')
library('zeallot')
library('optparse')
library('logr')
import::from(magrittr, '%>%')
import::from(tidyr, 'pivot_longer')
import::from(dplyr, 'group_by', 'ungroup', 'filter', 'top_n', 'slice_head', 'select')
import::from(DropletUtils, 'barcodeRanks')
import::from(SingleR, 'SingleR', 'plotScoreHeatmap')

import::from(file.path(wd, 'R', 'tools', 'df_tools.R'),
    'fillna', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    'read_10x', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'multiple_replacement', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'plotting.R'),
    'savefig', 'plot_violin', 'plot_waterfall', 'plot_scatter', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'single_cell_tools.R'),
    'celldex_switch', .character_only=TRUE)
import::from(file.path(wd, 'R', 'functions', 'quality_control.R'),
    'compute_thresholds', 'compute_cell_counts',
    'draw_qc_plots', 'draw_predictions', 'draw_clusters', 'draw_gene_of_interest',
    .character_only=TRUE)


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-dir"),
                default='data/ballesteros-2020/input',
                metavar='data/ballesteros-2020/input',
                type="character",
                help="directory of directories (one up) containing standard 10x files: barcodes.tsv, genes.tsv, and matrix.mtx"),
    
    make_option(c("-o", "--output-dir"),
                default='data/ballesteros-2020/output',
                metavar='data/ballesteros-2020/output',
                type="character",
                help="set the output directory"),

    make_option(c("-j", "--config"), default='config.json',
                metavar='config.json', type="character",
                help="json file containing group information"),
    
    make_option(c("-d", "--ndim"), default=30,
                metavar="30", type="integer",
                help="choose the number of dimensions for the data integration"),

    make_option(c("-c", "--celldex"),
                default='ImmGen', metavar='ImmGen', type="character",
                help="Choose from: ['ENCODE', 'HPCA', 'DICE', 'ImmGen', 'Monaco', 'MouseRNAseq', 'Hemato']"),

    make_option(c("-e", "--ensembl"),  default=FALSE,
                metavar="FALSE", action="store_true", type="logical",
                help="Used when querying celldex"),

    make_option(c("-g", "--gene-of-interest"), default="Dnase1l1",
                metavar="Dnase1l1", type="character",
                help="choose a gene"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
troubleshooting <- opt[['troubleshooting']]
output_dir <- multiple_replacement(opt[['input-dir']], c('input'='output'))
figures_dir <- multiple_replacement(opt[['input-dir']], c('data'='figures', 'input'='output'))
gene <- opt[["gene-of-interest"]]

# Start Log
start_time <- Sys.time()
log <- log_open(paste0("integrate_data-",
                       strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Main

config <- fromJSON(file=file.path(wd, opt[['input-dir']], opt[['config']]))

for (group_name in names(config)) {

    loop_start_time <- Sys.time()
    log_print(paste(loop_start_time, 'Reading Seurat objects...'))

    sample_names <- config[[group_name]]
    expr_mtxs <- new.env()
    threshold_data <- new.env()


    # ----------------------------------------------------------------------
    # Read Data
    
    for (sample_name in sample_names) {
        log_print(paste(Sys.time(), 'Reading...', sample_name))

        tryCatch({

            # ----------------------------------------------------------------------
            # Preprocessing

            expr_mtx <- read_10x(file.path(wd, opt[['input-dir']], sample_name))
            tmp_seurat_obj <- CreateSeuratObject(counts = expr_mtx, min.cells = 3)
            tmp_seurat_obj$sample_name <- sample_name  # consider using orig.ident
            tmp_seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
                tmp_seurat_obj, pattern = "^mt-"
            )
            c(tmp_threshold_data, upper_rna_count,
              upper_rna_feature, upper_pct_mt) %<-% compute_thresholds(
                  tmp_seurat_obj,
                  sample_name=sample_name
            )
            threshold_data[[sample_name]] <- tmp_threshold_data

            draw_qc_plots(
                tmp_seurat_obj,
                dirpath=file.path(wd, figures_dir, 'individual', sample_name, 'qc'),
                sample_name=sample_name,
                group.by='sample_name',
                threshold_data=tmp_threshold_data,
                troubleshooting=troubleshooting,
                showfig=TRUE
            )

            # ----------------------------------------------------------------------
            # Filter and Normalize

            tmp_seurat_obj <- subset(tmp_seurat_obj,
                subset = ((nCount_RNA > 1000) & (nCount_RNA <= upper_rna_count) &
                          (nFeature_RNA > 300) & (nFeature_RNA <= upper_rna_feature) &
                          (percent.mt <= upper_pct_mt))
            )

            # Filter genes expressed by less than 3 cells
            # See: https://matthieuxmoreau.github.io/EarlyPallialNeurogenesis/html-Reports/Quality_Control.html
            num_cells_per_gene <- rowSums(tmp_seurat_obj[['RNA']]@counts > 0)
            genes_to_keep <- names(num_cells_per_gene[num_cells_per_gene >= 3])
            tmp_seurat_obj <- tmp_seurat_obj[genes_to_keep, ]


            # Filter barcodes below knee point
            # Note: Disabled due to the high amount of false negatives
            # barcodes_to_keep <- rownames(barcode_data[(barcode_data['total'] >= knee), ])
            # tmp_seurat_obj <- tmp_seurat_obj[, barcodes_to_keep]

            # Normalize
            tmp_seurat_obj <- tmp_seurat_obj %>%
                NormalizeData(verbose = FALSE) %>%
                FindVariableFeatures(
                    verbose = FALSE, selection.method = "vst",
                    nfeatures = 2000
                )

            expr_mtxs[[sample_name]] <- tmp_seurat_obj
        },

        # pass
        error = function(condition) {
            log_print(paste("WARNING: SAMPLE", sample_name, "NOT PROCESSED!!!"))
            log_print(paste("Error message: ", conditionMessage(condition)))
            NULL
        })
    }

    expr_mtxs <- as.list(expr_mtxs)
    threshold_data <- do.call("rbind", as.list(threshold_data))


    # ----------------------------------------------------------------------
    # Integrate Data

    log_print(paste(Sys.time(), 'Integrating data...'))

    features_list <- SelectIntegrationFeatures(object.list = expr_mtxs, nfeatures = 2000)
    anchors <- FindIntegrationAnchors(object.list = expr_mtxs, anchor.features = features_list)
    seurat_obj <- IntegrateData(anchorset = anchors)
    DefaultAssay(seurat_obj) <- "integrated"

    draw_qc_plots(
        seurat_obj,
        dirpath=file.path(wd, figures_dir, 'integrated', group_name, 'qc'),
        sample_name=group_name,
        group.by='sample_name',
        threshold_data=threshold_data,
        troubleshooting=troubleshooting,
        showfig=TRUE
    )

    # Run the standard workflow
    seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
    seurat_obj <- RunPCA(seurat_obj, npcs = opt[['ndim']], verbose = FALSE)  # required for RunUMAP
    seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:opt[['ndim']])
    seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:opt[['ndim']])
    seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)


    # ----------------------------------------------------------------------
    # Label Clusters using SingleR
    
    log_print(paste(Sys.time(), 'Labeling clusters...'))

    DefaultAssay(seurat_obj) <- "RNA"
    ref_data <- celldex_switch[[opt$celldex]](ensembl=opt[['ensembl']])
    predictions <- SingleR(
        test=as.SingleCellExperiment(seurat_obj),
        assay.type.test=1,
        ref=ref_data,
        labels=ref_data[['label.main']]
    )
    seurat_obj$cell_type <- predictions[['labels']]
    Idents(seurat_obj) <- seurat_obj$cell_type

    # export predictions as RData
    filepath=file.path(wd, output_dir, 'rdata', paste0('integrated-', group_name, '.RData'))
    if (!troubleshooting) {
        if (!dir.exists(dirname(filepath))) { dir.create(dirname(filepath), recursive=TRUE) }
        save(seurat_obj, file=filepath)
    }

    # import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    #     'load_rdata', .character_only=TRUE) 
    # seurat_obj <- load_rdata(
    #     filepath=file.path(wd, output_dir, 'integrated', paste0(group_name,'.RData'))
    # )


    # ----------------------------------------------------------------------
    # Plot

    plot_scatter(seurat_obj, group.by='cell_type')
    if (!troubleshooting) {
        savefig(file.path(wd, figures_dir, 'integrated', group_name, 'qc',
                          paste0('scatter-features_vs_counts-cell_type-', group_name, '.png')),
        troubleshooting=troubleshooting)
    }

    draw_predictions(
        predictions,
        dirpath=file.path(wd, figures_dir, 'integrated', group_name, 'clustering'),
        group_name=group_name,
        troubleshooting=troubleshooting,
        showfig=TRUE
    )


    # ----------------------------------------------------------------------
    # Gene of Interest

    draw_clusters(
        seurat_obj,
        dirpath=file.path(wd, figures_dir, 'integrated', group_name, 'expression'),
        group_name=group_name,
        troubleshooting=troubleshooting,
        showfig=TRUE
    )

    cell_counts <- compute_cell_counts(seurat_obj, gene=gene, ident='cell_type')

    # save as csv
    filepath=file.path(wd, output_dir, 'expression',
                       paste0(group_name, '-cell_type-', tolower(gene), '.csv'))
    if (!troubleshooting) {
        if ( !dir.exists(dirname(filepath)) ) { dir.create(dirname(filepath), recursive=TRUE) }
        write.table(cell_counts, file = filepath, row.names = FALSE, sep = ',')
    }

    draw_gene_of_interest(
        seurat_obj,
        gene=gene,
        dirpath=file.path(wd, figures_dir, 'integrated', group_name, 'expression'),
        group_name=group_name,
        troubleshooting=troubleshooting,
        showfig=TRUE
    )

    value_cols = c('num_cells_neg', 'num_cells_pos')
    cell_counts_long <- cell_counts %>%
        select(all_of(c('cell_type', value_cols))) %>%
        pivot_longer(cols=value_cols)
    cell_counts_long[gene] <- sapply(cell_counts_long['name'], function(x) gsub('num_cells_', '', x))
    
    ggplot(data=cell_counts_long[order(-cell_counts_long$value, decreasing = FALSE), ],
           aes( x=reorder(.data[['cell_type']], .data[['value']], decreasing=TRUE),
                y=.data[['value']],
                fill=.data[[gene]] )
        ) +
        geom_bar(stat="identity") +
        labs(title = paste0('Number of ', gene, '+ Cells'),
             x = NULL,
             y = "Number of Cells") +
        scale_x_discrete(guide = guide_axis(angle = 45)) +
        theme(legend.position = "bottom")
    
    if (!troubleshooting) {
        savefig(file.path(wd, figures_dir, 'integrated', group_name, 'expression', tolower(gene),
                              paste0('histogram-cell_type-', group_name, '-', tolower(gene), '.png')),
                height=800, width=1200, dpi=400, troubleshooting=troubleshooting)
    }


    # ----------------------------------------------------------------------
    # Plot subset

    # Keep significant populations (percent cells > 3%)
    num_cells_per_label <- table(predictions['labels'])  # value counts
    num_cells_per_label <- num_cells_per_label[sort.list(num_cells_per_label, decreasing=TRUE)]  # sort
    pct_cells_per_label <- num_cells_per_label/sum(num_cells_per_label)
    populations_to_keep <- names(num_cells_per_label[(pct_cells_per_label > 0.03)])
    seurat_obj_subset <- seurat_obj[, seurat_obj$cell_type %in% populations_to_keep]

    draw_qc_plots(
        seurat_obj_subset,
        dirpath=file.path(wd, figures_dir, 'integrated', group_name, 'qc', 'subset'),
        prefix='subset-',
        sample_name=group_name,
        group.by='cell_type',
        threshold_data=NULL,
        troubleshooting=troubleshooting,
        showfig=TRUE
    )

    draw_clusters(
        seurat_obj_subset,
        dirpath=file.path(wd, figures_dir, 'integrated', group_name, 'expression'),
        prefix='subset-',
        group_name=group_name,
        split.by='sample_name',
        troubleshooting=troubleshooting,
        showfig=TRUE
    )

    # gene of interest
    seurat_obj_subset[[gene]] <- seurat_obj_subset[["RNA"]]@data[gene, ]
    plot_violin(seurat_obj_subset,
        cols=c(gene), group.by='cell_type',
        threshold_data=NULL, alpha=0.5)
    if (!troubleshooting) {
        savefig(file.path(wd, figures_dir, 'integrated', group_name, 'expression', tolower(gene),
                          paste0('violin-integrated-', group_name, '-subset-', tolower(gene), '.png')),
                height=800, width=800,
                troubleshooting=troubleshooting)
    }


    # ----------------------------------------------------------------------
    # Find Top Markers (Filtered data only)

    Idents(seurat_obj_subset) <- seurat_obj_subset$cell_type

    # see: https://satijalab.org/seurat/articles/pbmc3k_tutorial#finding-differentially-expressed-features-cluster-biomarkers
    markers <- FindAllMarkers(
        seurat_obj_subset,
        logfc.threshold = 0.25,  # increased to improve speed; we only need the top 12 markers
        min.pct = 0.1,
        test.use = "wilcox", only.pos = TRUE)
    top_markers <- markers %>%
        group_by(cluster) %>% filter(avg_log2FC > 1) %>%
        slice_head(n = 12) %>% ungroup()

    # Plot heatmap
    seurat_obj_subset <- ScaleData(seurat_obj_subset, verbose = FALSE)  # see: https://github.com/satijalab/seurat/issues/2960
    DoHeatmap(seurat_obj_subset, features = top_markers$gene, label=FALSE, raster=TRUE)
    if (!troubleshooting) {
        savefig(file.path(wd, figures_dir, 'integrated', group_name, 'clustering',
                     paste0('heatmap-top_markers-', group_name, '.png')),
                height=400*length(populations_to_keep), width=1200, dpi=400,
                troubleshooting=troubleshooting)
    }

    log_print(paste("Loop completed in:", difftime(Sys.time(), loop_start_time)))
}

end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
