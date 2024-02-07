## Code used to analyze scRNAseq data
## Integrates the datasets based on groups set in config.json

## References:
## https://satijalab.org/seurat/archive/v4.3/integration_introduction
## https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
## https://www.10xgenomics.com/analysis-guides/common-considerations-for-quality-control-filters-for-single-cell-rna-seq-data

wd = dirname(this.path::here())  # wd = '~/github/R/seaFish'
suppressPackageStartupMessages(library('Seurat'))
library('zeallot')  # %<-%
library('optparse')
library('logr')
import::from(magrittr, '%>%')
import::from(rjson, 'fromJSON', 'toJSON')
import::from(dplyr, 'group_by', 'ungroup', 'filter', 'top_n', 'slice_head')
import::from(SingleR, 'SingleR')

import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    'read_10x', 'savefig', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'multiple_replacement', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'single_cell_tools.R'),
    'celldex_switch', .character_only=TRUE)
import::from(file.path(wd, 'R', 'functions', 'computations.R'),
    'compute_thresholds', 'compute_cell_counts', .character_only=TRUE)
import::from(file.path(wd, 'R', 'functions', 'draw_plots.R'),
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

    make_option(c("-s", "--slice"), default="",
                metavar="", type="character",
                help="enter comma separated list of indices"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
troubleshooting <- opt[['troubleshooting']]
gene <- opt[["gene-of-interest"]]

output_dir <- multiple_replacement(opt[['input-dir']], c('input'='output'))
figures_dir <- multiple_replacement(opt[['input-dir']], c('data'='figures', 'input'='output'))


# Start Log
start_time <- Sys.time()
log <- log_open(paste0("integrate_data-",
                       strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Parse config

config_file <- file.path(wd, opt[['input-dir']], opt[['config']])
if (file.exists(config_file)) {
    log_print(paste(Sys.time(), 'Reading config...'))
    config <- fromJSON(file=config_file)
} else {
    log_print(paste(Sys.time(), 'Generating config...'))
    group_names <- list.dirs(file.path(wd, opt[['input-dir']]), full.names=FALSE)
    group_names <- group_names[(group_names != '')]  # filter empty
    config <- as.list(setNames(group_names, group_names))
    
    if (!troubleshooting) {
        write(toJSON(config), file=config_file)
    }
}

# subset
if (opt[['slice']] != '') {
    config <- config[ eval(parse( text=paste0('c(', opt[['slice']], ')') )) ]
}


# ----------------------------------------------------------------------
# Main

log_print(paste(Sys.time(), 'Groups found...', paste(names(config), collapse=', ' )))

for (group_name in names(config)) {

    loop_start_time <- Sys.time()
    log_print(paste(loop_start_time, 'Processing group:', group_name))

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

            # For the filters, thresholds are adjusted automatically
            # Lower bound: low quality cells
            # Upper bound: potential doublets
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
    
    # ----------------------------------------------------------------------
    # Integrate Data

    if (length(expr_mtxs)>1) {

        log_print(paste(Sys.time(), 'Integrating data...'))

        threshold_data <- do.call("rbind", as.list(threshold_data))

        # integrate
        features_list <- SelectIntegrationFeatures(object.list = expr_mtxs, nfeatures = 2000)
        anchors <- FindIntegrationAnchors(object.list = expr_mtxs, anchor.features = features_list)
        seurat_obj <- IntegrateData(anchorset = anchors)  # reduction = "pca" may run faster?
        DefaultAssay(seurat_obj) <- "integrated"

        draw_qc_plots(
            seurat_obj,
            dirpath=file.path(wd, figures_dir, "integrated", group_name, 'qc'),
            prefix='integrated-',
            sample_name=group_name,
            group.by='sample_name',
            threshold_data=threshold_data,
            troubleshooting=troubleshooting,
            showfig=TRUE
        )

        integrated_label <- "integrated"
        rm(tmp_seurat_obj)
        rm(expr_mtxs)

    } else if (length(expr_mtxs)==1) {

        log_print(paste(Sys.time(), 'Processing individual dataset...'))

        seurat_obj <- expr_mtxs[[1]]
        integrated_label <- "individual"
        group_name <- sample_name

        rm(tmp_seurat_obj)
        rm(expr_mtxs)

    } else {

        log_print(paste(Sys.time(), 'No files processed, skipping loop...'))

        next()
    }


    # ----------------------------------------------------------------------
    # Label Clusters using SingleR
    
    # Run the standard workflow
    seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
    seurat_obj <- RunPCA(seurat_obj, npcs = opt[['ndim']], verbose = FALSE)  # required for RunUMAP
    seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:opt[['ndim']])
    seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:opt[['ndim']])
    seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

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

    draw_predictions(
        predictions,
        dirpath=file.path(wd, figures_dir, integrated_label, group_name, 'labels'),
        group_name=group_name,
        troubleshooting=troubleshooting,
        showfig=TRUE
    )


    # ----------------------------------------------------------------------
    # Save Data

    log_print(paste(Sys.time(), 'Saving Seurat object...'))

    prefix <- ifelse(integrated_label=='integrated', 'integrated-', '')

    # seurat_obj
    filepath=file.path(wd, output_dir, 'rdata',
        paste0('seurat_obj-', prefix, group_name, '.RData'))
    if (!troubleshooting) {
        if (!dir.exists(dirname(filepath))) { dir.create(dirname(filepath), recursive=TRUE) }
        save(seurat_obj, file=filepath)
    }

    # expression csv
    cell_counts <- compute_cell_counts(seurat_obj, gene=gene, ident='cell_type')

    filepath=file.path(wd, output_dir, 'expression', tolower(gene), integrated_label,
        paste0('cell_type-', prefix, group_name, '-', tolower(gene), '.csv'))
    if (!troubleshooting) {
        if ( !dir.exists(dirname(filepath)) ) { dir.create(dirname(filepath), recursive=TRUE) }
        write.table(cell_counts, file = filepath, row.names = FALSE, sep = ',')
    }


    # ----------------------------------------------------------------------
    # Plot Results

    log_print(paste(Sys.time(), 'Plotting full results...'))

    # filepath=file.path(wd, output_dir, 'rdata', paste0('integrated-', group_name, '.RData'))
    # import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    #     'load_rdata', .character_only=TRUE) 
    # seurat_obj <- load_rdata(filepath=filepath)

    draw_clusters(
        seurat_obj,
        dirpath=file.path(wd, figures_dir, integrated_label, group_name, 'expression'),
        prefix=ifelse(integrated_label=='integrated', 'integrated-', ''),
        group_name=group_name,
        troubleshooting=troubleshooting,
        showfig=TRUE
    )

    draw_gene_of_interest(
        seurat_obj,
        gene=gene,
        dirpath=file.path(wd, figures_dir, integrated_label, group_name, 'expression', tolower(gene)),
        prefix=ifelse(integrated_label=='integrated', 'integrated-', ''),
        group_name=group_name,
        troubleshooting=troubleshooting,
        showfig=TRUE
    )


    # ----------------------------------------------------------------------
    # Plot Subset

    log_print(paste(Sys.time(), 'Plotting subset...'))

    # Keep significant populations (number of cells > 50)
    num_cells_per_label <- as.data.frame(table(seurat_obj$cell_type))  # value counts
    populations_to_keep <- num_cells_per_label[num_cells_per_label['Freq'] > 50, 'Var1']
    seurat_obj_subset <- seurat_obj[, seurat_obj$cell_type %in% populations_to_keep]

    draw_qc_plots(
        seurat_obj_subset,
        dirpath=file.path(wd, figures_dir, integrated_label, group_name, 'qc', 'subset'),
        prefix=ifelse(integrated_label=='integrated', 'integrated-', ''),
        suffix='-subset',
        sample_name=group_name,
        group.by='cell_type',
        threshold_data=NULL,
        troubleshooting=troubleshooting,
        showfig=TRUE
    )

    draw_clusters(
        seurat_obj_subset,
        dirpath=file.path(wd, figures_dir, integrated_label, group_name, 'expression'),
        prefix=ifelse(integrated_label=='integrated', 'integrated-', ''),
        suffix='-subset',
        celldex_dataset=opt[['celldex']],
        group_name=group_name,
        # split.by='sample_name',
        troubleshooting=troubleshooting,
        showfig=TRUE
    )

    draw_gene_of_interest(
        seurat_obj_subset,
        gene=gene,
        dirpath=file.path(wd, figures_dir, integrated_label, group_name, 'expression', tolower(gene), 'subset'),
        prefix=ifelse(integrated_label=='integrated', 'integrated-', ''),
        suffix='-subset',
        group_name=group_name,
        troubleshooting=troubleshooting,
        showfig=TRUE
    )


    # ----------------------------------------------------------------------
    # Find Top Markers in subset

    log_print(paste(Sys.time(), 'Identifying top markers...'))

    Idents(seurat_obj_subset) <- seurat_obj_subset$cell_type

    # Takes a while, especially if there are many groups
    # see: https://satijalab.org/seurat/articles/pbmc3k_tutorial#finding-differentially-expressed-features-cluster-biomarkers
    markers <- FindAllMarkers(
        seurat_obj_subset,
        logfc.threshold = 0.25,  # increased to improve speed; we only need the top 12 markers
        min.pct = 0.1,
        test.use = "wilcox",
        only.pos = TRUE
    )
    top_markers <- markers %>%
        group_by(cluster) %>% filter(avg_log2FC > 1) %>%
        slice_head(n = 12) %>% ungroup()

    # Plot
    seurat_obj_subset <- ScaleData(seurat_obj_subset, verbose = FALSE)  # see: https://github.com/satijalab/seurat/issues/2960
    DoHeatmap(seurat_obj_subset, features = top_markers$gene, label=FALSE, raster=TRUE)
    dirpath <- file.path(wd, figures_dir, integrated_label, group_name, 'labels')
    savefig(file.path(dirpath, paste0('heatmap-top_markers-', group_name, '-subset.png')),
            height=400*length(populations_to_keep), width=1200, dpi=400,
            troubleshooting=troubleshooting)


    log_print(paste("Loop completed in:", difftime(Sys.time(), loop_start_time)))
}


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
