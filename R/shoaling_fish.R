## Shoals are groups of fish that aggregate together, but are not yet synchronized.
## This script performs data integration, clustering, and unbiased labeling of clusters.
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
import::from(DropletUtils, 'emptyDrops')

import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    'read_10x', 'savefig', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'multiple_replacement', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'text_tools.R'),
    'random_hash', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'single_cell_tools.R'),
    'celldex_switch', .character_only=TRUE)

import::from(file.path(wd, 'R', 'functions', 'computations.R'),
    'compute_thresholds', .character_only=TRUE)
import::from(file.path(wd, 'R', 'functions', 'draw_plots.R'),
    'draw_qc_plots', 'draw_predictions', .character_only=TRUE)
import::from(file.path(wd, 'R', 'functions', 'export_results.R'),
    'export_clustering_results', .character_only=TRUE)

# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-dir"), default='data/ballesteros-2020/input',
                metavar='data/ballesteros-2020/input', type="character",
                help=paste("directory of directories (one up) containing standard 10x files:",
                           "barcodes.tsv, genes.tsv, and matrix.mtx")),
    
    make_option(c("-o", "--output-dir"), default='output',
                metavar='output', type="character",
                help="set the output directory"),

    make_option(c("-j", "--config"), default='config.json',
                metavar='config.json', type="character",
                help="json file containing group information"),
    
    make_option(c("-c", "--celldex"), default='ImmGen',
                metavar='ImmGen', type="character",
                help="Choose from: ['ENCODE', 'HPCA', 'DICE', 'ImmGen', 'Monaco', 'MouseRNAseq', 'Hemato']"),

    make_option(c("-e", "--ensembl"),  default=FALSE,
                metavar="FALSE", action="store_true", type="logical",
                help=paste("Celldex option. Logical scalar indicating whether to convert",
                           "row names to Ensembl IDs. Genes without a mapping to a",
                           "non-duplicated Ensembl ID are discarded.")),

    make_option(c("-m", "--markers"), default=FALSE,
                metavar="FALSE", action="store_true", type="logical",
                help="turn this on to plot heatmap of top markers"),

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

# file io
output_dirname <- opt[['output-dir']]
output_dir <- ifelse(
     grepl('input', opt[['input-dir']]),
     multiple_replacement(opt[['input-dir']], c('input'=output_dirname)),
     file.path(opt[['input-dir']], 'output')
)
figures_dir <- multiple_replacement(
    opt[['input-dir']], c('data'='figures', 'input'=output_dirname)
)


# Start Log
start_time <- Sys.time()
log <- log_open(paste0("shoaling_fish-", random_hash(), '-',
                       strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Parse config

config_file <- file.path(wd, opt[['input-dir']], opt[['config']])
if (file.exists(config_file)) {
    log_print(paste(Sys.time(), 'Reading config...'))
    config <- fromJSON(file=config_file)
    group_names <- names(config)
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

for (group_name in group_names) {
    loop_start_time <- Sys.time()
    log_print(paste(loop_start_time, 'Processing group:', group_name))

    # ----------------------------------------------------------------------
    # Read Data

    sample_names <- config[[group_name]]

    seurat_objs <- new.env()
    threshold_data <- new.env()  
    for (sample_name in sample_names) {
        log_print(paste(Sys.time(), 'Reading...', sample_name))

        tryCatch({

            # ----------------------------------------------------------------------
            # Preprocessing

            raw_mtx <- read_10x(file.path(wd, opt[['input-dir']], sample_name))
            c(num_features, num_cells) %<-% dim(raw_mtx)
            log_print(paste(
                Sys.time(),'Matrix dims: [', num_features, 'features x', num_cells, 'cells]'
            ))

            # Filter empty drops
            tryCatch({
                drop_stats <- emptyDrops(raw_mtx)
                filt_mtx <- raw_mtx[ ,
                    (drop_stats[['FDR']] <= 0.05) &
                    (is.na(drop_stats[['FDR']])==FALSE)
                ]
                c(num_features, num_cells) %<-% dim(filt_mtx)

                log_print(paste(
                    Sys.time(), 'Filtered Matrix dims: [', num_features, 'features x', num_cells, 'cells]'
                ))
            },
            error = function(condition) {
                log_print(paste("Matrix already filtered. Using raw_mtx as filtered."))
                log_print(paste("Message:", conditionMessage(condition)))
                filt_mtx <<- raw_mtx
            })

            # Convert to Seurat Object
            withCallingHandlers({
                tmp_seurat_obj <- CreateSeuratObject(counts = filt_mtx, min.cells = 3)
            }, warning = function(w) {
                if ( any(grepl("Non-unique features", w)) ) {
                    invokeRestart("muffleWarning")
                }
            })
            if (!troubleshooting) {
                rm(raw_mtx)
                rm(filt_mtx)
                gc()
            }

            tmp_seurat_obj$sample_name <- sample_name
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
                group.by='sample_name',
                threshold_data=tmp_threshold_data,
                dirpath=file.path(wd, figures_dir, 'individual', sample_name, 'qc'),
                file_basename=sample_name,
                troubleshooting=troubleshooting,
                showfig=troubleshooting
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
            # Note: Disabled due to overfiltering of false negatives
            # barcodes_to_keep <- rownames(barcode_data[(barcode_data['total'] >= knee), ])
            # tmp_seurat_obj <- tmp_seurat_obj[, barcodes_to_keep]

            # Normalize
            tmp_seurat_obj <- tmp_seurat_obj %>%
                NormalizeData(verbose = FALSE) %>%
                FindVariableFeatures(
                    verbose = FALSE, selection.method = "vst",
                    nfeatures = 2000
                )

            seurat_objs[[sample_name]] <- tmp_seurat_obj

            if (!troubleshooting) {
                rm(tmp_seurat_obj)
                gc()
            }
        },

        # pass
        error = function(condition) {
            log_print(paste("WARNING: SAMPLE", sample_name, "NOT PROCESSED!!!"))
            log_print(paste("Error message: ", conditionMessage(condition)))
            sample_names <<- sample_names[(sample_names != sample_name)]
        })
    }

    seurat_objs <- as.list(seurat_objs)
    threshold_data <- do.call("rbind", as.list(threshold_data))


    # ----------------------------------------------------------------------
    # Integrate Data

    if (length(seurat_objs)>1) {

        log_print(paste(Sys.time(), 'Integrating data...'))

        # integrate
        features_list <- SelectIntegrationFeatures(object.list = seurat_objs, nfeatures = 2000)
        withCallingHandlers({
            anchors <- FindIntegrationAnchors(object.list = seurat_objs, anchor.features = features_list)
        }, warning = function(w) {
            if ( any(grepl("Some cell names are duplicated", w)) ) {
                invokeRestart("muffleWarning")
            }
        })
        seurat_obj <- IntegrateData(anchorset = anchors)  # reduction = "pca" may run faster?
        DefaultAssay(seurat_obj) <- "integrated"
        is_integrated <- TRUE

    } else if (length(seurat_objs)==1) {

        log_print(paste(Sys.time(), 'Processing individual dataset...'))
        seurat_obj <- seurat_objs[[1]]
        group_name <- sample_name
        is_integrated <- FALSE

    } else {

        log_print(paste(Sys.time(), 'No files processed, skipping loop...'))
        next()
    }

    multiplicity <- ifelse(is_integrated, 'integrated', 'individual')
    prefix <- ifelse(is_integrated, 'integrated-', '')


    # ----------------------------------------------------------------------
    # Label Clusters using SingleR
    
    # Run the standard workflow
    ndim = 30  # standard
    seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
    seurat_obj <- RunPCA(seurat_obj, npcs = ndim, verbose = FALSE)  # required for RunUMAP
    seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:ndim)
    seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:ndim)
    seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
    DefaultAssay(seurat_obj) <- "RNA"


    log_print(paste(Sys.time(), 'Labeling clusters...'))

    withCallingHandlers({
        ref_data <- celldex_switch[[opt$celldex]](ensembl=opt[['ensembl']])
    }, warning = function(w) {
        if ( any(grepl("was built under R version", w)) ) {
            invokeRestart("muffleWarning")
        }
    })

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
        dirpath=file.path(wd, figures_dir, multiplicity, group_name, 'labels'),
        file_basename=group_name,
        troubleshooting=troubleshooting,
        showfig=troubleshooting
    )

    # ----------------------------------------------------------------------
    # Save Data

    log_print(paste(Sys.time(), 'Saving Seurat object...'))

    # seurat_obj
    filepath=file.path(wd, output_dir, 'rdata',
        paste0('seurat_obj-', prefix, group_name, '.RData'))
    if (!troubleshooting) {
        if (!dir.exists(dirname(filepath))) { dir.create(dirname(filepath), recursive=TRUE) }
        save(seurat_obj, file=filepath)
    }
    

    # ----------------------------------------------------------------------
    # Export Integrated Results

    if (is_integrated) {

        log_print(paste(Sys.time(), 'Exporting integrated full...'))

        export_clustering_results(
            seurat_obj,
            threshold_data=threshold_data,
            group.by='sample_name',
            sample_name=group_name,
            multiplicity=multiplicity,
            data_dir=output_dir,
            figures_dir=figures_dir,
            file_basename=paste0(prefix, group_name),
            troubleshooting=troubleshooting,
            showfig=troubleshooting
        )

        log_print(paste(Sys.time(), 'Exporting integrated subset...'))

        # Keep significant populations (number of cells > 50)
        num_cells_per_label <- as.data.frame(table(seurat_obj$cell_type))  # value counts
        populations_to_keep <- num_cells_per_label[num_cells_per_label['Freq'] > 50, 'Var1']
        seurat_obj_subset <- seurat_obj[, seurat_obj$cell_type %in% populations_to_keep]

        export_clustering_results(
            seurat_obj_subset,
            group.by='cell_type',
            sample_name=group_name,
            multiplicity=multiplicity,
            data_dir=output_dir,
            figures_dir=figures_dir,
            figures_subdir='subset',
            file_basename=paste0(prefix, group_name, '-subset'),
            troubleshooting=troubleshooting,
            showfig=troubleshooting
        )
    }


    # ----------------------------------------------------------------------
    # Export Individual Results

    log_print(paste(Sys.time(), 'Exporting individual results...'))

    seurat_objs <- SplitObject(seurat_obj, split.by = "sample_name")

    for (sample_name in sample_names) {

        export_clustering_results(
            seurat_objs[[sample_name]],
            sample_name=sample_name,
            multiplicity='individual',
            data_dir=output_dir,
            figures_dir=figures_dir,
            file_basename=sample_name,
            troubleshooting=troubleshooting,
            showfig=troubleshooting
        )
    }
    

    # ----------------------------------------------------------------------
    # Find Top Markers in subset
    # Takes a while, especially if there are many groups
    # See: https://satijalab.org/seurat/articles/pbmc3k_tutorial#finding-differentially-expressed-features-cluster-biomarkers
    # Disabled in normal runthrough to save time


    if (opt[['markers']]) {

        log_print(paste(Sys.time(), 'Identifying top markers...'))

        Idents(seurat_obj_subset) <- seurat_obj_subset$cell_type

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

        # draw heatmap
        seurat_obj_subset <- ScaleData(seurat_obj_subset, verbose = FALSE)  # see: https://github.com/satijalab/seurat/issues/2960
        DoHeatmap(seurat_obj_subset, features = top_markers$gene, label=FALSE, raster=TRUE)
        dirpath <- file.path(wd, figures_dir, multiplicity, group_name, 'labels')
        savefig(file.path(dirpath, paste0('heatmap-top_markers-', group_name, '-subset.png')),
                height=400*length(populations_to_keep),
                width=400*length(populations_to_keep)+400, dpi=400,
                troubleshooting=troubleshooting)
    }


    log_print(paste("Loop completed in:", difftime(Sys.time(), loop_start_time)))
}


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
