## Code used to analyze scRNAseq data
## Integrates the datasets based on the information in config.json
## Based on: https://satijalab.org/seurat/articles/integration_introduction.html

wd = dirname(this.path::here())  # wd = '~/github/R/seaFish'
suppressPackageStartupMessages(library('Seurat'))
library('rjson')
library('optparse')
library('logr')
import::from('magrittr', '%>%', .character_only=TRUE)
import::from('DropletUtils', 'write10xCounts', .character_only=TRUE)
import::from('SingleR', 'SingleR', .character_only=TRUE)
import::from('ggplot2', 'ggsave', 'ggtitle', .character_only=TRUE)
import::from(file.path(wd, 'R', 'utils', 'file_io.R'),
    'read_10x', .character_only=TRUE)
import::here(file.path(wd, 'R', 'utils', 'list_tools.R'),
    'multiple_replacement', .character_only=TRUE)
source(file.path(wd, 'R', 'functions', 'bio.R'))  # celldex_switch, run_standard_clustering


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

    make_option(c("-c", "--celldex"),
                default='ImmGen', metavar='ImmGen', type="character",
                help="Choose from: ['ENCODE', 'HPCA', 'DICE', 'ImmGen', 'Monaco', 'MouseRNAseq', 'Hemato']"),

    make_option(c("-e", "--ensembl"),  default=FALSE,
                metavar="FALSE", action="store_true", type="logical",
                help="When querying celldex"),

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
figures_dir <- multiple_replacement(
    opt[['input-dir']], c('data'='figures', 'input'='output')
)

# Start Log
start_time <- Sys.time()
log <- log_open(paste0("cluster_groups-",
                       strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Main

config <- fromJSON(file=file.path(wd, opt[['input-dir']], opt[['config']]))

for (group_name in names(config)) {

    # ----------------------------------------------------------------------
    # Read Seurat Objects

    log_print(paste(Sys.time(), 'Reading Seurat objects...'))

    sample_names <- config[[group_name]]
    expr_mtxs <- new.env()

    for (sample_name in sample_names) {

        log_print(paste(Sys.time(), 'Reading...', sample_name))

        tryCatch({            
            expr_mtx <- read_10x(file.path(wd, opt[['input-dir']], sample_name))  # read seurat

            tmp_seurat_obj <- CreateSeuratObject(counts = expr_mtx, min.cells = 3) %>% 
                NormalizeData(verbose = FALSE) %>% 
                FindVariableFeatures(verbose = FALSE)

            tmp_seurat_obj$treatment <- sample_name
            tmp_seurat_obj <- subset(tmp_seurat_obj, subset = (
                (nCount_RNA < 20000) & (nCount_RNA > 1000) & (nFeature_RNA > 1000))
            )  # filter

            # seurat_obj <- run_standard_clustering(seurat_obj, ndim=40)  # cluster
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
    # Integrate data, then cluster

    log_print(paste(Sys.time(), 'Integrating data...'))

    integration_features <- SelectIntegrationFeatures(object.list = expr_mtxs)
    integration_anchors <- FindIntegrationAnchors(
        object.list = expr_mtxs,
        anchor.features = integration_features
    )
    seurat_obj <- IntegrateData(anchorset = integration_anchors)
    seurat_obj <- run_standard_clustering(seurat_obj, ndim=30)

    # save this for plotting later
    if (!troubleshooting) {
        if (!dir.exists(file.path(wd, opt[['output-dir']], 'integrated'))) {
            dir.create(file.path(wd, opt[['output-dir']], 'integrated'), recursive=TRUE)
        }
        save(seurat_obj,
             file=file.path(
                wd, opt[['output-dir']], 'integrated', paste0(group_name,'.RData')
            ))
    }

    # ----------------------------------------------------------------------
    # Label clusters

    log_print(paste(Sys.time(), 'Labeling clusters...'))

    # note: this overwrites the labels generated by the standard clustering workflow
    sce_counts <- as.SingleCellExperiment(seurat_obj)
    get_data=celldex_switch[[opt$celldex]]
    ref_data <- get_data(ensembl = opt[['ensembl']])  # default: ensembl=FALSE
    predictions <- SingleR(
        test=sce_counts,
        assay.type.test=1,
        ref=ref_data,
        labels=ref_data[['label.main']]
    )
    sce_counts[["cell_type"]] <- predictions[['labels']]


    # ----------------------------------------------------------------------
    # Plot UMAP

    log_print(paste(Sys.time(), 'Plotting...'))

    # not needed yet
    # plotScoreHeatmap(predictions)

    seurat_counts <- as.Seurat(sce_counts, counts = NULL)
    fig <- DimPlot(
        seurat_counts, reduction = "UMAP", 
        group.by = "cell_type",
        # split.by = "treatment",  # facet if necessary
        label = TRUE
    ) + ggtitle(group_name)

    if (!troubleshooting) {
        if (!dir.exists(file.path(wd, figures_dir, 'integrated'))) {
            dir.create(file.path(wd, figures_dir, 'integrated'), recursive=TRUE)
        }
        ggsave(file.path(wd, figures_dir, 'integrated',
                         paste0('umap-integrated-', group_name, '.png')),
               height=800, width=1000, dpi=300, units="px", scaling=0.5)
    }

    # Plot UMAP for specific gene of interest
    FeaturePlot(seurat_obj, 
                reduction = "umap", 
                features = opt[['gene-of-interest']],
                pt.size = 0.4, 
                order = TRUE,
                # split.by = "treatment",  # no need to facet
                min.cutoff = 'q10',
                label = FALSE) + ggtitle(paste(opt[['gene-of-interest']], 'in', group_name))
    if (!troubleshooting) {
        ggsave(
            file.path(
                wd, figures_dir, 'integrated',
                paste0('umap-integrated-', group_name, '-', tolower(opt[['gene-of-interest']]), '.png')
            ),
            height=800, width=800,
            dpi=300, units="px", scaling=0.5
        )
    }

}

end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
