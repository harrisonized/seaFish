## Code used to analyze scRNAseq data
## Cluster each sample individually
## Based on: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

wd = dirname(this.path::here())  # wd = '~/github/R/seaFish'
suppressPackageStartupMessages(library('Seurat'))
library('ggplot2')
library('optparse')
library('logr')
import::from('magrittr', '%>%', .character_only=TRUE)
import::from('SingleR', 'SingleR', .character_only=TRUE)
# import::from('DropletUtils', 'emptyDrops', 'write10xCounts', .character_only=TRUE)
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
log <- log_open(paste0("cluster_samples-",
                       strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Main

samples <- list.dirs(file.path(wd, opt[['input-dir']]), full.names=FALSE, recursive=FALSE)
log_print(paste(
    Sys.time(), 'Found', length(samples), 'samples...', paste(samples, collapse=', ')
))

for (sample_name in samples) {

    tryCatch({

        # ----------------------------------------------------------------------
        # Read Data

        log_print(paste(Sys.time(), 'Reading', sample_name, '...'))

        expr_mtx <- read_10x(file.path(wd, opt[['input-dir']], sample_name))
        seurat_obj <- CreateSeuratObject(counts = expr_mtx, min.cells = 3) %>% 
            NormalizeData(verbose = FALSE) %>% 
            FindVariableFeatures(verbose = FALSE)
        seurat_obj$sample_name <- sample_name

        # QC plot 1
        VlnPlot(seurat_obj, c("nCount_RNA", "nFeature_RNA"), pt.size = 0.1)
        if (!troubleshooting) {
            ggsave(file.path(wd, figures_dir, sample_name, 'qc',
                       paste0('violin-', sample_name, '-raw_qc.png')),
                   height=800, width=1200, dpi=300, units="px", scaling=0.5)
        }

        # Note: emptyDrops cannot be run twice
        # For the ballesteros-2020 dataset, this has already been run
        # drop_stats <- DropletUtils::emptyDrops(expr_mtx)

        # filter data
        seurat_obj <- subset(seurat_obj,
            subset = ((nCount_RNA < 20000) & (nCount_RNA > 1000) & (nFeature_RNA > 1000))
        )

        # QC plot 2
        ggplot(seurat_obj@meta.data, aes(nCount_RNA, nFeature_RNA)) +
            geom_point(alpha = 0.7, size = 0.5) +
            labs(x = "Total UMI counts per cell", y = "Number of genes detected")
        if (!troubleshooting) {
            ggsave(file.path(wd, figures_dir, sample_name, 'qc',
                       paste0('scatter-', sample_name, '-num_genes_vs_counts.png')),
                   height=800, width=1200, dpi=300, units="px", scaling=0.5)
        }


        # ----------------------------------------------------------------------
        # Cluster

        log_print(paste(Sys.time(), 'Running standard clustering workflow...'))

        seurat_obj <- run_standard_clustering(seurat_obj, ndim=30)

        # Plot UMAP
        DimPlot(seurat_obj,
                reduction = "umap",
                split.by = "orig.ident",
                label = TRUE) +
            theme(strip.text.x = element_blank(), plot.title = element_text(hjust = 0.5)) +
            ggtitle("Unlabeled Clusters")
        if (!troubleshooting) {
            ggsave(file.path(wd, figures_dir, sample_name,
                             paste0('umap-', sample_name, '-clusters.png')),
                   height=800, width=800, dpi=300, units="px", scaling=0.5)
        }

        # Plot UMAP for specific gene of interest
        gene = opt[['gene-of-interest']]
        FeaturePlot(seurat_obj,
            reduction = "umap",
            # split.by = "orig.ident",
            features = c(gene),
            min.cutoff = 'q10',
            pt.size = 0.4, order = TRUE, label = FALSE
        )
        if (!troubleshooting) {
            ggsave(file.path(wd, figures_dir, sample_name,
                       paste0('umap-', sample_name, '-', tolower(gene), '.png')),
                   height=800, width=800, dpi=300, units="px", scaling=0.5)
        }

        # ----------------------------------------------------------------------
        # Label clusters

        log_print(paste(Sys.time(), 'Labeling clusters...'))

        # Note: this overwrites the labels generated by the standard clustering workflow
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

        # plotScoreHeatmap(predictions)

        # Plot UMAP
        seurat_obj <- as.Seurat(sce_counts, counts = NULL)
        DimPlot(seurat_obj,
                reduction = "UMAP", 
                group.by = "cell_type",
                # split.by = "sample_name",  # facet if necessary
                label = TRUE
            ) + ggtitle("Labeled Clusters")
        if (!troubleshooting) {
            ggsave(file.path(wd, figures_dir, sample_name,
                       paste0('umap-', sample_name, '-', tolower(opt[['celldex']]), '_labeled.png')),
                   height=800, width=1000, dpi=300, units="px", scaling=0.5)
        }
    },

    # pass
    error = function(condition) {
        log_print(paste("WARNING: SAMPLE", sample_name, "NOT PROCESSED!!!"))
        log_print(paste("Error message: ", conditionMessage(condition)))
        NULL
    })

}

end_time <- Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
