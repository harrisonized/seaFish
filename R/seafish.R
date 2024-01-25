## Code used to analyze scRNAseq data
## Integrates the datasets based on the information in config.json
## Based on: https://satijalab.org/seurat/archive/v4.3/integration_introduction

wd = dirname(this.path::here())  # wd = '~/github/R/seaFish'
suppressPackageStartupMessages(library('Seurat'))
library("ggplot2")
library('rjson')
library('optparse')
library('logr')
import::from(magrittr, '%>%')
import::from(grid, 'grid.newpage', 'grid.draw')
import::from(SingleR, 'SingleR', 'plotScoreHeatmap')
# import::from(DropletUtils, 'emptyDrops', 'write10xCounts')
import::from(ggplot2,
    'theme', 'element_text', 'element_blank', 'ggtitle', 'ggsave'
)
import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    'read_10x', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'multiple_replacement', .character_only=TRUE)
import::from(file.path(wd, 'R', 'functions', 'single_cell.R'),
    'celldex_switch', 'run_standard_analysis_workflow', .character_only=TRUE)


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
ndim=opt[['ndim']]

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
    # Read Data

    loop_start_time <- Sys.time()
    log_print(paste(loop_start_time, 'Reading Seurat objects...'))

    sample_names <- config[[group_name]]
    
    expr_mtxs <- new.env()
    for (sample_name in sample_names) {

        log_print(paste(Sys.time(), 'Reading...', sample_name))

        tryCatch({

            # ----------------------------------------------------------------------
            # Preprocessing

            # See: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
            # According to this tutorial, the steps are run in this order:
            # 1. Create QC violin plot
            # 2. Filter first
            # 3. Normalize after filtering

            expr_mtx <- read_10x(file.path(wd, opt[['input-dir']], sample_name))

            tmp_seurat_obj <- CreateSeuratObject(counts = expr_mtx, min.cells = 3)
            tmp_seurat_obj$sample_name <- sample_name  # consider using orig.ident


            # ----------------------------------------------------------------------
            # Quality Control Filters

            tmp_seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
                tmp_seurat_obj, pattern = "^mt-"
            )

            # QC plot
            VlnPlot(tmp_seurat_obj,
                   c("nCount_RNA", "nFeature_RNA", "percent.mt"),
                   ncol = 3, pt.size = 0.1)
            if (!troubleshooting) {
                if (!dir.exists(file.path(wd, figures_dir, sample_name, 'qc'))) {
                    dir.create(file.path(wd, figures_dir, sample_name, 'qc'), recursive=TRUE)
                }
                ggsave(
                    file.path(wd, figures_dir, sample_name, 'qc',
                              paste0('violin-', sample_name, '-qc.png')),
                    height=800, width=1200, dpi=300, units="px", scaling=0.5
                )
            }

            # QC plot 2
            ggplot(tmp_seurat_obj@meta.data,
                   aes(.data[['nCount_RNA']], .data[['nFeature_RNA']])) +
                   geom_point(alpha = 0.7, size = 0.5) +
                   labs(x = "Counts Per Cell", y = "Genes Per Cell")
            if (!troubleshooting) {
                ggsave(file.path(wd, figures_dir, sample_name, 'qc',
                           paste0('scatter-', sample_name, '-num_genes_vs_counts.png')),
                       height=800, width=1200, dpi=300, units="px", scaling=0.5)
            }

            ## TODO:
            ## 1. DropletUtils::barcodeRanks(tmp_seurat_obj)
            ## 2. filter below knee point
            ## see: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8848315/

            #' Filters
            #' See: https://www.10xgenomics.com/analysis-guides/common-considerations-for-quality-control-filters-for-single-cell-rna-seq-data
            #' Lower bound to filter low quality cells
            #' Upper bound to filter potential doublets (use median+3*stdev)
            upper_rna_count <- median(tmp_seurat_obj$nCount_RNA) + 3*sd(tmp_seurat_obj$nCount_RNA)
            upper_rna_feature <- median(tmp_seurat_obj$nFeature_RNA) + 3*sd(tmp_seurat_obj$nFeature_RNA)
            upper_pct_mt <- median(tmp_seurat_obj$percent.mt) + 3*sd(tmp_seurat_obj$percent.mt)

            # Filter
            tmp_seurat_obj <- subset(
                tmp_seurat_obj,
                subset = (
                    (nCount_RNA > 1000) & (nCount_RNA <= upper_rna_count) &
                    (nFeature_RNA > 1000) & (nFeature_RNA <= upper_rna_feature) &
                    (percent.mt <= upper_pct_mt)
                )
            )

            # Normalize after filtering
            tmp_seurat_obj <- tmp_seurat_obj %>%
                NormalizeData(verbose = FALSE) %>% 
                FindVariableFeatures(
                    verbose = FALSE,
                    selection.method = "vst",
                    nfeatures = 2000
                )  # same function used in SelectIntegrationFeatures


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

    log_print(paste(Sys.time(), 'Integrating data...'))

    features <- c(SelectIntegrationFeatures(
        object.list = expr_mtxs,
        nfeatures = 2000  # these are applied onto the "integrated" assay
    ))  # add any genes here, including the gene-of-interest does not really make a difference
    anchors <- FindIntegrationAnchors(
        object.list = expr_mtxs,
        anchor.features = features
    )
    seurat_obj <- IntegrateData(anchorset = anchors)


    # ----------------------------------------------------------------------
    # Cluster assignment

    DefaultAssay(seurat_obj) <- "integrated"
    
    # Run the standard workflow for visualization and clustering
    seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
    seurat_obj <- RunPCA(seurat_obj, npcs = ndim, verbose = FALSE)  # required for RunUMAP to work
    seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:ndim)
    seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:ndim)
    seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)


    log_print(paste(Sys.time(), 'Labeling clusters...'))

    ref_data <- celldex_switch[[opt$celldex]](ensembl=opt[['ensembl']])
    predictions <- SingleR(
        test=as.SingleCellExperiment(seurat_obj),
        assay.type.test=1,
        ref=ref_data,
        labels=ref_data[['label.main']]
    )
    seurat_obj$cell_type <- predictions[['labels']]

    # switch back after cluster assignment
    DefaultAssay(seurat_obj) <- "RNA"
    Idents(seurat_obj) <- seurat_obj@meta.data$sample_name

    # save
    if (!troubleshooting) {
        if (!dir.exists(file.path(wd, opt[['output-dir']], 'integrated'))) {
            dir.create(file.path(wd, opt[['output-dir']], 'integrated'), recursive=TRUE)
        }
        save(seurat_obj,
             file=file.path(wd, opt[['output-dir']], 'integrated',
                            paste0(group_name,'.RData')) )
    }

    # ----------------------------------------------------------------------
    # Visualization

    log_print(paste(Sys.time(), 'Plotting...'))

    # Plot unlabeled clusters
    DimPlot(seurat_obj,
            reduction = "umap",
            split.by = "sample_name",
            label = TRUE) +
        theme(strip.text.x = element_blank(), plot.title = element_text(hjust = 0.5)) +
        ggtitle("Unlabeled Clusters")
    if (!troubleshooting) {
        if (!dir.exists(file.path(wd, figures_dir, 'integrated', group_name))) {
            dir.create(file.path(wd, figures_dir, 'integrated', group_name), recursive=TRUE)
        }
        ggsave(file.path(wd, figures_dir, 'integrated', group_name,
                         paste0('umap-integrated-', group_name, '-', 'unlabeled.png')),
               height=800, width=1000, dpi=300, units="px", scaling=0.5)
    }

    # Plot CellDex labeled clusters
    DimPlot(seurat_obj,
            reduction = "umap", 
            group.by = "cell_type",
            # split.by = "sample_name",  # facet if necessary
            label = TRUE
        ) + ggtitle(group_name)
    if (!troubleshooting) {
        ggsave(file.path(wd, figures_dir, 'integrated', group_name,
                         paste0('umap-integrated-', group_name, '-', tolower(opt[['celldex']]), '_labeled.png')),
               height=800, width=1000, dpi=300, units="px", scaling=0.5)
    }

    # Plot CellDex QC Heatmap
    fig <- plotScoreHeatmap(predictions)
    if (!troubleshooting) {
        png(file.path(wd, figures_dir, 'integrated', group_name,
                      paste0('heatmap-', group_name, '-', 'predictions.png')),
            width=2400, height=1600, res=360, units='px'
        )
        grid::grid.newpage()
        grid::grid.draw(fig$gtable)
        dev.off()
    }

    # ----------------------------------------------------------------------
    # Gene of Interest

    gene = opt[['gene-of-interest']]
    
    # Plot UMAP
    FeaturePlot(seurat_obj, 
                reduction = "umap", 
                features = gene,
                pt.size = 0.4, 
                order = TRUE,
                # split.by = "sample_name",  # no need to facet
                min.cutoff = 'q10',
                label = FALSE) + ggtitle(paste(gene, 'in', group_name))
    if (!troubleshooting) {
        ggsave(file.path(wd, figures_dir, 'integrated', group_name,
                   paste0('umap-integrated-', group_name, '-', tolower(gene), '.png')),
               height=800, width=800, dpi=300, units="px", scaling=0.5)
    }

    # Plot Ridge
    RidgePlot(
        seurat_obj,
        features=gene,
        group.by = "cell_type"
    ) + ggtitle(paste(gene, 'in', group_name))
    if (!troubleshooting) {
        ggsave(file.path(wd, figures_dir, 'integrated', group_name,
                   paste0('ridgeline-integrated-', group_name, '-', tolower(gene), '.png')),
               height=800, width=1200, dpi=300, units="px", scaling=0.5)
    }

    log_print(paste("Loop completed in:", difftime(Sys.time(), loop_start_time)))
}

end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
