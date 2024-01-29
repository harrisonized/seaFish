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
suppressPackageStartupMessages(library('plotly'))
library("ggplot2")
library('rjson')
library('optparse')
library('logr')
import::from(magrittr, '%>%')
import::from(dplyr, 'group_by', 'top_n', 'ungroup', 'slice_head')
import::from(DropletUtils, 'barcodeRanks')
import::from(SingleR, 'SingleR', 'plotScoreHeatmap')
import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    'read_10x', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'multiple_replacement', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'plotting.R'),
    'savefig', 'plot_violin', 'plot_waterfall', .character_only=TRUE)
import::from(file.path(wd, 'R', 'functions', 'single_cell.R'),
    'celldex_switch', .character_only=TRUE)


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
log <- log_open(paste0("integrate_data-",
                       strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Subroutines


#' Draw QC plots
#' 
#' @description
#' 
draw_qc_plots <- function(
    seurat_obj,
    dirpath,
    sample_name = 'sample',
    threshold_data = NULL,
    troubleshooting=FALSE
) {

    # Figure 1. Standard QC Metrics
    plot_violin(seurat_obj, group.by='sample_name', threshold_data=threshold_data)
    savefig(file.path(dirpath, paste0('violin-qc-', sample_name, '.png')),
            makedir=TRUE, troubleshooting=troubleshooting)

    # Figure 2. Waterfall Plot of Read Counts
    # See: https://sydneybiox.github.io/SingleCellPlus/qc.html#3_qc1:_waterfall_plot
    barcode_ranks <- barcodeRanks(seurat_obj)
    plot_waterfall(barcode_ranks)
    savefig(file.path(dirpath, paste0('waterfall-counts_per_cell-', sample_name, '.png')),
            troubleshooting=troubleshooting)


    # Figure 3. Waterfall Plot of Gene Representation
    # See: https://ucdavis-bioinformatics-training.github.io/2017_2018-single-cell-RNA-sequencing-Workshop-UCD_UCB_UCSF/day2/scRNA_Workshop-PART2.html
    cells_per_gene <- data.frame(
        num_cells=sort(
            rowSums(seurat_obj[['RNA']]@counts>=2),
            decreasing=TRUE)
    )
    cells_per_gene['rank'] = 1:nrow(seurat_obj[['RNA']]@counts)

    ggplot(cells_per_gene, aes(x = rank, y=num_cells)) +
        geom_point(size = 0.2) +
        labs(title = "Cells Per Gene ( >=2 )",
             x = "Gene Rank",
             y = "Number of Cells")
    savefig(file.path(dirpath, paste0('waterfall-cells_per_gene-', sample_name, '.png')),
            troubleshooting=troubleshooting)


    # Figure 4. Genes vs. Read Counts per cell
    ggplot(seurat_obj@meta.data,
           aes(.data[['nCount_RNA']], .data[['nFeature_RNA']])) +
           geom_point(alpha = 0.7, size = 0.5) +
           labs(title = "Read Counts vs Sequencing Depth",
                x = "Read Counts",
                y = "Number of Genes")
    savefig(file.path(dirpath, paste0('violin-reads_vs_seqdepth-', sample_name, '.png')),
            troubleshooting=troubleshooting)

}


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

            expr_mtx <- read_10x(file.path(wd, opt[['input-dir']], sample_name))
            tmp_seurat_obj <- CreateSeuratObject(counts = expr_mtx, min.cells = 3)
            tmp_seurat_obj$sample_name <- sample_name  # consider using orig.ident
            tmp_seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
                tmp_seurat_obj, pattern = "^mt-"
            )

            upper_rna_count <- median(tmp_seurat_obj$nCount_RNA) + 3*sd(tmp_seurat_obj$nCount_RNA)
            upper_rna_feature <- median(tmp_seurat_obj$nFeature_RNA) + 3*sd(tmp_seurat_obj$nFeature_RNA)
            upper_pct_mt <- min(median(tmp_seurat_obj$percent.mt) + 3*sd(tmp_seurat_obj$percent.mt), 5)
            
            threshold_data <- data.frame(
                sample_name=sample_name,
                col=c('nCount_RNA', 'nCount_RNA', 'nFeature_RNA', 'nFeature_RNA', 'percent.mt'),
                threshold=c(1000, upper_rna_count, 300, upper_rna_feature, upper_pct_mt)
            )

            draw_qc_plots(
                tmp_seurat_obj,
                dirpath=file.path(wd, figures_dir, sample_name, 'qc'),
                sample_name=sample_name,
                threshold_data = threshold_data,
                troubleshooting=troubleshooting
            )

            # ----------------------------------------------------------------------
            # Filter

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


            # ----------------------------------------------------------------------
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

    log_print(paste(Sys.time(), 'Integrating data...'))

    features_list <- SelectIntegrationFeatures(object.list = expr_mtxs, nfeatures = 2000)
    anchors <- FindIntegrationAnchors(object.list = expr_mtxs, anchor.features = features_list)
    seurat_obj <- IntegrateData(anchorset = anchors)
    DefaultAssay(seurat_obj) <- "integrated"

    # Run the standard workflow
    seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
    seurat_obj <- RunPCA(seurat_obj, npcs = ndim, verbose = FALSE)  # required for RunUMAP
    seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:ndim)
    seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:ndim)
    seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
    

    # ----------------------------------------------------------------------
    # Label Clusters
    
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

    # save object
    if (!troubleshooting) {
        if (!dir.exists(file.path(wd, opt[['output-dir']], 'integrated'))) {
            dir.create(file.path(wd, opt[['output-dir']], 'integrated'), recursive=TRUE)}
        save(seurat_obj,
             file=file.path(wd, opt[['output-dir']], 'integrated', paste0(group_name,'.RData')) )
    }


    # ----------------------------------------------------------------------
    # Plot

    # Plot unlabeled clusters
    DimPlot(seurat_obj,
            reduction = "umap",
            group.by = "sample_name",
            # split.by = "sample_name",
            label = TRUE) +
        theme(strip.text.x = element_blank(), plot.title = element_text(hjust = 0.5)) +
        ggtitle("Unlabeled Clusters")
    savefig(file.path(wd, figures_dir, 'integrated', group_name,
                      paste0('umap-integrated-', group_name, '-', 'unlabeled.png')),
            width=1000, makedir=TRUE, troubleshooting=troubleshooting)

    # Plot CellDex labeled clusters
    DimPlot(seurat_obj,
            reduction = "umap", 
            group.by = "cell_type",
            # split.by = "sample_name",  # facet if necessary
            label = TRUE
        ) + ggtitle(group_name)
    savefig(file.path(wd, figures_dir, 'integrated', group_name,
                      paste0('umap-integrated-', group_name, '-', tolower(opt[['celldex']]), '_labeled.png')),
            width=1000, troubleshooting=troubleshooting)

    # Plot CellDex QC Heatmap
    fig <- plotScoreHeatmap(predictions)
    savefig(file.path(wd, figures_dir, 'integrated', group_name,
                      paste0('heatmap-', group_name, '-', 'predictions.png')),
            fig=fig, height=1600, width=2400, dpi=300, troubleshooting=troubleshooting)


    # ----------------------------------------------------------------------
    # Remove rare populations

    # Number of cells per label
    num_cells_per_label <- table(predictions['labels'])
    num_cells_per_label <- num_cells_per_label[sort.list(num_cells_per_label, decreasing=TRUE)]
    

    # ----------------------------------------------------------------------
    # Plot

    fig <- ggplot( data=data.frame(num_cells_per_label), aes(x=labels, y=Freq)) +
        geom_bar(stat="identity", fill="steelblue") +
        labs(title = "Number of Cells Per Label",
                     x = "Label",
                     y = "Number of Cells") +
        scale_x_discrete(guide = guide_axis(angle = 60)) +
        theme_minimal()
    savefig(file.path(wd, figures_dir, 'integrated', group_name,
                      paste0('bar-', group_name, '-', 'num_cells_per_label.png')),
            height=800, width=1000, troubleshooting=troubleshooting)


    pct_cells_per_label <- num_cells_per_label/sum(num_cells_per_label)
    populations_to_keep <- names(num_cells_per_label[(pct_cells_per_label > 0.03)])
    seurat_obj_filt <- seurat_obj[, seurat_obj$cell_type %in% populations_to_keep]

    # Plot CellDex labeled clusters
    DimPlot(seurat_obj,
            reduction = "umap", 
            group.by = "cell_type",
            split.by = "sample_name",  # facet if necessary
            label = TRUE
        ) + ggtitle(group_name)
    savefig(file.path(wd, figures_dir, 'integrated', group_name,
                      paste0('umap-integrated-', group_name, '-', tolower(opt[['celldex']]), '_labeled-filtered.png')),
            width=1000, troubleshooting=troubleshooting)

    # Genes and clusters Heatmap
    Idents(seurat_obj) <- seurat_obj$cell_type

    # see: https://satijalab.org/seurat/articles/pbmc3k_tutorial#finding-differentially-expressed-features-cluster-biomarkers
    markers <- FindAllMarkers(
        seurat_obj,
        test.use = "wilcox",
        logfc.threshold = 0.25,  # speed
        min.pct = 0.1,
        only.pos = TRUE
    )
    top_markers <- markers %>%
        group_by(cluster) %>%
        dplyr::filter(avg_log2FC > 1) %>%
        slice_head(n = 12) %>%
        ungroup()

    # See: https://github.com/satijalab/seurat/issues/2960
    seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
    DoHeatmap(seurat_obj, features = top_markers$gene)
    savefig(file.path(wd, figures_dir, 'integrated', group_name,
                         paste0('heatmap-', group_name, '-', 'clusters.png')),
            width=1000, troubleshooting=troubleshooting)


    log_print(paste("Loop completed in:", difftime(Sys.time(), loop_start_time)))
}

end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
