## Code used to analyze scRNAseq data
## Based on: https://satijalab.org/seurat/articles/integration_introduction.htmls
##
## This code is adapted from DIY_scRNAseq.R, which is freely available here: 
## https://diytranscriptomics.com/scripts
## See: https://github.com/harrisonized/diy-transcriptomics

wd = dirname(this.path::here())  # wd = '~/github/R/harrisonRTools'
suppressMessages(library('Seurat'))
library('Matrix')
# suppressMessages(library('DropletUtils'))
suppressMessages(library('dplyr'))
library('readr')
library('ggplot2')
library('optparse')
library('logr')
source(file.path(wd, 'R', 'functions', 'file_io.R'))  # read_10x
source(file.path(wd, 'R', 'functions', 'bio.R'))  # run_standard_clustering


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--inputs"),
                default=paste('data/scrnaseq-ballesteros/10x_counts/spleen1/counts',
                              'data/scrnaseq-ballesteros/10x_counts/spleen2/counts', sep=':'),
                metavar=paste('data/scrnaseq-ballesteros/10x_counts/spleen1/counts',
                              'data/scrnaseq-ballesteros/10x_counts/spleen2/counts', sep=':'),
                type="character",
                help="Colon separated list: path/to/dir1:path/to/dir2:path/to/dir3"),

    make_option(c("-o", "--output-dir"),
                default="figures/scrnaseq-ballesteros/integrated/spleen",
                metavar="figures/scrnaseq-ballesteros/integrated/spleen", type="character",
                help="set the output directory for the figures"),

    make_option(c("-l", "--labels"),  default='spleen1:spleen2', metavar='spleen1:spleen2',
                type="character", help="Colon separated list of labels"),

    make_option(c("-g", "--gene-of-interest"), default="Dnase1l1",
                metavar="Dnase1l1", type="character",
                help="choose a gene"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
troubleshooting = opt[['troubleshooting']]

# Start Log
start_time = Sys.time()
log <- log_open(paste0("scrnaseq_eda-",
                       strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Create Seurat objects, then cluster

log_print(paste(Sys.time(), 'Reading Seurat objects...'))

inputs = unlist(strsplit(opt[['inputs']], ':'))
labels = unlist(strsplit(opt[['labels']], ':'))
names(labels) = inputs

expr_mtxs = c()
for (input in inputs) {
    label = labels[[input]]

    log_print(paste(Sys.time(), 'Reading...', label))

    # Perform single umap
    expr_mtx <- read_10x(file.path(wd, input))  # read seurat
    seurat_obj <- CreateSeuratObject(counts = expr_mtx, min.cells = 3) %>% 
        NormalizeData(verbose = FALSE) %>% 
        FindVariableFeatures(verbose = FALSE)
    seurat_obj$treatment <- label
    seurat_obj <- subset(seurat_obj, subset = (
        (nCount_RNA < 20000) & (nCount_RNA > 1000) & (nFeature_RNA > 1000))
    )  # filter

    seurat_obj <- run_standard_clustering(seurat_obj, ndim=40)  # cluster

    # Plot Unintegrated
    DimPlot(seurat_obj, reduction = "umap", split.by = "orig.ident", label = TRUE)
    if (!troubleshooting) {
        ggsave(file.path(wd, opt[['output-dir']], paste0('umap-single-', label, '.png')),
               height=750, width=1200, dpi=300, units="px", scaling=0.5)
    }

    expr_mtxs <- append(expr_mtxs, seurat_obj)
}


# ----------------------------------------------------------------------
# Integrate data, then cluster

log_print(paste(Sys.time(), 'Integrating data...'))

integration_features <- SelectIntegrationFeatures(object.list = expr_mtxs)
integration_anchors <- FindIntegrationAnchors(
    object.list = expr_mtxs,
    anchor.features = integration_features
)
integrated_seurat <- IntegrateData(anchorset = integration_anchors)
integrated_seurat <- run_standard_clustering(integrated_seurat, ndim=30)

log_print(paste(Sys.time(), 'Plotting...'))

# Plot UMAP with clusters highlighted
DimPlot(integrated_seurat, reduction = "umap", label = TRUE)
if (!troubleshooting) {
    ggsave(file.path(wd, opt[['output-dir']], 'umap-integrated.png'),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}

# Plot UMAP for specific gene of interest
FeaturePlot(integrated_seurat, 
            reduction = "umap", 
            features = opt[['gene-of-interest']],
            pt.size = 0.4, 
            order = TRUE,
            split.by = "treatment",
            min.cutoff = 'q10',
            label = FALSE)
if (!troubleshooting) {
    ggsave(file.path(wd, opt[['output-dir']], paste0('umap-', tolower(opt[['gene-of-interest']]), '.png')),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
