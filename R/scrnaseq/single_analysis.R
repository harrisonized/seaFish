## Code used to analyze scRNAseq data
## Based on: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
##
## This code is adapted from DIY_scRNAseq.R, which is freely available here: 
## https://diytranscriptomics.com/scripts
## See: https://github.com/harrisonized/diy-transcriptomics


wd = dirname(dirname(this.path::here()))  # wd = '~/github/R/harrisonRTools'
suppressMessages(library('Seurat'))
library('Matrix')
# suppressMessages(library('DropletUtils'))
suppressMessages(library('dplyr'))
library('readr')
library('ggplot2')
library('optparse')
library('logr')
source(file.path(wd, 'R', 'functions', 'file_io.R'))  # read_10x


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-dir"),
                default='data/scrnaseq-ballesteros/10x_counts/spleen1/counts',
                metavar='data/scrnaseq-ballesteros/10x_counts/spleen1/counts',
                type="character",
                help="directory containing standard 10x output: barcodes.tsv, genes.tsv, and matrix.mtx"),
    
    make_option(c("-o", "--output-dir"),
                default="figures/scrnaseq-ballesteros/10x_counts/spleen1",
                metavar="figures/scrnaseq-ballesteros/10x_counts/spleen1", type="character",
                help="set the output directory for the figures"),

    make_option(c("-g", "--gene-of-interest"), default="Dnase1l1",
                metavar="Dnase1l1", type="character",
                help="choose a gene"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
dataset = basename(opt[['output-dir']])
troubleshooting = opt[['troubleshooting']]

# Start Log
start_time = Sys.time()
log <- log_open(paste0("scrnaseq_single_analysis-",
                       strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Read Data

log_print(paste(Sys.time(), 'Reading data...'))

expr_mtx <- read_10x(file.path(wd, opt[['input-dir']]))

# data already filtered
# this will throw an error
# drop_stats <- DropletUtils::emptyDrops(expr_mtx)

seurat_obj <- CreateSeuratObject(counts = expr_mtx, min.cells = 3) %>% 
    NormalizeData(verbose = FALSE) %>% 
    FindVariableFeatures(verbose = FALSE)

# QC plot
VlnPlot(seurat_obj, c("nCount_RNA", "nFeature_RNA"), pt.size = 0.1)
if (!troubleshooting) {
    ggsave(file.path(wd, opt[['output-dir']], paste0('violin-', dataset, '-raw_qc.png')),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# ----------------------------------------------------------------------
# Filter data

seurat_obj <- subset(
    seurat_obj,
    subset = ((nCount_RNA < 20000) & 
              (nCount_RNA > 1000) & 
              (nFeature_RNA > 1000))
)

# another QC plot
ggplot(seurat_obj@meta.data, aes(nCount_RNA, nFeature_RNA)) +
    geom_point(alpha = 0.7, size = 0.5) +
    labs(x = "Total UMI counts per cell", y = "Number of genes detected")

if (!troubleshooting) {
    ggsave(file.path(wd, opt[['output-dir']], paste0('scatter-', dataset, '-num_genes_vs_counts.png')),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


# ----------------------------------------------------------------------
# Standard Clustering Workflow

log_print(paste(Sys.time(), 'Running standard clustering workflow...'))

seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
seurat_obj <- RunPCA(seurat_obj, npcs = 40, verbose = FALSE)
seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:40)
seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:40)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)


# ----------------------------------------------------------------------
# Plot UMAPs

# Plot UMAP with clusters highlighted
DimPlot(seurat_obj, reduction = "umap", split.by = "orig.ident", label = TRUE)
if (!troubleshooting) {
    ggsave(file.path(wd, opt[['output-dir']], paste0('umap-', dataset, '-clusters.png')),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}

# Plot UMAP for specific gene of interest
gene = opt[['gene-of-interest']]  # gene = 'IGHM'
FeaturePlot(
    seurat_obj, reduction = "umap", features = c(gene),
    pt.size = 0.4,  min.cutoff = 'q10',
    # split.by = "orig.ident",
    order = TRUE, label = FALSE
)
if (!troubleshooting) {
    ggsave(file.path(wd, opt[['output-dir']], paste0('umap-', dataset, '-', tolower(gene), '.png')),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
