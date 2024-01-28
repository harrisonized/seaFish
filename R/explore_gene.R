## Code used to analyze scRNAseq data
## Integrates the datasets based on the information in config.json
## Based on: https://satijalab.org/seurat/archive/v4.3/integration_introduction

wd = dirname(this.path::here())  # wd = '~/github/R/seaFish'
suppressPackageStartupMessages(library('Seurat'))
suppressPackageStartupMessages(library('plotly'))
library("ggplot2")
library('optparse')
library('logr')
import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    'load_rdata', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'multiple_replacement', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'plotting.R'),
    'savefig', .character_only=TRUE)


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-dir"),
                default='data/ballesteros-2020/output/integrated',
                metavar='data/ballesteros-2020/output/integrated',
                type="character",
                help="directory containing Seurat objects stored as .RData files"),
    
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
gene = opt[['gene-of-interest']]

# Start Log
start_time <- Sys.time()
log <- log_open(paste0("explore_gene-",
                       strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Main

filepaths <- list.files(file.path(wd, opt[['input-dir']]), full.names=TRUE)
log_print(paste(Sys.time(), 'Found', length(filepaths), 'files.'))

for (filepath in filepaths) {

    loop_start_time <- Sys.time()
    filename <- tools::file_path_sans_ext(basename(filepath))
    log_print(paste(Sys.time(), 'Processing', filename, '...'))

    seurat_obj <- load_rdata(filepath)

    # ----------------------------------------------------------------------
    # Gene of Interest Expression

    # Plot UMAP
    FeaturePlot(seurat_obj, 
                reduction = "umap", 
                features = gene,
                pt.size = 0.4, 
                order = TRUE,
                # split.by = "sample_name",  # no need to facet
                min.cutoff = 'q10',
                label = FALSE) + ggtitle(paste(gene, 'in', filename))

    savefig(file.path(wd, figures_dir, filename,
                      paste0('umap-integrated-', filename, '-', tolower(gene), '.png')),
            width=800, makedir=TRUE, troubleshooting=troubleshooting)

    # plot violin
    # note: ridge looks poor if there are not many positive cells
    VlnPlot(
        seurat_obj,
        features=gene,
        group.by = "cell_type"
    ) + ggtitle(paste(gene, 'in', filename))
    savefig(file.path(wd, figures_dir, filename,
                   paste0('violin-integrated-', filename, '-', tolower(gene), '.png')),
            height=800, width=1200, troubleshooting=troubleshooting)


    # TODO: bar chart of percent of each cell type with Dnase1l1 expression

    # ----------------------------------------------------------------------
    # Gene of Interest Correlations
    # TODO: check if this is correct

    # Volcano plot
    # See: https://github.com/satijalab/seurat/issues/4118
    Idents(seurat_obj, 
        WhichCells(object = seurat_obj,
                   expression = Dnase1l1 > 0,
                   slot = 'data')
    ) <- 'pos'
    
    # takes ~1 min
    markers <- FindMarkers(
        seurat_obj,
        ident.1 = 'pos',
        test.use = "wilcox",
        logfc.threshold = 0.01,
        min.pct = 0.25
    )
    markers['gene'] <- rownames(markers)

    # plot
    fig <- ggplot(markers) +
        aes(y=-log10(p_val_adj), x=avg_log2FC, text = paste("Symbol:", gene)) +
        geom_point(size=0.5) +
        geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", linewidth=1) +
        geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=1) +
        geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=1) +
        annotate("rect", xmin = 1, xmax = 2, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#BE684D") +
        annotate("rect", xmin = -1, xmax = -2, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#2C467A") +
        labs(title="Volcano plot") +
        theme_bw()
    ggplotly(fig)

    log_print(paste("Loop completed in:", difftime(Sys.time(), loop_start_time)))
}

end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
