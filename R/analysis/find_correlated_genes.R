## Code used to analyze scRNAseq data
## Plots a volcano plot of genes related to the gene-of-interest

wd = dirname(dirname(this.path::here()))  # wd = '~/github/R/seaFish'
suppressPackageStartupMessages(library('Seurat'))
library("ggplot2")
library('optparse')
library('logr')
import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    'load_rdata', 'savefig', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'multiple_replacement', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'plotting.R'),
    'plot_volcano', .character_only=TRUE)


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-dir"), default='data/ballesteros-2020/output/rdata',
                metavar='data/ballesteros-2020/output/rdata', type="character",
                help="directory containing Seurat objects stored as .RData files"),
    
    make_option(c("-o", "--output-dir"),
                default='data/ballesteros-2020/output',
                metavar='data/ballesteros-2020/output',
                type="character",
                help="set the output directory"),

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
    opt[['input-dir']], c('output/rdata'='output', 'data'='figures')
)
gene = opt[['gene-of-interest']]

# Start Log
start_time <- Sys.time()
log <- log_open(paste0("find_correlated_genes-",
                       strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Main

filepaths <- list.files(file.path(wd, opt[['input-dir']]), full.names=TRUE)
log_print(paste(Sys.time(), 'Found', length(filepaths), 'files.'))

for (filepath in filepaths) {

    loop_start_time <- Sys.time()
    filename <- basename(filepath)
    log_print(paste(Sys.time(), 'Processing', filename, '...'))

    seurat_obj <- load_rdata(filepath)

    # ----------------------------------------------------------------------
    # Gene of Interest Correlations
    # See: https://github.com/satijalab/seurat/issues/4118

    cells <- WhichCells(object = seurat_obj, expression = Dnase1l1 > 0, slot = 'data')
    Idents(seurat_obj, cells) <- 'pos'
    
    # takes ~1 min
    markers <- FindMarkers(
        seurat_obj,
        ident.1 = 'pos',
        logfc.threshold = 0.01,  # if this is high, lower significance genes are not computed
        min.pct = 0.25,
        test.use = "wilcox",
        only.pos = FALSE
    )
    markers['gene'] <- rownames(markers)

    # Plot
    fig <- plot_volcano(markers)
    dirpath <- file.path(
        wd, figures_dir,  # figures/ballesteros-2020/output
        multiple_replacement(filename, c('-'='/', '.RData'='')),  # integrated/bm/
        'expression', tolower(gene), 'volcano'  # expression/dnase1l1/volcano
    )
    savefig(file.path(dirpath,
        paste0('volcano-', tools::file_path_sans_ext(filename), '-', tolower(gene), '.png')),
        troubleshooting=troubleshooting
    )

    log_print(paste("Loop completed in:", difftime(Sys.time(), loop_start_time)))
}

end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
