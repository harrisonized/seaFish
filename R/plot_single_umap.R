## Bare minimum to plot single umap
## To be converted to shiny app


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


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-file"),
                default='data/scrnaseq-ballesteros/integrated/spleen/integrated_seurat.RData',
                metavar='data/scrnaseq-ballesteros/integrated/spleen/integrated_seurat.RData',
                type="character",
                help="output of scrnaseq_integrated_analysis.R"),

    make_option(c("-o", "--output-dir"),
                default="data/scrnaseq-ballesteros/integrated/spleen",
                metavar="data/scrnaseq-ballesteros/integrated/spleen", type="character",
                help="set the output directory for the data"),

    make_option(c("-f", "--figures-dir"),
                default="figures/scrnaseq-ballesteros/integrated/spleen",
                metavar="figures/scrnaseq-ballesteros/integrated/spleen", type="character",
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
log <- log_open(paste0("plot_single_umap-",
                       strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Read Data

log_print(paste(Sys.time(), 'Reading data...'))

integrated_seurat <- load_rdata(file.path(wd, opt[['input-file']]))

log_print(paste(Sys.time(), 'Plotting...'))

# ----------------------------------------------------------------------
# Plot Data

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
    ggsave(file.path(wd, opt[['figures-dir']],
                     paste0('umap-', dataset, '-integrated-', tolower(opt[['gene-of-interest']]), '.png')),
           height=750, width=1200, dpi=300, units="px", scaling=0.5)
}


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
