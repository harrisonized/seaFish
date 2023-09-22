## Code used to analyze scRNAseq data

wd = dirname(this.path::here())  # wd = '~/github/R/harrisonRTools'
suppressMessages(library('Seurat'))
library('Matrix')
# suppressMessages(library('DropletUtils'))
library('tibble')
suppressMessages(library('dplyr'))
library('rjson')
library('R2HTML')
library('readr')
suppressMessages(library('textshape'))
library('scales')
library('ggplot2')
suppressMessages(library('DT'))
library('optparse')
library('logr')
source(file.path(wd, 'R', 'functions', 'file_io.R'))  # read_10x


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-dir"),
                default='data/ballesteros-scrnaseq/10x_counts/spleen1/raw',
                metavar='data/ballesteros-scrnaseq/10x_counts/spleen1/raw',
                type="character",
                help="directory containing standard 10x output: barcodes.tsv, genes.tsv, and matrix.mtx"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
troubleshooting = opt['troubleshooting'][[1]]

# Start Log
start_time = Sys.time()
log <- log_open(paste0("scrnaseq_filter_drops",
                       strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Read Data

log_print(paste(Sys.time(), 'Reading data...'))

expr_mtx <- read_10x(
    file.path(wd, opt['input-dir'][[1]]),
    matrix_file='GSM4239566_Spleen1_mm10_matrix.mtx.gz',
    genes_file='GSM4239566_Spleen1_mm10_genes.tsv.gz',
    barcodes_file='GSM4239566_barcodes_Spleen1.tsv.gz'
)



end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
