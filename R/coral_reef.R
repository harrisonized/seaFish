## Get an overview of all the datasets in the form of a dotplot

wd = dirname(this.path::here())  # wd = '~/github/R/seaFish'
# suppressPackageStartupMessages(library('ggplot2'))
library('optparse')
library('logr')

import::from(stringr, 'str_match')
import::from(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'multiple_replacement', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    'append_many_csv', 'savefig', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'plotting.R'),
    'plot_dotplot', .character_only=TRUE)

# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-dir"), default='data',
                metavar='data/ballesteros-2020/input', type="character",
                help=paste("directory of directories (one up) containing standard 10x files:",
                           "barcodes.tsv, genes.tsv, and matrix.mtx")),
    
    make_option(c("-o", "--output-dir"), default='data/output',
                metavar='output', type="character",
                help="set the output directory"),

    make_option(c("-f", "--figures-dir"), default='figures',
                metavar='figures', type="character",
                help="set the figures directory"),

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
gene <- opt[["gene-of-interest"]]

# file io
input_dir <- opt[['input-dir']]
output_dir <- opt[['output-dir']]
figures_dir <- opt[['figures-dir']]
multiplicity <- 'integrated'


# Start Log
start_time <- Sys.time()
log <- log_open(paste0("coral_reef-",
    strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Main

# read files
all_files <- file.path('data', list.files(
    file.path(wd, 'data'), pattern=paste0(tolower(gene), '.csv'),
    recursive = TRUE, full.names = FALSE
))
relevant_files <- grep(
    pattern=paste(
        '[[:alnum:][:space:]-]', 'output', 'expression', tolower(gene), multiplicity,
        sep='/'
    ), all_files, value = TRUE)
df <- append_many_csv(file.path(wd, relevant_files))


# extract metadata
df[['dataset']] <- sapply(df[['filepath']], 
    function(x) stringr::str_match(x, paste(wd, 'data', "(.*?)", "output*", sep='/'))[[2]]
)
df[['sample_name']] <- sapply(df[['filepath']], 
    function(x) stringr::str_match(basename(x),
        paste0('cell_type-integrated-', "(.*?)", '-', tolower(gene), '.csv'))[[2]]
)
df[['id']] <- paste(df[['dataset']], df[['sample_name']])  # x axis text


# plot
fig <- plot_dotplot(
    df[(df[['num_cells_total']] >= 50), ],  # filter
    title='Overview'
)
if (!troubleshooting) {
    savefig(file.path(wd, figures_dir, paste0('dotplot-overview.png')),
            height=1600, width=3200, dpi=800, 
            makedir=FALSE, troubleshooting=troubleshooting)    
}


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
