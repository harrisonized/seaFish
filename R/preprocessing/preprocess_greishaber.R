## Clean up data before they go into the pipeline

wd = dirname(dirname(this.path::here()))  # wd = '~/github/R/seaFish'
library('optparse')
library('logr')
import::from(reshape2, 'melt')

import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    'read_counts_data', 'write_gz', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'df_tools.R'),
    'reset_index', .character_only=TRUE)

# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-dir"), default='data/greishaber-2021-neutrotime/downloads',
                metavar='data/greishaber-2021-neutrotime/downloads', type="character",
                help=paste("directory of files to be converted")),
    
    make_option(c("-o", "--output-dir"), default='output',
                metavar='output', type="character",
                help="set the output directory"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
troubleshooting <- opt[['troubleshooting']]

# file io
output_dir <- file.path(opt[['input-dir']], opt[['output-dir']])

# Start Log
start_time <- Sys.time()
log <- log_open(paste0("preprocess_greishaber-",
                       strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))

# ----------------------------------------------------------------------
# Main

log_print(paste(Sys.time(), 'Reading data...'))

metadata <- read.delim(file.path(wd, opt[['input-dir']],
    "GSM5029341_inflammation_dataset3_readme.txt.gz")
)
partitions <- melt(table(metadata[, c('tissue', 'stimulus')]))
partitions <- reset_index(partitions[(partitions['value'] != 0), ], drop=TRUE)
counts <- read_counts_data(
    file.path(wd, opt[['input-dir']], "GSM5029341_inflammation_dataset3.txt.gz"),
    as.matrix=FALSE
)

log_print(paste(Sys.time(), 'Writing data...'))
log_print(paste(Sys.time(), 'Found', nrow(partitions), 'partitions'))

for (i in 1:nrow(partitions)){

    tissue <- partitions[i, 'tissue']
    stimulus <- partitions[i, 'stimulus']

    log_print(paste(Sys.time(), stimulus, tissue))

    barcodes <- metadata[
        (metadata[['tissue']]==tissue) &
        (metadata[['stimulus']]==stimulus),
        'barcode'
    ]

    write_gz(
        counts[, barcodes],
        file.path(wd, output_dir, paste0(tissue, '-', stimulus),
                  paste0('GSM5029341_', tissue, '_', stimulus, '.txt'))
    )
}

end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
