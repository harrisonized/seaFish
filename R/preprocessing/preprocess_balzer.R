## Clean up data before they go into the pipeline
# sample_names were hardcoded in the following script:
# https://github.com/ms-balzer/IRI_adaptive_maladaptive_kidney_regeneration/blob/main/R/Seurat__allcells.R

wd = dirname(dirname(this.path::here()))  # wd = '~/github/R/seaFish'
library('optparse')
library('logr')
import::from(rjson, 'fromJSON')

import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    'write_gz', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'index_last_occurrence', .character_only=TRUE)

# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-dir"), default='data/balzer-2022-aki/downloads',
                metavar='data/balzer-2022-aki/downloads', type="character",
                help=paste("directory of files to be converted")),
    
    make_option(c("-o", "--output-dir"), default='input',
                metavar='input', type="character",
                help="set the output directory"),

    make_option(c("-j", "--config"), default='config.json',
                metavar='config.json', type="character",
                help="json file containing group information"),

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
log <- log_open(paste0("preprocess_balzer-",
                       strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Config File

config_file <- file.path(wd, opt[['input-dir']], opt[['config']])
name_from_id <- list2env(fromJSON(file=config_file))


# ----------------------------------------------------------------------
# Main

log_print(paste(Sys.time(), 'Reading data...'))
counts <- readRDS(file.path(wd, opt[['input-dir']],
    "GSE180420_EXPORTinhib_counts.rds")
)

# barcodes are provided in the following format: "norm1_AAACCTGAGATATGCA"
prefixes <- sapply(colnames(counts), function(x) strsplit(x, split='_')[[1]][[1]])
col_indexes <- index_last_occurrence(prefixes)

# split the data up
for (i in 1:length(col_indexes)) {
    sample_id <- names(col_indexes)[i]
    sample_name <- name_from_id[[sample_id]]

    start <- ifelse(i==1, 1, col_indexes[[i-1]]+1)
    stop <- col_indexes[[i]]
    tmp <- counts[, start:stop, drop=FALSE]

    # check uniqueness
    prefixes_found <- unique(sapply(colnames(tmp),
        function(x) strsplit(x, split='_')[[1]][[1]]
    ))
    log_print(paste(Sys.time(), 'Processing...', paste(prefixes_found, collapse=', ' )))
    
    log_print(paste(Sys.time(), 'Writing data...'))
    write_gz(
        tmp,
        file.path(wd, output_dir, sample_name,
                  paste0('GSE180420_EXPORT_counts_', sample_name, '.txt'))
    )
}

end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
