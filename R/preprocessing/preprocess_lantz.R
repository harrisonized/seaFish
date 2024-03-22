## Clean up data before they go into the pipeline

wd = dirname(dirname(this.path::here()))  # wd = '~/github/R/seaFish'
library('optparse')
library('logr')

import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    'read_counts_data', 'write_gz', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'df_tools.R'),
    'reset_index', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'dict_zip', .character_only=TRUE)

# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-dir"), default='data/lantz-2020-efferocytosis/downloads',
                metavar='data/greishaber-2021-neutrotime/downloads', type="character",
                help=paste("directory of files to be converted")),

    make_option(c("-m", "--metadata"), default='metadata.csv',
                metavar='metadata.csv', type="character",
                help=paste("file that matches index to group")),

    make_option(c("-o", "--output-dir"), default='input',
                metavar='input', type="character",
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

metadata <- read.csv(file.path(wd, opt[['input-dir']], opt[['metadata']]))
id_to_name <- dict_zip(metadata[['id']], metadata[['sample_name']])
counts <- read_counts_data(
    file.path(wd, opt[['input-dir']], "GSE156234_aggregated_raw_counts.tsv.gz"),
    as.matrix=FALSE
)

for (id in 1:nrow(metadata)) {
    name <- id_to_name[[as.character(id)]]
    log_print(paste(Sys.time(), 'Processing...', name))

    barcodes <- colnames(counts)[grepl(id, colnames(counts))]
    tmp <- counts[, barcodes]
    tmp <- tmp[rowSums(tmp)>0, ]

    write_gz(
        tmp,
        file.path(wd, output_dir, name,
                  paste0('GSE156234_', name, '.txt'))
    )
}

end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
