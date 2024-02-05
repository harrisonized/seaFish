## Get the top markers from celldex data
## See: https://bioconductor.org/books/3.12/OSCA/marker-detection.html

wd = dirname(dirname(this.path::here()))  # wd = '~/github/R/seaFish'
library('optparse')
library('logr')
import::from(progress, 'progress_bar')
import::from(scran, 'findMarkers')
import::from(file.path(wd, 'R', 'tools', 'df_tools.R'),
    'reset_index', .character_only=TRUE)
import::from(file.path(wd, 'R', 'functions', 'single_cell.R'),
    'celldex_switch', .character_only=TRUE)


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-c", "--celldex"),
                default='ImmGen', metavar='ImmGen', type="character",
                help="Choose from: ['ENCODE', 'HPCA', 'DICE', 'ImmGen', 'Monaco', 'MouseRNAseq', 'Hemato']"),
    
    make_option(c("-o", "--output-dir"),
                default='data/celldex',
                metavar='data/celldex',
                type="character",
                help="set the output directory"),

    make_option(c("-e", "--ensembl"),  default=FALSE,
                metavar="FALSE", action="store_true", type="logical",
                help="When querying celldex"),


    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
troubleshooting <- opt[['troubleshooting']]

# Start Log
start_time <- Sys.time()
log <- log_open(paste0("cluster_samples-",
                       strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Main

log_print('Getting data...')

ref_data <- celldex_switch[[opt$celldex]](ensembl=opt[['ensembl']])  # SummarizedExperiment
markers <- findMarkers(ref_data, groups=ref_data[['label.main']])  # List of dataframes
cell_types <- names(markers)  # cell_type <- "Neutrophils"

log_print('Exporting data...')

pb <- progress_bar$new(total = length(cell_types))
for (cell_type in cell_types) {
    pb$tick()

    df <- reset_index(data.frame(markers[[cell_type]]), index_name='gene')  # subset
    df <- df[(df[['summary.logFC']] >= 2) & (df[['Top']] <= 12), ]  # filter
    df <- df[order(df[["summary.logFC"]], decreasing=TRUE), ]  # sort

    # save
    if (!troubleshooting) {
        if (!dir.exists(file.path(wd, opt['output-dir'], opt['celldex']))) {
            dir.create(file.path(wd, opt['output-dir'], opt['celldex']), recursive=TRUE)
        }
        write.table(
            df,
            file = file.path(wd, opt['output-dir'], opt['celldex'],
                             paste0('top_genes-', cell_type, '.csv')),
            row.names = FALSE, sep = ','
        )
    }
}

end_time <- Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
