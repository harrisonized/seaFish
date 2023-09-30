wd = dirname(dirname(this.path::here()))  # wd = '~/github/R/harrisonRTools'
library('optparse')
library('logr')
source(file.path(wd, 'R', 'functions', 'file_io.R'))


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
troubleshooting = opt[['troubleshooting']]

# Start Log
start_time = Sys.time()
log <- log_open(paste0("run_all_scrnaseq_assign_clusters-",
                       strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Run scirpt


io_df <- data.frame(
    'input_file' = c(
        'data/scrnaseq-ballesteros/integrated/bm/integrated_seurat.RData',
        'data/scrnaseq-ballesteros/integrated/lung/integrated_seurat.RData',
        'data/scrnaseq-ballesteros/integrated/pbzt/integrated_seurat.RData',
        'data/scrnaseq-ballesteros/integrated/spleen/integrated_seurat.RData'
    ),
    'output_dir' = c(
        "figures/scrnaseq-ballesteros/integrated/bm",
        "figures/scrnaseq-ballesteros/integrated/lung",
        "figures/scrnaseq-ballesteros/integrated/pbzt",
        "figures/scrnaseq-ballesteros/integrated/spleen"
    ),
    stringsAsFactors = FALSE
)

for (i in 1:nrow(io_df)) {

    dirname = basename(io_df[i, 'output_dir'])
    log_print(paste(Sys.time(), 'Processing...', dirname))

    cmd = paste(
        'Rscript', file.path(wd, 'R', 'scrnaseq_assign_clusters.R'),
        '-i', io_df[i, 'input_file'],
        '-o', io_df[i, 'output_dir']
    )

    # run
    if (!troubleshooting) {
        system(cmd)
    } else {
        print(cmd)
    }

}

end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
