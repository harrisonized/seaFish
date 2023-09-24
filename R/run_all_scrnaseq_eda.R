wd = dirname(this.path::here())  # wd = '~/github/R/harrisonRTools'
library('optparse')
library('logr')
source(file.path(wd, 'R', 'functions', 'file_io.R'))


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-d", "--data-dir"),
                default='data/scrnaseq-ballesteros/10x_counts',
                metavar='data/scrnaseq-ballesteros/10x_counts', type="character",
                help="set folder containing all the 10x_counts directories"),

    make_option(c("-o", "--output-dir"),
                default="figures/scrnaseq-ballesteros",
                metavar="figures/scrnaseq-ballesteros", type="character",
                help="set the output directory for the figures"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
troubleshooting = opt[['troubleshooting']]
data_dir = opt[['data-dir']]
fig_dir = opt[['output-dir']]

# Start Log
start_time = Sys.time()
log <- log_open(paste0("run_all_scrnaseq_eda-",
                       strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Run scirpt

log_print(paste(Sys.time(), 'Data dir...', data_dir))

dirnames = basename(list_files(file.path(wd, data_dir), recursive=FALSE))

for (dirname in dirnames){
    log_print(paste(Sys.time(), 'Processing...', dirname))
    
    cmd = paste(
        'Rscript', file.path(wd, 'R', 'scrnaseq_eda.R'),
        '-i', file.path(data_dir, dirname, 'counts'),
        '-o', paste0(fig_dir, "/", dirname)
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
