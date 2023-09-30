wd = dirname(dirname(this.path::here()))  # wd = '~/github/R/harrisonRTools'
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
log <- log_open(paste0("run_all_scrnaseq_integrated_analysis-",
                       strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Run scirpt

log_print(paste(Sys.time(), 'Data dir...', data_dir))


io_df <- data.frame(
    'inputs' = c(
        paste('data/scrnaseq-ballesteros/10x_counts/bm1/counts',
              'data/scrnaseq-ballesteros/10x_counts/bm2/counts', sep=':'),
        paste('data/scrnaseq-ballesteros/10x_counts/lung1/counts',
              'data/scrnaseq-ballesteros/10x_counts/lung2/counts', sep=':'),
        paste('data/scrnaseq-ballesteros/10x_counts/pbzt5/counts',
              'data/scrnaseq-ballesteros/10x_counts/pbzt13_1/counts',
              # 'data/scrnaseq-ballesteros/10x_counts/pbzt13_2/counts',  # PBZT13_2_genes has 3204 rows when 3388 expected
               sep=':'),  
        paste('data/scrnaseq-ballesteros/10x_counts/spleen1/counts',
              'data/scrnaseq-ballesteros/10x_counts/spleen2/counts', sep=':')
        ),
    'output_dir' = c(
        "data/scrnaseq-ballesteros/integrated/bm",
        "data/scrnaseq-ballesteros/integrated/lung",
        "data/scrnaseq-ballesteros/integrated/pbzt",
        "data/scrnaseq-ballesteros/integrated/spleen"
    ),

    'figures_dir' = c(
        "figures/scrnaseq-ballesteros/integrated/bm",
        "figures/scrnaseq-ballesteros/integrated/lung",
        "figures/scrnaseq-ballesteros/integrated/pbzt",
        "figures/scrnaseq-ballesteros/integrated/spleen"
    ),
    'labels' = c(
        "bm1:bm2",
        "lung1:lung2",
        # "pbzt5:pbzt13_1:pbzt13_2",  # PBZT13_2_genes has 3204 rows when 3388 expected
        "pbzt5:pbzt13_1",
        "spleen1:spleen2"
    ),
    stringsAsFactors = FALSE
)

for (i in 1:nrow(io_df)) {

    dirname = basename(io_df[i, 'output_dir'])
    log_print(paste(Sys.time(), 'Processing...', dirname))

    cmd = paste(
        'Rscript', file.path(wd, 'R', 'scrnaseq_integrated_analysis.R'),
        '-i', io_df[i, 'inputs'],
        '-o', io_df[i, 'output_dir'],
        '-f', io_df[i, 'figures_dir'],
        '-l', io_df[i, 'labels']
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
