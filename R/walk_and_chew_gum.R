## Parallelize integrate_data.R

wd = dirname(this.path::here())  # wd = '~/github/R/seaFish'
library('optparse')
library('logr')
import::from(rjson, 'fromJSON', 'toJSON')
import::from(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'chunker', 'collate', 'multiple_replacement', .character_only=TRUE)


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-dir"),
                default='data/ballesteros-2020/input',
                metavar='data/ballesteros-2020/input',
                type="character",
                help="directory of directories (one up) containing standard 10x files: barcodes.tsv, genes.tsv, and matrix.mtx"),
    
    make_option(c("-o", "--output-dir"),
                default='data/ballesteros-2020/output',
                metavar='data/ballesteros-2020/output',
                type="character",
                help="set the output directory"),

    make_option(c("-j", "--config"), default='config.json',
                metavar='config.json', type="character",
                help="json file containing group information"),
    
    make_option(c("-d", "--ndim"), default=30,
                metavar="30", type="integer",
                help="choose the number of dimensions for the data integration"),

    make_option(c("-c", "--celldex"),
                default='ImmGen', metavar='ImmGen', type="character",
                help="Choose from: ['ENCODE', 'HPCA', 'DICE', 'ImmGen', 'Monaco', 'MouseRNAseq', 'Hemato']"),

    make_option(c("-e", "--ensembl"),  default=FALSE,
                metavar="FALSE", action="store_true", type="logical",
                help="Used when querying celldex"),

    make_option(c("-g", "--gene-of-interest"), default="Dnase1l1",
                metavar="Dnase1l1", type="character",
                help="choose a gene"),

    make_option(c("-s", "--slice"), default="",
                metavar="", type="character",
                help="enter comma separated list of indices"),

    make_option(c("-m", "--method"), default="collate",
                metavar="collate", type="character",
                help="choose from 'collate' or 'chunk'"),

    make_option(c("-n", "--num"), default=2,
                metavar="4", type="integer",
                help="enter a number"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
troubleshooting <- opt[['troubleshooting']]
gene <- opt[["gene-of-interest"]]


# Parse config
config_file <- file.path(wd, opt[['input-dir']], opt[['config']])
if (file.exists(config_file)) {
    
    log_print(paste(Sys.time(), 'Reading config...'))
    config <- fromJSON(file=config_file)

} else {
    
    log_print(paste(Sys.time(), 'Generating config...'))
    
    group_names <- list.dirs(file.path(wd, opt[['input-dir']]), full.names=FALSE)
    group_names <- group_names[(group_names != '')]  # filter empty
    config <- as.list(setNames(group_names, group_names))
    
    if (!troubleshooting) {
        write(toJSON(config), file=config_file)
    }
}

# divide data
idxs <- seq(1, length(config), by=1)
if (opt[['method']]=='collate') {
    grouped_idxs <- collate(idxs, opt[['num']])
} else if (opt[['method']]=='chunk') {
    print(idxs)
    grouped_idxs <- chunker(idxs, opt[['num']])
}

# Start Log
start_time <- Sys.time()
log <- log_open(paste0("walk_and_chew_gum-",
                       strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))
log_print(paste(Sys.time(), 'Groups found...', paste(names(config), collapse=', ' )))


# ----------------------------------------------------------------------
# Integrate Data

for (idx in 1:length(grouped_idxs)) {

    if (troubleshooting) {
        run <- print
    } else {
        run <- system
    }

    run(paste(
        'Rscript R/integrate_data.R',
            '-i', opt[['input-dir']],
            '-o', opt[['output-dir']],
            '-j', opt[['config']],
            '-d', opt[['ndim']],
            '-c', opt[['celldex']],
            ifelse(opt[['ensembl']], '-e', ''),
            '-g', opt[['gene-of-interest']],
            '-s', '"', as.character(paste0(grouped_idxs[idx])), '"',
            ifelse(opt[['troubleshooting']], '-t', '')
    ))
}

end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
