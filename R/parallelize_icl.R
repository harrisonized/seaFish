## This script is an orchestrator that synchronously runs integrate_cluster_label.R.

wd = dirname(this.path::here())  # wd = '~/github/R/seaFish'
library('parallel')
library('optparse')
library('logr')
import::from(rjson, 'fromJSON', 'toJSON')
import::from(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'multiple_replacement', 'collate', 'chunker', .character_only=TRUE)


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-dir"), default='data/ballesteros-2020/input',
                metavar='data/ballesteros-2020/input', type="character",
                help="directory of directories (one up) containing standard 10x files: barcodes.tsv, genes.tsv, and matrix.mtx"),
    
    make_option(c("-o", "--output-dir"), default='output',
                metavar='output', type="character",
                help="set the output directory"),

    make_option(c("-j", "--config"), default='config.json',
                metavar='config.json', type="character",
                help="json file containing group information"),
    
    make_option(c("-c", "--celldex"),
                default='ImmGen', metavar='ImmGen', type="character",
                help="Choose from: ['ENCODE', 'HPCA', 'DICE', 'ImmGen', 'Monaco', 'MouseRNAseq', 'Hemato']"),

    make_option(c("-e", "--ensembl"),  default=FALSE,
                metavar="FALSE", action="store_true", type="logical",
                help=paste("Celldex option. Logical scalar indicating whether to convert",
                           "row names to Ensembl IDs. Genes without a mapping to a",
                           "non-duplicated Ensembl ID are discarded.")),

    make_option(c("-m", "--markers"), default=FALSE,
                metavar="FALSE", action="store_true", type="logical",
                help="turn this on to plot heatmap of top markers"),

    make_option(c("-f", "--human"), default=FALSE,
                metavar="FALSE", action="store_true", type="logical",
                help="default organism is mouse, turn this on for human datasets"),

    make_option(c("-s", "--slice"), default="",
                metavar="", type="character",
                help="enter comma separated list of indices"),

    make_option(c("-z", "--group-option"), default="collate",
                metavar="collate", type="character",
                help="choose from 'collate' or 'chunk'"),

    make_option(c("-n", "--num-instances"), default=10,
                metavar="10", type="integer",
                help="enter number of instances, can go up to the number of cores"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
troubleshooting <- opt[['troubleshooting']]
gene <- opt[["gene-of-interest"]]

# Start Log
start_time <- Sys.time()
log <- log_open(paste0("parallelize_icl-",
                       strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# System

n_cores <- detectCores()
log_print(paste('Number of cores:', n_cores))


# ----------------------------------------------------------------------
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
if (opt[['slice']] != '') {
    idxs <- eval(parse( text=opt[['slice']] ))
    idxs <- idxs[idxs <= length(config)]
} else {
    idxs <- seq(1, length(config), by=1)
}

if (opt[['group-option']]=='collate') {
    grouped_idxs <- collate(idxs, opt[['num-instances']])
} else if (opt[['group-option']]=='chunk') {
    grouped_idxs <- chunker(idxs, opt[['num-instances']])
} else {
    stop("Choose a valid method: 'collate' or 'chunk'")
}


# ----------------------------------------------------------------------
# Generate commands

log_print(paste(Sys.time(), 'Groups found...', paste(names(config), collapse=', ' )))

commands <- new.env()
for (idx in 1:length(grouped_idxs)) {

    command <- paste(
        'Rscript R/sea_turtle.R',
            '-i', opt[['input-dir']],
            '-o', opt[['output-dir']],
            '-j', opt[['config']],
            '-c', opt[['celldex']],
            ifelse(opt[['ensembl']], '-e', ''),
            ifelse(opt[['markers']], '-m', ''),
            ifelse(opt[['human']], '-f', ''),
            '-s', '"', as.character(paste0(grouped_idxs[idx])), '"',
            ifelse(opt[['troubleshooting']], '-t', '')
    )
    log_print(command)
    commands[[as.character(idx)]] <- command
}
commands <- as.list(commands)


# ----------------------------------------------------------------------
# Run in parallel

if (!troubleshooting) {
    invisible(mclapply(
        commands, function(x) system(x),
        mc.cores = min( n_cores, length(commands) )
    ))
}


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
