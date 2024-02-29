## Every fall, thousands of sea turtles arrive on the shores of Costa Rica to
## synchronously lay millions of eggs in an event known as the arribada,
## Spanish for "arrival by sea." This script is an orchestrator that
## synchronously runs the other scripts.

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

    make_option(c("-g", "--gene-of-interest"), default="Dnase1l1",
                metavar="Dnase1l1", type="character",
                help="choose a gene"),

    make_option(c("-s", "--slice"), default="",
                metavar="", type="character",
                help="enter comma separated list of indices"),

    make_option(c("-l", "--load-savepoint"), default=FALSE,
                metavar="FALSE", action="store_true", type="logical",
                help="use this to skip reading in the data and doing the data integration, use this for testing"),

    make_option(c("-m", "--method"), default="collate",
                metavar="collate", type="character",
                help="choose from 'collate' or 'chunk'"),

    make_option(c("-n", "--num"), default=4,
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

# Start Log
start_time <- Sys.time()
log <- log_open(paste0("arribada-",
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

if (opt[['method']]=='collate') {
    grouped_idxs <- collate(idxs, opt[['num']])
} else if (opt[['method']]=='chunk') {
    grouped_idxs <- chunker(idxs, opt[['num']])
} else {
    stop("Choose a valid method: 'collate' or 'chunk'")
}


# ----------------------------------------------------------------------
# Generate commands

log_print(paste(Sys.time(), 'Groups found...', paste(names(config), collapse=', ' )))

commands <- new.env()
for (idx in 1:length(grouped_idxs)) {

    command <- paste(
        'Rscript R/finding_nemo.R',
            '-i', opt[['input-dir']],
            '-o', opt[['output-dir']],
            '-j', opt[['config']],
            '-c', opt[['celldex']],
            ifelse(opt[['ensembl']], '-e', ''),
            '-g', opt[['gene-of-interest']],
            '-s', '"', as.character(paste0(grouped_idxs[idx])), '"',
            ifelse(opt[['load-savepoint']], '-l', ''),
            ifelse(opt[['troubleshooting']], '-t', '')
    )
    log_print(command)
    commands[[as.character(idx)]] <- command
}
commands <- as.list(commands)


# ----------------------------------------------------------------------
# Run in parallel

if (!troubleshooting) {
    invisible(
        mclapply(commands, function(x) system(x), mc.cores = min(n_cores, length(commands)))
    )
}


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
