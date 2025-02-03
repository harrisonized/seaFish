## Find your gene of interest

wd = dirname(this.path::here())  # wd = '~/github/R/seaFish'
suppressPackageStartupMessages(library('Seurat'))
library('zeallot')  # %<-%
library('optparse')
library('logr')
import::from(stringr, 'str_extract')
import::from(rjson, 'fromJSON')

import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    'load_rdata', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'multiple_replacement', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'text_tools.R'),
    'random_hash', .character_only=TRUE)
import::from(file.path(wd, 'R', 'functions', 'export_results.R'),
    'export_gene_of_interest', .character_only=TRUE)


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-dir"), default='data/ballesteros-2020/input',
                metavar='data/ballesteros-2020/input', type="character",
                help=paste("directory of directories (one up) containing standard 10x files:",
                           "barcodes.tsv, genes.tsv, and matrix.mtx")),
    
    make_option(c("-o", "--output-dir"), default='output',
                metavar='output', type="character",
                help="set the output directory"),

    make_option(c("-j", "--config"), default='config.json',   # figure out how to deal with this
                metavar='config.json', type="character",
                help="json file containing group information"),

    make_option(c("-g", "--gene-of-interest"), default="Dnase1l1",
                metavar="Dnase1l1", type="character",
                help="choose a gene"),

    make_option(c("-s", "--slice"), default="",
                metavar="", type="character",
                help="enter comma separated list of indices"),

    make_option(c("-v", "--volcano"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="use this to draw volcano plots, which require additional time to compute"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
troubleshooting <- opt[['troubleshooting']]
gene <- opt[["gene-of-interest"]]

# file io
output_dirname <- opt[['output-dir']]
output_dir <- ifelse(
     grepl('input', opt[['input-dir']]),
     multiple_replacement(opt[['input-dir']], c('input'=output_dirname)),
     file.path(opt[['input-dir']], 'output')
)
figures_dir <- multiple_replacement(
    opt[['input-dir']], c('data'='figures', 'input'=output_dirname)
)


# Start Log
start_time <- Sys.time()
log <- log_open(paste0("finding_nemo-", random_hash(), '-',
                       strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Parse config

config_file <- file.path(wd, opt[['input-dir']], opt[['config']])
if (file.exists(config_file)) {
    log_print(paste(Sys.time(), 'Reading config...'))
    config <- fromJSON(file=config_file)
    group_names <- names(config)
    multiplicities <- sapply(config, function(x) ifelse(length(x)>1, 'integrated-', ''))
    filenames <- paste0("seurat_obj-", paste0(multiplicities, group_names), ".RData")
    names(filenames) <- group_names
} else {
    filenames <- list.files(file.path(wd, output_dir, 'rdata'))
    group_names <- sub('.RData', '',
        sapply(strsplit(filenames, '-'), function(x) tail(x, n=1))
    )
    names(filenames) <- group_names
}

# subset
if (opt[['slice']] != '') {
    slice <- eval(parse( text=paste0('c(', opt[['slice']], ')') ))
    c(config, group_names, filenames) %<-% list(config[slice], group_names[slice], filenames[slice])
}


# ----------------------------------------------------------------------
# Main

log_print(paste(Sys.time(), 'Groups found...', paste(names(config), collapse=', ' )))

for (group_name in group_names) {
    loop_start_time <- Sys.time()
    log_print(paste(loop_start_time, 'Processing group:', group_name))

    tryCatch({

        # ----------------------------------------------------------------------
        # Read Data

        filename <- filenames[[group_name]]
        filepath=file.path(wd, output_dir, 'rdata', filename)
        seurat_obj <- load_rdata(filepath=filepath)

        multiplicity <- str_extract(filename, "(integrated|individual)")
        multiplicity <- ifelse(is.na(multiplicity), 'individual', 'integrated')
        is_integrated <- ifelse(multiplicity=='integrated', TRUE, FALSE)
        prefix <- ifelse(is_integrated, 'integrated-', '')
        
        # ----------------------------------------------------------------------
        # Export Results

        if (is_integrated) {

            log_print(paste(Sys.time(), 'Exporting integrated results...'))

            export_gene_of_interest(
                seurat_obj,
                gene=gene,
                sample_name=group_name,
                multiplicity=multiplicity,
                data_dir=output_dir,
                figures_dir=figures_dir,
                file_basename=paste0(prefix, group_name),
                include_table=TRUE,
                include_volcano=opt[['volcano']],
                include_pseudo_bulk=TRUE,
                troubleshooting=troubleshooting,
                showfig=troubleshooting
            )

            log_print(paste(Sys.time(), 'Exporting integrated subset...'))

            # Keep significant populations (number of cells > 50)
            num_cells_per_label <- as.data.frame(table(seurat_obj$cell_type))  # value counts
            populations_to_keep <- num_cells_per_label[num_cells_per_label['Freq'] > 50, 'Var1']
            seurat_obj_subset <- seurat_obj[, seurat_obj$cell_type %in% populations_to_keep]

            export_gene_of_interest(
                seurat_obj_subset,
                gene=gene,
                sample_name=group_name,
                multiplicity=multiplicity,
                data_dir=output_dir,
                figures_dir=figures_dir,
                figures_subdir='subset',
                file_basename=paste0(prefix, group_name, '-subset'),
                include_table=FALSE,
                include_volcano=FALSE,
                include_pseudo_bulk=FALSE,
                troubleshooting=troubleshooting,
                showfig=troubleshooting
            )
        }


        # ----------------------------------------------------------------------
        # Export Individual Results

        log_print(paste(Sys.time(), 'Exporting individual results...'))

        seurat_objs <- SplitObject(seurat_obj, split.by = "sample_name")
        sample_names <- c(unique(seurat_obj[['sample_name']]))[[1]]
        for (sample_name in sample_names) {

            export_gene_of_interest(
                seurat_objs[[sample_name]],
                gene=gene,
                sample_name=sample_name,
                multiplicity='individual',
                data_dir=output_dir,
                figures_dir=figures_dir,
                file_basename=sample_name,
                include_table=TRUE,
                include_volcano=opt[['volcano']],
                include_pseudo_bulk=FALSE,
                troubleshooting=troubleshooting,
                showfig=troubleshooting
            )
        }
    },

    # pass
    error = function(condition) {
        log_print(paste("WARNING: GROUP", group_name, "NOT PROCESSED!!!"))
        log_print(paste("Error message: ", conditionMessage(condition)))
    })

    log_print(paste("Loop completed in:", difftime(Sys.time(), loop_start_time)))
}


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
