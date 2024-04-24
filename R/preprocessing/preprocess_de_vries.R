## Clean up data before they go into the pipeline

wd = dirname(dirname(this.path::here()))  # wd = '~/github/R/seaFish'
library('optparse')
library('logr')
suppressPackageStartupMessages(library('biomaRt'))
import::from(Seurat, "UpdateSeuratObject")
import::from(rjson, 'fromJSON')

import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    'write_gz', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'multiple_replacement', .character_only=TRUE)

# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-dir"), default='data/de-vries-2020-candida-human/downloads',
                metavar='data/de-vries-2020-candida-human/downloads', type="character",
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
output_dir <- multiple_replacement(opt[['input-dir']], c('downloads'=opt[['output-dir']]))

# Start Log
start_time <- Sys.time()
log <- log_open(paste0("preprocess_balzer-",
                       strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Gene labels

# reqiures an internet connection
ensembl <- useEnsembl("ensembl", "hsapiens_gene_ensembl")
gene_id_map <- getBM(
    c("ensembl_gene_id", "external_gene_name"),
    mart=ensembl,
    useCache=TRUE
)

# ----------------------------------------------------------------------
# Main

log_print(paste(Sys.time(), 'Reading data...'))
seurat_obj <- readRDS(file.path(wd, opt[['input-dir']],
    "seurat_object_final_anonymized.Rds")
)
seurat_obj <- UpdateSeuratObject(seurat_obj)

# split seurat obj
sample_names <- unique(seurat_obj@meta.data[['orig.ident']])
for (sample_name in sample_names) {
    log_print(paste(Sys.time(), 'Processing...', sample_name))
    seurat_obj_subset <- subset(seurat_obj, subset = orig.ident == sample_name)
    
    counts <- seurat_obj_subset@assays$RNA@counts
    ensembl_gene_id <- rownames(counts)
    df <- data.frame(ensembl_gene_id)  # df with 'ensemble_gene_id' col

    # get gene names from gene_id_map
    df <- merge(
        df, gene_id_map,
        by.x = 'ensembl_gene_id',
        by.y = 'ensembl_gene_id',
        all.x = TRUE, all.y = FALSE,
        sort = FALSE
    )
    df[['gene_name']] <- ifelse(
        is.na(df[['external_gene_name']]) | (df[['external_gene_name']]==''),
        df[['ensembl_gene_id']],
        df[['external_gene_name']]
    )
    rownames(counts) <- df[['gene_name']]

    # write
    if (!troubleshooting) {
        write_gz(
            counts,
            file.path(wd, output_dir, sample_name,
                      paste0('counts_', sample_name, '.txt'))
        )  
    }
}

end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
