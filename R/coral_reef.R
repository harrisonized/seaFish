## Get an overview of all the datasets in the form of a dotplot

wd = dirname(this.path::here())  # wd = '~/github/R/seaFish'
# suppressPackageStartupMessages(library('ggplot2'))
library('optparse')
library('logr')

import::from(rjson, 'fromJSON', 'toJSON')
import::from(stringr, 'str_match')
import::here(tidyr, 'pivot_wider')
import::from(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'items_in_a_not_b', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    'append_many_csv', 'savefig', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'df_tools.R'),
    'set_index', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'plotting.R'),
    'plot_dotplot', 'plot_heatmap', .character_only=TRUE)


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(
    make_option(c("-i", "--input-dir"), default='data',
                metavar='data/ballesteros-2020/input', type="character",
                help=paste("directory of directories (one up) containing standard 10x files:",
                           "barcodes.tsv, genes.tsv, and matrix.mtx")),
    
    make_option(c("-o", "--output-dir"), default='data/output',
                metavar='output', type="character",
                help="set the output directory"),

    make_option(c("-f", "--figures-dir"), default='figures',
                metavar='figures', type="character",
                help="set the figures directory"),

    make_option(c("-j", "--config"), default='config.json',
                metavar='config.json', type="character",
                help="json file containing group information"),

    make_option(c("-l", "--height"), default=4200,
                metavar="3600", type="integer",
                help="height in px"),

    make_option(c("-w", "--width"), default=7000,
                metavar="7000", type="integer",
                help="width in px"),

    make_option(c("-r", "--order"), default='xlabel_order.csv',
                metavar='xlabel_order.csv', type="character",
                help="csv file containing a list"),

    make_option(c("-g", "--gene-of-interest"), default="Dnase1l1",
                metavar="Dnase1l1", type="character",
                help="choose a gene"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
troubleshooting <- opt[['troubleshooting']]
gene <- opt[["gene-of-interest"]]

# file io
input_dir <- opt[['input-dir']]
output_dir <- opt[['output-dir']]
figures_dir <- opt[['figures-dir']]

# Start Log
start_time <- Sys.time()
log <- log_open(paste0("coral_reef-",
    strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))


# ----------------------------------------------------------------------
# Parse config

config_file <- file.path(wd, opt[['input-dir']], opt[['config']])
if (file.exists(config_file)) {
    log_print(paste(Sys.time(), 'Reading config...'))
    config <- fromJSON(file=config_file)
    group_names <- names(config)
} else {
    log_print(paste(Sys.time(), 'Generating config...'))
    group_names <- list.dirs(
        file.path(wd, opt[['input-dir']]),
        full.names=FALSE, recursive=FALSE
    )
    group_names <- items_in_a_not_b(group_names, c('', 'archive', 'broken', 'data'))  # filter empty
    config <- as.list(setNames(rep('integrated', length(group_names)), group_names))
    config[['_ignore']] <- list()
    if (!troubleshooting) {
        write(toJSON(config), file=config_file)
    }
}

# split out ignore_list
ignore_list <- config[['_ignore']]
config <- config[1:length(config)-1]


# find files based on config
search_dirs <- sapply(
    names(config), function(x) 
    file.path(x, 'output/expression', tolower(gene), config[[x]]),
    USE.NAMES=TRUE
)
all_files <- unlist(lapply(
    file.path(wd, 'data', search_dirs), function(x)
        list.files(
        x, pattern=paste0(tolower(gene), '.csv'),
        recursive = FALSE, full.names = TRUE)
))
# remove ignored files
relevant_files <- names(items_in_a_not_b(
    setNames(basename(all_files), all_files),
    ignore_list
))


# ----------------------------------------------------------------------
# Read Data

raw_data <- append_many_csv(relevant_files)

# extract metadata
raw_data[['dataset']] <- sapply(raw_data[['filepath']], 
    function(x) stringr::str_match(x, paste('data', "(.*?)", "output*", sep='/'))[[2]]
)
raw_data[['sample_name']] <- sapply(raw_data[['filepath']], 
    function(x) stringr::str_match(basename(x),
        paste0('(cell_type-integrated|cell_type)-', "(.*?)", '-', tolower(gene), '.csv'))[[3]]
)
raw_data[['id']] <- paste(raw_data[['dataset']], raw_data[['sample_name']])  # x axis text


# ----------------------------------------------------------------------
# Dotplot

df <- raw_data[(raw_data[['num_cells_total']] >= 50), ]  # quality filter

# xlabel order
if (file.exists(file.path(wd, input_dir, opt[['order']]))) {
    xlabel_order <- read.csv(file.path(wd, 'data', opt[['order']]), header=FALSE)[['V1']]
    df <- df[(df[['id']] %in% xlabel_order), ]
    df[['id']] <- factor(df[['id']], levels = xlabel_order)
}

# ylabel order
ylabel_order <- c(
    "Macrophages", "Monocytes", "Neutrophils", "DC", "Basophils",  "Mast cells", "Stem cells",
    "B cells",  "B cells, pro", "T cells", "Tgd", "NK cells", "NKT", "ILC",
    "Endothelial cells", "Fibroblasts", "Epithelial cells", "Stromal cells"
)
df[['cell_type']] <- factor(df[['cell_type']], levels = ylabel_order)

# plot
fig <- plot_dotplot(
    df,
    title=expression(paste(italic("Dnase1l1"), " Expression"))
)

if (!troubleshooting) {
    savefig(file.path(wd, figures_dir, paste0('dotplot-overview.png')),
            height=3200, width=5000, dpi=320, scaling=2,
            makedir=FALSE, troubleshooting=troubleshooting)
}


# ----------------------------------------------------------------------
# Heatmap

tmp <- as.data.frame(pivot_wider(
    df,
    id_cols = 'cell_type',
    names_from='id',
    values_from='num_cells_total',
    values_fill = 0,
    names_sort=FALSE
))
tmp <- set_index(tmp, 'cell_type')
fig <- plot_heatmap(
    tmp, title='Number of Cells',
    xlabel_order=xlabel_order,
    annotations=TRUE, xaxis_angle=60
)
if (!troubleshooting) {
    savefig(file.path(wd, figures_dir, paste0('heatmap-overview.png')),
            height=1600+200, width=3200, dpi=500, 
            makedir=FALSE, troubleshooting=troubleshooting)    
}


end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()

quit()


# ----------------------------------------------------------------------
# Scratch

xlabel_order <- c(
    "ballesteros-2020 bm",  # bone marrow
    "ballesteros-2020 lung",
    "ballesteros-2020 spleen",
    "ballesteros-2020 pb",
    "edwards-2020-muscle-injury d0",  # muscle
    "edwards-2020-muscle-injury d3",
    "justynski-2023-apoptosis-skin 24hr_wound_beds",  # skin
    "justynski-2023-apoptosis-skin 48hr_wound_beds",

    # lung
    "moore-2023-lung-lps Homeostasis CD45",
    "moore-2023-lung-lps LPS Day 3 CD45",
    "moore-2023-lung-lps LPS Day 6 CD45",
    "tomlinson-2023-staph LAC-infected WT", 
    "tomlinson-2023-staph PBS-treated WT",
    "kasmani-2023-flu Day_3_Young",  
    "kasmani-2023-flu Day_9_Young",
    "curras-alonso-2023-radiation Ctrl_1",
    "curras-alonso-2023-radiation IR_17Gy_1M_1",
    "curras-alonso-2023-radiation IR_17Gy_3M_1",
    "curras-alonso-2023-radiation IR_17Gy_5M_1",

    # CD45+ spleen
    "dou-2024-pristane Wt_F",  
    "dou-2024-pristane Wt_M",
    "murao-2024-sepsis sham spleen",
    "murao-2024-sepsis CLP spleen",

    # peritoneum
    "murao-2024-sepsis sham peritoneum",
    "murao-2024-sepsis CLP peritoneum",
    "santeford-2021-thioglycollate Young_control",
    "santeford-2021-thioglycollate Old_control",
    "lantz-2020-efferocytosis WT_Ctrl",
    "lantz-2020-efferocytosis WT_2_hr",
    "lantz-2020-efferocytosis WT_6_hr",

    # kidney
    "rudman-melnick-2020-aki control",  
    "rudman-melnick-2020-aki day1",
    "rudman-melnick-2020-aki day2",
    "rudman-melnick-2020-aki day4",
    "rudman-melnick-2020-aki day7_2",
    "rudman-melnick-2020-aki day11_2",
    "rudman-melnick-2020-aki day14_2"
)

if (!troubleshooting) {
    write.table(xlabel_order,
        file.path(wd, 'data', paste0('_', opt[['order']])),
        sep=',', row.names=FALSE, col.names=FALSE
    )
}
