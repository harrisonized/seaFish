# library('Matrix')
wd = dirname(dirname(this.path::here()))
import::here(Matrix, 'readMM')
import::here(readr, 'read_tsv')
import::here(ggplot2, 'ggsave', 'last_plot')
import::here(grid, 'grid.newpage', 'grid.draw')
import::here(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'filter_list_for_match', .character_only=TRUE)

## Functions
## list_files
## load_rdata
## read_10x
## savefig


#' List all files with a specific extension
#' 
#' @description
#' This is a thin wrapper around [list.files()].
#' 
#' @references
#' \href{https://stackoverflow.com/questions/7187442/filter-a-vector-of-strings-based-on-string-matching}{StackOverflow post}
#' 
list_files <- function(dir_path, ext=NULL, recursive = TRUE) {
    all_files = list.files(dir_path, recursive = recursive, full.name=TRUE)

    if (!is.null(ext)) {
        # See: https://stackoverflow.com/questions/7187442/filter-a-vector-of-strings-based-on-string-matching
        return (all_files[tools::file_ext(all_files)==ext])
    } else {
        return (all_files)
    }
}


#' Loads an RData file and allows you to store it in a chosen variable
#' 
#' @description
#' Without this function, base R uses the filename as the variable name
#' 
#' @references
#' \href{https://stackoverflow.com/questions/5577221/can-i-load-a-saved-r-object-into-a-new-object-name}{StackOverflow post}
#' 
load_rdata <- function(filepath){
    load(filepath)
    return( get(ls()[ls() != "filepath"]) )
}


#' Read 10X
#'
#' @description Alternative to [Seurat::Read10X()] that enables you to specify the filenames
#' 
#' @references
#'\href{https://github.com/satijalab/seurat/issues/4096}{Github issue}
#'
#' @seealso [Seurat::Read10X()]
#' 
read_10x <- function(
    data_dir,
    matrix_file='matrix.mtx',
    genes_file='genes.tsv',
    barcodes_file='barcodes.tsv'
) {

    if (!file.exists(file.path(data_dir, matrix_file)) |
        !file.exists(file.path(data_dir, genes_file)) |
        !file.exists(file.path(data_dir, barcodes_file))) {

        filenames = basename(list_files(data_dir))
        matrix_file = filter_list_for_match(filenames, 'matrix')
        genes_file = unlist(lapply(
            c('genes', 'features'),
            function(pattern) filter_list_for_match(filenames, pattern))
        )
        barcodes_file = filter_list_for_match(filenames, 'barcodes')
    }

    expr_mtx <- Matrix::readMM(file.path(data_dir, matrix_file))
    genes <- read_tsv(file.path(data_dir, genes_file), col_names=FALSE, show_col_types = FALSE)
    barcodes <- read_tsv(file.path(data_dir, barcodes_file), col_names=FALSE, show_col_types = FALSE)

    colnames(expr_mtx) <- barcodes[['X1']]  # barcode sequence
    rownames(expr_mtx) <- genes[['X2']]  # gene names

    # Return dgCMatrix instead of dgTMatrix
    # See: https://slowkow.com/notes/sparse-matrix/#the-triplet-format-in-dgtmatrix
    expr_mtx <- as(expr_mtx, "CsparseMatrix")

    # Note: the above is equivalent to this, but is more explicit
    # expr_mtx <- Matrix::ReadMtx(
    #     mtx=file.path(data_dir, matrix_file),
    #     cells=file.path(data_dir, barcodes_file),
    #     features=file.path(data_dir, genes_file),
    #     feature.column=2
    # )

    return(expr_mtx)
}


#' Save Figure
#' 
#' @description Switch case to reduce the number of lines in the main script
#' 
savefig <- function(
    filepath,
    fig=NULL,
    height=800, width=1200, dpi=300, units="px", scaling=0.5,
    makedir=FALSE,
    troubleshooting=FALSE,
    lib='ggplot',  # choose: ggplot, grid
    default_ext = '.png'
) {
    if (!troubleshooting) {

        # make directory
        dirpath <- dirname(filepath)
        if (makedir && !dir.exists(dirpath)) {
            dir.create(dirpath, recursive=TRUE)
        }

        # add file extension if not included
        if (tools::file_ext(filepath)=='') {
            filepath <- paste0(filepath, default_ext)
        }

        if (lib=='ggplot') {
            if (!inherits(fig, "ggplot")) { fig <- last_plot()}

            withCallingHandlers({
                ggsave(
                    filepath,
                    plot=fig,
                    height=height, width=width, dpi=dpi, units=units, scaling=scaling
                )
            }, warning = function(w) {
                if ( any(grepl("rows containing non-finite values", w),
                         grepl("fewer than two data points", w)) ) {
                    invokeRestart("muffleWarning")
                }
            })

        } else if (lib=='grid') {
            png(filepath,
                height=height, width=width, res=dpi, units=units
            )
            grid.newpage()
            grid.draw(fig$gtable)
            dev.off()
        } else {
            warning(paste0("lib='", lib, "' not found"))
        }
    }
}
