## Functions
## run_standard_clustering

#' This is repeated across multiple analyses
#' See: https://satijalab.org/seurat/articles/integration_introduction.html#perform-an-integrated-analysis
#'
#' @export
run_standard_clustering <- function(seurat_obj, ndim=40) {
    seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
    seurat_obj <- RunPCA(seurat_obj, npcs = ndim, verbose = FALSE)
    seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:ndim)
    seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:ndim)
    seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
    return(seurat_obj)
}
