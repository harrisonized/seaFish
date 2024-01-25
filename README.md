## seaFish

With high "efficiency," we use this pipeline to search for "a fish in the sea." The purpose of this repository is to enable me to efficiently perform exploratory data analysis on publicly available scRNAseq datasets. This can serve as a valuable hypothesis generator to check our intuition before we invest time and resources into experimentation. This package uses Seurat v4, because Seurat v5 is still fairly new as of January 2024.

## Installation

Install the following packages in R:

```R
# regular packages
install.packages('dplyr')
install.packages('Matrix')
BiocManager::install('Seurat')
BiocManager::install('SingleR')
BiocManager::install('celldex')
BiocManager::install('DropletUtils')

# SeuratData
devtools::install_github('satijalab/seurat-data@d6a8ce6')
library('SeuratData')
SeuratData::options(timeout = 240)
SeuratData::InstallData("ifnb")
```

## Scripts

These are the available scripts, to be run in the order listed.

| Script Number | Script Name | Description |
| :--- | ------ | ----------- |
| 1 | R/cluster_samples.R | Individually clusters each sample in a dataset. |
| 2 | R/cluster_groups.R  | Clusters groups of samples in a dataset based on the config.json file. |
| 3 | R/shiny_umap.sh | Takes saved RData and plots them. |


## Copyright

This code is copyright by Harrison Wang in 2023.
