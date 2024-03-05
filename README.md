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
BiocManager::install('DESeq2')

# SeuratData
devtools::install_github('satijalab/seurat-data@d6a8ce6')
library('SeuratData')
SeuratData::options(timeout = 240)
SeuratData::InstallData("ifnb")
```

## Scripts

These are the available scripts.

| Script Number | Script Name | Description |
| :------------ | ----------- | ----------- |
| 1 | R/shoaling\_fish.R | This script performs basic qc filtering, integrates the data as defined by the config.json files,  and uses SingleR to label the clusters in an unbiased way. It outputs a variety of figures and the .RData file to be used downstream. |
| 2 | R/finding\_nemo.R | This script visualizes your gene of interest using the .RData file output by R/shoal\_formation.R. This is a relatively fast script. |
| 3 | R/arribada.R | Parallelizes shoal_formation, since this takes some time. |
| 4 | R/aquarium.R  | Shiny app. Still under development. |


## Copyright

This code is copyright by Harrison Wang in 2023.
