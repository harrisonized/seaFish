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
| 1 | R/finding\_nemo.R | Main script for analyzing scRNAseq dataset. A more detailed description is coming |
| 2 | R/arribada.R | Parallelize the data analysis. |
| 3 | R/aquarium.R  | Shiny app. Still under development. |


## Copyright

This code is copyright by Harrison Wang in 2023.
