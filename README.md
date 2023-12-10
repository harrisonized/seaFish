## seaFish

With high "efficiency," we use this pipeline to search for "a fish in the sea." The purpose of this repository is to enable me to efficiently perform exploratory data analysis on publicly available scRNAseq datasets. This can serve as a valuable hypothesis generator to check our intuition before we invest time and resources into experimentation.

This repository is currently undergoing major refactoring to allow it to handle multiple datasets.

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
