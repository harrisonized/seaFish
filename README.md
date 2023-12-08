## seaFish

With high "efficiency," we use this pipeline to search for "a fish in the sea." The purpose of this repository is to enable me to efficiently perform exploratory data analysis on publicly available scRNAseq datasets. This can serve as a valuable hypothesis generator to check our intuition before we invest time and resources into experimentation.

This repository is currently undergoing major refactoring to allow it to handle multiple datasets.

## Installation

Install the following packages in R:

```R
# regular packages
install.packages('Matrix')
install.packages('dplyr')
install.packages('shinydashboard')
BiocManager::install('DropletUtils')
```

## Scripts

These are the available scripts, to be run in the order listed.

| Script Number | Script Name | Description |
| :--- | ------ | ----------- |
| 1 | R/single_analysis.R |  |
| 2 | R/scripts/integrated_analysis.R  |  |
| 3 | R/scripts/plot\_single\_umap.sh |  |
| 4 | R/scripts/assign_clusters.sh |  |

How this should go:
Integrate if necessary, otherwise don't.


## Copyright

This code is copyright by Harrison Wang in 2023.
