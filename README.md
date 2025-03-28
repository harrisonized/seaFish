## seaFish

The purpose of this repository is to enable me to efficiently perform exploratory data analysis on publicly available scRNAseq datasets. This can serve as a valuable hypothesis generator to check our intuition before we invest time and resources into experimentation. This package currently uses Seurat v4, but will be updated to use Seurat v5.


## Scripts

| Script Number | Script Name | Description |
| :------------ | ----------- | ----------- |
| 1 | integrate\_cluster\_label\_v4.R | This script performs basic qc filtering, integrates the data as defined by the config.json files,  and uses SingleR to label the clusters in an unbiased way. It outputs a variety of figures and the .RData file to be used downstream. By default, this script assumes the input data is mouse. To switch to human, use the `-f` flag.  |
| 2 | plot\_gene\_of\_interest.R | This script visualizes the gene of interest using the .RData file output by integrate\_cluster\_label\_v4.R. This is a relatively fast script. |
| 3 | parallelize_icl.R | Parallelizes integrate\_cluster\_label\_v4.R, since this takes some time. |
| 4 | plot\_overview.R | Creates a dotplot visualization to provide an overview of all the processed datasets. Still under development. |
| 5 | dashboard.R  | Shiny app. Very basic visualizations at the moment. At some point once I finalize what figures I consider to be standard, I will build this out. |

## How To Use

Set up your data directory, then run each script from the command line:

```bash
cd path/to/seaFish
Rscript R/parallelize_icl.R
Rscript R/plot_gene_of_interest.R -v
```

Alternatively, if you are troubleshooting, you can copy/paste each line into RStudio.

## Package Requirements

I developed this package on an Apple M1 Pro. In order to compile some packages, such as Matrix, I had to install [gfortran-6.1](https://cran.r-project.org/bin/macosx/tools/). In order to use this, please check that you have the following packages.

| Package | Version | Install Command          | Comments |
| :------ | ------- | ------------------------ | -------- |
| dplyr   | 1.1.2   | install.packages()       | |
| Seurat  | 4.3.0.1 | remotes::install_version() | |
| Matrix  | 1.6-1   | remotes::install_version() | This is required for RunUMAP(). Higher versions break the `irlba` package. |
| SingleR | 2.2.0   | BiocManager::install()   | |
| celldex | 1.10.1  | BiocManager::install()   | |
| DropletUtils | 1.10.1 | BiocManager::install() | |
| DESeq2  | 1.40.2  | BiocManager::install()   | |
| SeuratData | 0.2.2 | devtools::install_github('satijalab/seurat-data@d6a8ce6') | |
| ifnb.SeuratData | 3.1.0 | SeuratData::InstallData("ifnb") | If the code exits out, set SeuratData::options(timeout = 240). |
| EnhancedVolcano | 1.18.0 | BiocManager::install() | |


## Copyright

This code is copyright by Harrison Wang in 2023.
