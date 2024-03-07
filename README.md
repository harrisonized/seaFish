## seaFish

With high "efficiency," we use this pipeline to search for "a fish in the sea." The purpose of this repository is to enable me to efficiently perform exploratory data analysis on publicly available scRNAseq datasets. This can serve as a valuable hypothesis generator to check our intuition before we invest time and resources into experimentation. This package uses Seurat v4, because Seurat v5 is still fairly new as of January 2024.


## Scripts

| Script Number | Script Name | Description |
| :------------ | ----------- | ----------- |
| 1 | R/shoaling\_fish.R | This script performs basic qc filtering, integrates the data as defined by the config.json files,  and uses SingleR to label the clusters in an unbiased way. It outputs a variety of figures and the .RData file to be used downstream. |
| 2 | R/finding\_nemo.R | This script visualizes your gene of interest using the .RData file output by R/shoal\_formation.R. This is a relatively fast script. |
| 3 | R/arribada.R | Parallelizes shoal_formation, since this takes some time. |
| 4 | R/aquarium.R  | Shiny app. Very basic visualizations at the moment. At some point once I finalize what figures I consider to be standard, I will build this out. |


## Package Requirements

I developed this package on an Apple M1 Pro. In order to compile some packages, such as Matrix, I had to install [gfortran-6.1](https://cran.r-project.org/bin/macosx/tools/). In order to use this, please check that you have the following packages.

| Package | Version | Install Command          | Comments |
| :------ | ------- | ------------------------ | -------- |
| dplyr   | 1.1.2   | install.packages()       |
| Seurat  | 4.3.0.1 | remotes::install_version() | |
| Matrix  | 1.6-1   | remotes::install_version() | This is required for RunUMAP(). Higher versions break the `irlba` package.
| SingleR | 2.2.0   | BiocManager::install()   | |
| celldex | 1.10.1  | BiocManager::install()   | |
| DropletUtils | 1.10.1 | BiocManager::install() | |
| DESeq2  | 1.40.2  | BiocManager::install()   | |
| SeuratData | 0.2.2 | devtools::install_github('satijalab/seurat-data@d6a8ce6') | |
| ifnb.SeuratData | 3.1.0 | SeuratData::InstallData("ifnb") | If the code exits out, set SeuratData::options(timeout = 240). |
| EnhancedVolcano | 1.18.0 | BiocManager::install() |


## Copyright

This code is copyright by Harrison Wang in 2023.
